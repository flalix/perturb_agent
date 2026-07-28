#!/usr/bin/env python3
"""
batch_diag.py — Diagnóstico e tratamento do efeito TCGA vs CPTAC3
                antes da deconvolução BayesPrism.

Três funções, na ordem em que devem ser usadas:

    check_confounding()      decide SE você pode corrigir
    antisense_blacklist()    corrige a causa (strandedness) em vez do sintoma
    validate_theta()         verifica, depois, se sobrou efeito de estudo

Requer: pandas, numpy, scipy, statsmodels, pyranges (só para a blacklist)
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2_contingency

log = logging.getLogger("batch_diag")

# Nomes equivalentes aceitos, para que o metadado bruto funcione sem conversão.
_STUDY_ALIASES = ("study", "program", "cohort", "project", "batch", "site")
_GROUP_ALIASES = ("group", "condition", "sample_type", "tissue", "status")
_SAMPLE_ALIASES = ("sample", "sample_id", "cols", "barcode", "case_id")
_CONTROL_SYNONYMS = {"normal", "adjacent", "nat", "healthy",
                     "non_tumor", "non-tumor", "benign"}


def _load(obj, index_alias: tuple[str, ...] = _SAMPLE_ALIASES) -> pd.DataFrame:
    """
    Aceita um DataFrame já carregado OU um caminho para csv/tsv.

    Era esta a peça que faltava: chamado pela linha de comando, o fire
    entrega a string do caminho, não o objeto.
    """
    print("loading ...")
    
    if isinstance(obj, pd.DataFrame):
        return obj
    path = Path(obj)
    if not path.exists():
        raise FileNotFoundError(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    df = pd.read_csv(path, sep=sep)
    first = str(df.columns[0]).lower()
    for c in df.columns:
        if str(c).lower() in index_alias:
            return df.set_index(c)
    if first.startswith("unnamed"):
        return df.set_index(df.columns[0])
    return df


def _resolve(df: pd.DataFrame, requested: str,
             aliases: tuple[str, ...], what: str) -> str:
    """Encontra a coluna pedida, ou a primeira equivalente conhecida."""
    if requested in df.columns:
        return requested
    lower = {str(c).lower(): c for c in df.columns}
    for a in aliases:
        if a in lower:
            log.info("Coluna de %s: usando %r", what, lower[a])
            return lower[a]
    raise KeyError(
        f"Nenhuma coluna de {what} encontrada. Procurei por {requested!r} e "
        f"por {aliases}. Colunas disponíveis: {list(df.columns)}"
    )


# =========================================================================== #
# 1. Confundimento estudo x grupo
# =========================================================================== #
def check_confounding(design,
                      study_col: str = "study",
                      group_col: str = "group") -> dict:
    """
    Decide se a correção de lote é sequer admissível.

    O ponto: ComBat_seq remove a variação atribuída ao lote. Se lote e
    grupo biológico forem parcialmente colineares, parte do efeito
    biológico é atribuída ao lote e removida junto — e o resultado parece
    limpo, o que é pior do que parecer errado.

    Em PAAD isso é um risco concreto: o TCGA-PAAD tem pouquíssimos normais
    de pâncreas, enquanto o CPTAC3 tem dezenas. Se os seus 21 controles
    forem majoritariamente de um único estudo, estudo ~ grupo.

    Veredito
    --------
    OK            cada estudo contribui com ambos os grupos em número
                  razoável -> ComBat_seq(group=) é seguro.
    PARCIAL       desbalanço forte -> corrigir encolhe o efeito real;
                  prefira modelar o estudo como covariável.
    CONFUNDIDO    algum estudo tem só um grupo -> NÃO corrija. O efeito de
                  estudo e o efeito de grupo são inseparáveis nos dados.
    """
    design = _load(design)
    study_col = _resolve(design, study_col, _STUDY_ALIASES, "estudo")
    group_col = _resolve(design, group_col, _GROUP_ALIASES, "grupo")

    grp = design[group_col].astype(str).str.strip().str.lower()
    grp = grp.replace({s: "control" for s in _CONTROL_SYNONYMS})

    tab = pd.crosstab(design[study_col], grp)
    log.info("Tabela cruzada estudo x grupo:\n%s", tab.to_string())

    zero_cells = int((tab == 0).sum().sum())
    min_cell = int(tab.values.min())

    # Cramér's V como medida de associação entre as duas variáveis.
    chi2 = chi2_contingency(tab.values + 0.5)[0]
    n = tab.values.sum()
    cramers_v = float(np.sqrt(chi2 / (n * (min(tab.shape) - 1))))

    if zero_cells > 0:
        verdict = "CONFUNDIDO"
        advice = ("Não use ComBat_seq. Modele o estudo como covariável no DE "
                  "e/ou restrinja o contraste ao estudo que tem os dois grupos.")
    elif min_cell < 5 or cramers_v > 0.7:
        verdict = "PARCIAL"
        advice = ("ComBat_seq vai encolher o efeito biológico. Prefira não "
                  "corrigir e usar ~ study + group no modelo.")
    else:
        verdict = "OK"
        advice = "ComBat_seq(group=) é seguro, se você ainda quiser corrigir."

    log.info("Cramér's V = %.3f | célula mínima = %d | veredito: %s",
             cramers_v, min_cell, verdict)
    log.info("%s", advice)

    return {"crosstab": tab, "cramers_v": cramers_v,
            "min_cell": min_cell, "verdict": verdict, "advice": advice}


# =========================================================================== #
# 2. Blacklist de genes com sobreposição antisenso
# =========================================================================== #
def antisense_blacklist(gtf_path: str | Path,
                        min_overlap_bp: int = 1,
                        feature: str = "exon") -> pd.Series:
    """
    Lista de genes cuja quantificação é ambígua em biblioteca não-direcional.

    Um gene entra na lista se algum éxon seu se sobrepõe a um éxon de outro
    gene na FITA OPOSTA. Nessas regiões o TCGA (não-stranded) não distingue
    a leitura do transcrito sense da do antisenso, e a contagem de ambos
    fica errada; o CPTAC3 (stranded) resolve. Remover esses genes elimina o
    viés na origem em vez de tentar ajustá-lo depois.

    Use o MESMO GTF que gerou as suas contagens (GENCODE v36 para o GDC
    atual). Anotações diferentes dão listas diferentes.
    """
    import pyranges as pr

    log.info("Lendo GTF: %s", gtf_path)
    gr = pr.read_gtf(str(gtf_path))
    ex = gr[gr.Feature == feature]
    ex = ex[["gene_id", "gene_name", "Strand"]]

    plus = ex[ex.Strand == "+"]
    minus = ex[ex.Strand == "-"]

    # join ignorando fita: sobreposições +/- são exatamente as ambíguas
    hits = plus.join(minus, strandedness=False, report_overlap=True)
    hits = hits.df
    hits = hits[hits.Overlap >= min_overlap_bp]

    flagged = pd.unique(np.concatenate([
        hits.gene_id.values, hits.gene_id_b.values
    ]))
    flagged = pd.Series(flagged, name="gene_id").str.replace(
        r"\.\d+$", "", regex=True
    ).drop_duplicates()

    total = ex.df.gene_id.nunique()
    log.info("Genes com sobreposição antisenso exônica: %d / %d (%.1f%%)",
             len(flagged), total, 100 * len(flagged) / total)
    return flagged


def apply_blacklist(counts: pd.DataFrame,
                    blacklist: pd.Series) -> pd.DataFrame:
    """counts: genes x amostras, índice = Ensembl sem versão."""
    before = counts.shape[0]
    out = counts.loc[~counts.index.isin(set(blacklist))]
    log.info("Blacklist aplicada: %d -> %d genes", before, out.shape[0])
    return out


# =========================================================================== #
# 3. Validação pós-deconvolução
# =========================================================================== #
def validate_theta(theta,
                   design,
                   study_col: str = "study",
                   group_col: str = "group",
                   n_perm: int = 999,
                   random_state: int = 0) -> pd.DataFrame:
    """
    O teste que realmente decide se a sua escolha funcionou.

    Restringe aos TUMORES (para tirar o grupo da equação) e pergunta se a
    composição celular estimada ainda separa por estudo. Biologicamente,
    duas coortes de PDAC não têm por que diferir sistematicamente em
    composição do microambiente; se θ separar, o que sobrou é técnico.

    Roda duas coisas:
      - PERMANOVA (pseudo-F por permutação) sobre a matriz de distâncias
        de Aitchison entre perfis composicionais;
      - teste por compartimento, para localizar ONDE está a diferença.

    Leitura: R² < ~0.05 e p > 0.05 = θ comparável entre estudos, pode
    prosseguir. R² alto = o efeito de estudo entrou na composição e
    qualquer LFC por tipo celular estará contaminado.
    """
    rng = np.random.default_rng(random_state)

    theta = _load(theta)
    design = _load(design)
    study_col = _resolve(design, study_col, _STUDY_ALIASES, "estudo")
    group_col = _resolve(design, group_col, _GROUP_ALIASES, "grupo")

    idx = design.index[
        design[group_col].astype(str).str.strip().str.lower() == "tumor"
    ]
    idx = idx.intersection(theta.index)
    th = theta.loc[idx]
    grp = design.loc[idx, study_col].astype(str)

    if grp.nunique() < 2:
        raise ValueError("Menos de dois estudos entre os tumores.")

    # Distância de Aitchison: θ é composicional (soma 1), distância
    # euclidiana direta sobre proporções é inadequada.
    clr = np.log(th.clip(lower=1e-6))
    clr = clr.sub(clr.mean(axis=1), axis=0)
    D = squareform(pdist(clr.values, metric="euclidean"))

    def pseudo_f(labels: np.ndarray) -> float:
        n = len(labels)
        ss_total = (D ** 2).sum() / (2 * n)
        ss_within = 0.0
        for lv in np.unique(labels):
            m = labels == lv
            k = m.sum()
            if k > 1:
                ss_within += (D[np.ix_(m, m)] ** 2).sum() / (2 * k)
        a = len(np.unique(labels))
        return ((ss_total - ss_within) / (a - 1)) / (ss_within / (n - a))

    obs_labels = grp.values
    f_obs = pseudo_f(obs_labels)
    null = np.array([pseudo_f(rng.permutation(obs_labels))
                     for _ in range(n_perm)])
    p_perm = (np.sum(null >= f_obs) + 1) / (n_perm + 1)

    ss_total = (D ** 2).sum() / (2 * len(obs_labels))
    ss_within = sum(
        (D[np.ix_(obs_labels == lv, obs_labels == lv)] ** 2).sum()
        / (2 * (obs_labels == lv).sum())
        for lv in np.unique(obs_labels) if (obs_labels == lv).sum() > 1
    )
    r2 = (ss_total - ss_within) / ss_total

    log.info("PERMANOVA (só tumores): pseudo-F = %.2f, R² = %.3f, p = %.4f",
             f_obs, r2, p_perm)
    log.info("Veredito: %s",
             "θ comparável entre estudos" if (r2 < 0.05 and p_perm > 0.05)
             else "EFEITO DE ESTUDO RESIDUAL em θ")

    # Por compartimento: Mann-Whitney nos dois maiores estudos.
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    top2 = grp.value_counts().index[:2]
    a, b = th[grp == top2[0]], th[grp == top2[1]]
    rows = []
    for ct in th.columns:
        u, p = mannwhitneyu(a[ct], b[ct], alternative="two-sided")
        rows.append({"cell_type": ct,
                     f"median_{top2[0]}": round(a[ct].median(), 4),
                     f"median_{top2[1]}": round(b[ct].median(), 4),
                     "p": p})
    res = pd.DataFrame(rows)
    res["padj"] = multipletests(res.p, method="fdr_bh")[1]
    res = res.sort_values("padj")
    log.info("Diferença de θ por compartimento entre estudos:\n%s",
             res.to_string(index=False))

    res.attrs.update({"permanova_F": f_obs, "permanova_R2": r2,
                      "permanova_p": p_perm})
    return res


# =========================================================================== #
def _cli() -> None:
    """CLI sem dependencias extras (argparse em vez de fire)."""
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("check_confounding")
    a.add_argument("--design", required=True)
    a.add_argument("--study-col", default="study")
    a.add_argument("--group-col", default="group")

    b = sub.add_parser("antisense_blacklist")
    b.add_argument("--gtf", required=True)
    b.add_argument("--out", default="antisense_blacklist.txt")

    c = sub.add_parser("validate_theta")
    c.add_argument("--theta", required=True)
    c.add_argument("--design", required=True)
    c.add_argument("--out", default="theta_by_study.csv")

    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.cmd == "check_confounding":
        check_confounding(args.design, args.study_col, args.group_col)
    elif args.cmd == "antisense_blacklist":
        bl = antisense_blacklist(args.gtf)
        bl.to_csv(args.out, index=False, header=False)
        log.info("Blacklist escrita em %s", args.out)
    elif args.cmd == "validate_theta":
        validate_theta(args.theta, args.design).to_csv(args.out, index=False)
        log.info("Escrito em %s", args.out)


if __name__ == "__main__":
    _cli()