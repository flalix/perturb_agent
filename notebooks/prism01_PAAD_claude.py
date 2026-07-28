#!/usr/bin/env python3
"""
bp_postprocess.py — Pós-processamento do BayesPrism para o problema PAAD.

O objetivo: reclassificar classical vs basal-like usando SOMENTE o
compartimento maligno (Z_malignant), removendo o confundimento com pureza
tumoral que domina o bulk.

Requer: pandas, numpy, scipy, scikit-learn, statsmodels
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, zscore
from sklearn.decomposition import NMF

log = logging.getLogger("bp_post")

# --------------------------------------------------------------------------- #
# Assinaturas de subtipo. Moffitt et al. 2015 Nat Genet (top-25 de cada).
# --------------------------------------------------------------------------- #
MOFFITT_BASAL = [
    "VGLL1", "UCA1", "S100A2", "LY6D", "SPRR3", "SPRR1B", "LEMD1", "KRT15",
    "CTSL2", "DHRS9", "AREG", "CST6", "SERPINB3", "KRT6C", "KRT6A", "SERPINB4",
    "FAM83A", "SCEL", "FGFBP1", "KRT7", "KRT17", "GPR87", "TNS4", "SLC2A1", "ANXA8L2",
]
MOFFITT_CLASSICAL = [
    "BTNL8", "FAM3D", "ATAD4", "AGR3", "CTSE", "LOC400573", "LYZ", "TFF2",
    "TFF1", "ANXA10", "LGALS4", "PLA2G10", "CEACAM6", "VSIG2", "TSPAN8",
    "ST6GALNAC1", "AGR2", "TFF3", "CYP3A7", "MYO1A", "CLRN3", "KRT20",
    "CDH17", "SPINK4", "REG4",
]

# Moffitt tumor-stroma: úteis como CONTROLE NEGATIVO. Depois da deconvolução,
# a assinatura estromal ativada NÃO deve mais separar as amostras dentro do
# Z_malignant. Se ainda separar, sobrou vazamento estromal.
MOFFITT_STROMA_ACTIVATED = [
    "SPARC", "COL10A1", "COL5A2", "COL1A2", "WNT2", "CTHRC1", "INHBA",
    "GJB2", "COL1A1", "SDC1", "FAP", "POSTN", "THBS2", "SULF1", "COMP",
]


# --------------------------------------------------------------------------- #
def load_results(outdir: str | Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Lê frações, CVs e as matrizes Z."""
    outdir = Path(outdir)
    theta = pd.read_csv(outdir / "fractions_type.csv", index_col=0)
    cv = pd.read_csv(outdir / "fractions_type_cv.csv", index_col=0)
    Z = {
        p.stem[2:]: pd.read_csv(p, index_col=0)
        for p in sorted(outdir.glob("Z_*.csv"))
    }
    log.info("theta: %s | compartimentos Z: %s", theta.shape, list(Z))
    return theta, cv, Z


def normalize_Z(Z: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """
    Z vem em escala de contagens e sua soma por amostra é proporcional à
    ABUNDÂNCIA daquele compartimento. Normalizar por linha é obrigatório:
    sem isso você reintroduz a pureza como primeiro componente e desfaz
    todo o propósito da deconvolução.
    """
    cpm = Z.div(Z.sum(axis=1), axis=0) * 1e6
    return np.log2(cpm + pseudocount)


def signature_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Média dos z-scores dos genes da assinatura (amostras x genes)."""
    present = [g for g in genes if g in expr.columns]
    missing = set(genes) - set(present)
    if missing:
        log.warning("%d/%d genes ausentes: %s",
                    len(missing), len(genes), sorted(missing)[:8])
    if len(present) < 5:
        raise ValueError("Assinatura com <5 genes presentes — verifique o namespace.")
    sub = expr[present]
    # z-score entre amostras, por gene
    z = sub.apply(lambda col: zscore(col, nan_policy="omit"), axis=0)
    return z.mean(axis=1)


# --------------------------------------------------------------------------- #
def subtype_from_tumor_compartment(
    Z_mal: pd.DataFrame,
    n_top_genes: int = 2000,
) -> pd.DataFrame:
    """
    Escores classical/basal e chamada de subtipo a partir do Z tumoral.

    Devolve DataFrame com basal, classical, delta (basal - classical),
    subtype e o cluster hierárquico independente.
    """
    expr = normalize_Z(Z_mal)

    basal = signature_score(expr, MOFFITT_BASAL)
    classical = signature_score(expr, MOFFITT_CLASSICAL)
    stroma = signature_score(expr, MOFFITT_STROMA_ACTIVATED)

    out = pd.DataFrame(
        {"basal": basal, "classical": classical, "stroma_leak": stroma}
    )
    out["delta"] = out.basal - out.classical
    out["subtype"] = np.where(out.delta > 0, "basal-like", "classical")

    # Clusterização não-supervisionada independente, para checar se a
    # estrutura emerge sozinha e não só por projeção nas assinaturas.
    hv = expr.var(axis=0).nlargest(n_top_genes).index
    X = expr[hv].apply(lambda c: zscore(c, nan_policy="omit"), axis=0).fillna(0)
    link = linkage(pdist(X.values, metric="correlation"), method="ward")
    out["hclust_k2"] = fcluster(link, t=2, criterion="maxclust")

    # Concordância entre a chamada por assinatura e a clusterização cega.
    tab = pd.crosstab(out.subtype, out.hclust_k2)
    log.info("Concordância assinatura x hclust:\n%s", tab)

    return out


def check_purity_deconfounding(
    theta: pd.DataFrame,
    result: pd.DataFrame,
    malignant_col: str = "malignant",
) -> pd.DataFrame:
    """
    O teste que decide se a deconvolução funcionou.

    No bulk cru, o escore basal correlaciona fortemente com pureza tumoral
    (é o artefato que você quer eliminar). Depois da deconvolução, o delta
    calculado sobre Z_malignant deve ser aproximadamente ORTOGONAL à
    fração maligna. Se |r| ainda for alto (>~0.4), sobrou confundimento:
    tipicamente porque o prior maligno tinha um único estado, ou porque
    tipos celulares reais faltaram na referência e seu sinal foi absorvido
    pelo compartimento tumoral.
    """
    purity = theta[malignant_col].reindex(result.index)
    rows = []
    for col in ["basal", "classical", "delta", "stroma_leak"]:
        ok = purity.notna() & result[col].notna()
        r, p = pearsonr(purity[ok], result.loc[ok, col])
        rows.append({"score": col, "r_com_pureza": round(r, 3),
                     "p": f"{p:.2e}",
                     "veredito": "OK" if abs(r) < 0.4 else "CONFUNDIDO"})
    df = pd.DataFrame(rows)
    log.info("Deconfundimento por pureza:\n%s", df.to_string(index=False))
    return df


def malignant_programs(Z_mal: pd.DataFrame, n_programs: int = 5) -> tuple:
    """
    NMF sobre o compartimento maligno — o análogo pós-deconvolução da sua
    clusterização original, mas agora sem TME. Programas soft-assigned
    lidam melhor do que clusters rígidos com o fato de que classical e
    basal coexistem na mesma amostra (Chan-Seng-Yue et al. 2020).
    """
    expr = normalize_Z(Z_mal)
    hv = expr.var(axis=0).nlargest(3000).index
    X = expr[hv]
    X = (X - X.min().min())  # NMF exige não-negatividade

    model = NMF(n_components=n_programs, init="nndsvda",
                max_iter=2000, random_state=0)
    W = pd.DataFrame(model.fit_transform(X.values), index=X.index,
                     columns=[f"program_{i+1}" for i in range(n_programs)])
    H = pd.DataFrame(model.components_, index=W.columns, columns=X.columns)

    W_norm = W.div(W.sum(axis=1), axis=0)  # uso relativo de programa
    top_genes = {c: H.loc[c].nlargest(30).index.tolist() for c in H.index}
    return W_norm, H, top_genes


# --------------------------------------------------------------------------- #
def main(outdir: str = "bp_output", malignant: str = "malignant") -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    theta, cv, Z = load_results(outdir)

    if malignant not in Z:
        raise KeyError(f"Z_{malignant} não encontrado. Disponíveis: {list(Z)}")

    result = subtype_from_tumor_compartment(Z[malignant])
    qc = check_purity_deconfounding(theta, result, malignant_col=malignant)
    W, H, top = malignant_programs(Z[malignant])

    out = Path(outdir)
    result.join(theta.add_prefix("frac_")).to_csv(out / "subtype_calls.csv")
    qc.to_csv(out / "qc_purity_deconfounding.csv", index=False)
    W.to_csv(out / "malignant_program_usage.csv")
    pd.DataFrame({k: pd.Series(v) for k, v in top.items()}).to_csv(
        out / "malignant_program_top_genes.csv", index=False
    )

    # Frações com alta incerteza — não use em testes downstream.
    unstable = (cv > 0.5).sum().sort_values(ascending=False)
    log.info("Amostras com CV(theta) > 0.5 por tipo:\n%s", unstable.to_string())

    log.info("Escrito em %s", out)


if __name__ == "__main__":
    import fire  # pip install fire; ou troque por argparse
    fire.Fire(main)