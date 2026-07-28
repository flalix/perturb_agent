#!/usr/bin/env python3
"""
bp_prepare.py — Preparação de inputs para BayesPrism (PAAD bulk deconvolution).

Gera, em `outdir`:
    bulk.tsv                  amostras x genes (contagens brutas, inteiros)
    sc_matrix.mtx             células x genes, esparsa (MatrixMarket)
    sc_genes.txt              símbolos dos genes (colunas de sc_matrix)
    sc_cells.txt              barcodes (linhas de sc_matrix)
    cell_type_labels.txt      rótulo grosseiro por célula  (~8-12 níveis)
    cell_state_labels.txt     rótulo fino por célula       (estados/clusters)

Requer: scanpy, anndata, pandas, numpy, scipy
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

log = logging.getLogger("bp_prepare")


# --------------------------------------------------------------------------- #
# 1. Bulk
# --------------------------------------------------------------------------- #
def load_bulk(
    counts_path: str | Path,
    id2symbol: dict[str, str] | None = None,
    min_samples_expressed: int = 10,
) -> pd.DataFrame:
    """
    Lê a matriz de contagens brutas do bulk e devolve amostras x genes.

    IMPORTANTE: BayesPrism exige CONTAGENS BRUTAS. Não passe TPM, FPKM,
    log-CPM, nem VST/rlog. O modelo é multinomial/Poisson sobre contagens —
    dados transformados quebram a inferência silenciosamente (theta parece
    razoável, o Z sai errado).

    Parameters
    ----------
    counts_path
        TSV/CSV com genes nas linhas e amostras nas colunas (layout típico
        do TCGA `*.htseq_counts.tsv` / STAR-counts do GDC).
    id2symbol
        Mapa Ensembl -> símbolo HGNC. Se None, mantém os IDs como estão.
        BayesPrism alinha genes automaticamente entre referência e mistura,
        mas os NAMESPACES precisam bater (Ensembl vs Ensembl, ou símbolo
        vs símbolo). Misturar os dois zera a interseção.
    """
    log.info("Lendo bulk de %s", counts_path)
    bulk = pd.read_csv(counts_path, sep=None, engine="python", index_col=0)

    # Remove sufixo de versão do Ensembl (ENSG00000141510.16 -> ENSG00000141510).
    bulk.index = bulk.index.astype(str).str.replace(r"\.\d+$", "", regex=True)

    # Descarta as linhas de sumário do HTSeq, se presentes.
    bulk = bulk[~bulk.index.str.startswith("__")]

    if id2symbol is not None:
        symbols = bulk.index.map(id2symbol)
        bulk = bulk.loc[symbols.notna()]
        bulk.index = symbols[symbols.notna()]
        # Genes duplicados após o mapeamento: soma as contagens (não média —
        # continuamos no espaço de contagens).
        bulk = bulk.groupby(level=0).sum()

    # Arredonda: alguns pipelines (RSEM, Salmon) devolvem contagens fracionárias.
    if not np.issubdtype(bulk.values.dtype, np.integer):
        log.warning("Contagens não inteiras detectadas — arredondando.")
        bulk = bulk.round().astype(np.int64)

    keep = (bulk > 0).sum(axis=1) >= min_samples_expressed
    log.info("Bulk: %d genes -> %d após filtro de expressão", len(keep), keep.sum())
    bulk = bulk.loc[keep]

    return bulk.T  # amostras x genes


# --------------------------------------------------------------------------- #
# 2. Referência scRNA-seq
# --------------------------------------------------------------------------- #
def build_reference(
    h5ad_path: str | Path,
    type_key: str = "cell_type",
    patient_key: str = "patient",
    malignant_label: str = "malignant",
    n_malignant_clusters: int = 4,
    max_cells_per_state: int = 300,
    random_state: int = 0,
) -> ad.AnnData:
    """
    Carrega a referência de célula única e constrói `cell.type` / `cell.state`.

    Referências publicadas de PDAC (escolha uma, ou combine):
        Peng et al. 2019, Cell Res     — CRA001160  (24 PDAC, ~41k células)
        Steele et al. 2020, Nat Cancer — GSE155698  (16 PDAC + 3 adjacentes)
        Raghavan et al. 2021, Cell     — inclui organoides/metástase
        Lin et al. 2020                — GSE154778  (metastático)

    O ponto mais importante desta função
    ------------------------------------
    `cell.state.labels` das células malignas DEVEM ser subclusters
    DENTRO DE CADA PACIENTE. Se você rotular todas as células tumorais
    como um único estado "malignant", o prior do compartimento tumoral vira
    um único perfil médio e o posterior Z_tumor é encolhido em direção a ele —
    exatamente o eixo classical/basal que você quer recuperar é o que se perde.
    Com subclusters por paciente, o prior cobre o espectro classical<->basal e
    o Gibbs consegue mover cada amostra dentro dele.

    `cell.type.labels`, em contraste, deve ser GROSSEIRO. É nele que o theta
    final é somado. Tipos com perfis quase colineares (p.ex. myCAF vs iCAF)
    ficam identificáveis a nível de state, não de type.
    """
    log.info("Lendo referência de %s", h5ad_path)
    adata = sc.read_h5ad(h5ad_path)

    # BayesPrism quer contagens brutas. Se .raw ou layers['counts'] existir, use.
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    elif adata.raw is not None:
        adata = adata.raw.to_adata()
    _assert_counts(adata.X)

    adata.var_names_make_unique()

    # ---- cell.type: grosseiro ----
    adata.obs["cell_type_bp"] = adata.obs[type_key].astype(str)

    # ---- cell.state: fino ----
    adata.obs["cell_state_bp"] = adata.obs["cell_type_bp"].copy()

    is_mal = adata.obs["cell_type_bp"] == malignant_label
    if is_mal.sum() == 0:
        raise ValueError(
            f"Nenhuma célula com cell_type == {malignant_label!r}. "
            f"Níveis disponíveis: {sorted(adata.obs[type_key].unique())}"
        )

    log.info("Subclusterizando %d células malignas por paciente", is_mal.sum())
    for pat in adata.obs.loc[is_mal, patient_key].unique():
        mask = is_mal & (adata.obs[patient_key] == pat)
        if mask.sum() < 50:  # poucas células: mantém como um estado só
            adata.obs.loc[mask, "cell_state_bp"] = f"{malignant_label}_{pat}"
            continue

        sub = adata[mask].copy()
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sc.pp.highly_variable_genes(sub, n_top_genes=2000)
        sub = sub[:, sub.var.highly_variable].copy()
        sc.pp.scale(sub, max_value=10)
        sc.tl.pca(sub, n_comps=min(30, sub.n_obs - 1, sub.n_vars - 1))
        sc.pp.neighbors(sub, n_neighbors=15)
        sc.tl.leiden(sub, key_added="mal_sub", flavor="igraph", n_iterations=2)

        # Colapsa para no máximo n_malignant_clusters estados por paciente:
        # muitos estados com poucas células deixam o prior ruidoso.
        lab = sub.obs["mal_sub"].astype(str)
        top = lab.value_counts().index[:n_malignant_clusters]
        lab = lab.where(lab.isin(top), other="rest")

        adata.obs.loc[mask, "cell_state_bp"] = [
            f"{malignant_label}_{pat}_{c}" for c in lab
        ]

    # ---- Subamostragem por estado ----
    # BayesPrism v2.2 aceita dgCMatrix esparsa, mas o custo do Gibbs cresce
    # com o número de estados x genes, não de células. Subamostrar mantém
    # a memória sob controle sem perder o prior.
    rng = np.random.default_rng(random_state)
    keep_idx: list[np.ndarray] = []
    for state, idx in pd.Series(np.arange(adata.n_obs)).groupby(
        adata.obs["cell_state_bp"].values
    ):
        arr = idx.to_numpy()
        if len(arr) > max_cells_per_state:
            arr = rng.choice(arr, max_cells_per_state, replace=False)
        keep_idx.append(arr)
    keep = np.sort(np.concatenate(keep_idx))

    log.info("Referência: %d -> %d células, %d estados, %d tipos",
             adata.n_obs, len(keep),
             adata.obs["cell_state_bp"].nunique(),
             adata.obs["cell_type_bp"].nunique())

    adata = adata[keep].copy()

    # Estados com <20 células degradam o prior; BayesPrism recomenda >=20-50.
    small = adata.obs["cell_state_bp"].value_counts()
    small = small[small < 20]
    if len(small):
        log.warning("Estados com <20 células (considere fundir): %s",
                    ", ".join(small.index[:10]))

    return adata


def _assert_counts(X) -> None:
    """Falha cedo se o usuário passou dados transformados."""
    sample = X[:50].toarray() if sp.issparse(X) else np.asarray(X[:50])
    if sample.size and (not np.allclose(sample, np.round(sample)) or sample.min() < 0):
        raise ValueError(
            "A matriz não parece contar contagens brutas (valores não inteiros "
            "ou negativos). BayesPrism precisa de counts, não de log/CPM/scaled."
        )


# --------------------------------------------------------------------------- #
# 3. Escrita
# --------------------------------------------------------------------------- #
def export(bulk: pd.DataFrame, ref: ad.AnnData, outdir: str | Path) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    overlap = len(set(bulk.columns) & set(ref.var_names))
    log.info("Genes em comum bulk<->referência: %d", overlap)
    if overlap < 5000:
        raise ValueError(
            f"Apenas {overlap} genes em comum. Quase certamente é incompatibilidade "
            "de namespace (Ensembl vs símbolo HGNC) ou de build. Resolva antes de rodar."
        )

    bulk.to_csv(outdir / "bulk.tsv", sep="\t")

    X = ref.X if sp.issparse(ref.X) else sp.csr_matrix(ref.X)
    scipy.io.mmwrite(str(outdir / "sc_matrix.mtx"), X.astype(np.int32))

    for name, values in [
        ("sc_genes.txt", ref.var_names),
        ("sc_cells.txt", ref.obs_names),
        ("cell_type_labels.txt", ref.obs["cell_type_bp"]),
        ("cell_state_labels.txt", ref.obs["cell_state_bp"]),
    ]:
        pd.Series(list(values)).to_csv(
            outdir / name, index=False, header=False
        )

    log.info("Inputs escritos em %s", outdir)


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bulk", required=True, help="TSV/CSV de contagens (genes x amostras)")
    ap.add_argument("--sc", required=True, help=".h5ad da referência scRNA-seq")
    ap.add_argument("--outdir", default="bp_input")
    ap.add_argument("--type-key", default="cell_type")
    ap.add_argument("--patient-key", default="patient")
    ap.add_argument("--malignant-label", default="malignant")
    ap.add_argument("--map", default=None,
                    help="CSV com colunas ensembl,symbol para converter o bulk")
    ap.add_argument("--max-cells-per-state", type=int, default=300)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    id2symbol = None
    if args.map:
        m = pd.read_csv(args.map)
        id2symbol = dict(zip(m.ensembl, m.symbol))

    bulk = load_bulk(args.bulk, id2symbol=id2symbol)
    ref = build_reference(
        args.sc,
        type_key=args.type_key,
        patient_key=args.patient_key,
        malignant_label=args.malignant_label,
        max_cells_per_state=args.max_cells_per_state,
    )
    export(bulk, ref, args.outdir)


if __name__ == "__main__":
    main()
