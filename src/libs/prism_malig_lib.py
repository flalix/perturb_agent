"""
prism/malignant_state.py
========================

Sample-level clustering of the *malignant compartment* recovered by BayesPrism,
with perturbation anchoring against Tahoe-100M and (optionally) scGPT.

Pipeline
--------
    export_bayesprism.R                     # R -> flat TSVs
        |
    load_bayesprism_export()                # -> MalignantState
        |
    prepare_malignant_matrix()              # compartment CPM, share filter,
        |                                   # purity decoupling, HVG
    consensus_cluster()                     # k selection via PAC / silhouette
        |
    cluster_signatures()                    # one-vs-rest signed signatures
        |
    +-- score_clusters_vs_tahoe()           # CMap-style WTCS / NCS
    +-- scgpt_embed_pseudobulk()            # foundation-model embedding space

Design notes that matter more than the code
-------------------------------------------
1.  `get.exp()` returns posterior-mean counts *attributable to* the malignant
    compartment. Its row sums are ~ theta_mal * library_size. Clustering the
    raw matrix therefore recovers tumour purity, not malignant state. We CPM
    within the compartment before anything else.

2.  BayesPrism shrinks Z_gs toward the scRNA-seq reference for genes with
    little compartment-specific evidence. For those genes the between-sample
    variance is a projection of the reference, not of the sample. We drop them
    with a malignant-share + effective-count filter (`min_share`, `min_counts`).
    Skipping this is the single easiest way to manufacture clusters.

3.  Purity is still a latent axis after CPM (higher theta_mal => less shrinkage
    => more sample-specific Z). `decouple_purity=True` regresses each gene on
    theta_mal and clusters the residuals. Always run both and compare; if the
    structure only exists in the un-decoupled version, it is a purity axis.

4.  Tahoe-100M is 47 cancer cell lines (pancreas is one of the better-covered
    organs), 379 compounds, 24 h, Parse 3'. A deconvolved malignant compartment
    is genuinely the closest bulk object to a cell line -- no TME, no stroma --
    which is what makes this comparison defensible at all. It still says
    nothing about stromal/immune programs.
"""

from __future__ import annotations

import glob
import os
import re
import warnings
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from pathlib import Path

from scipy.cluster.hierarchy import fcluster, linkage, cophenet
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from libs.Basic import pdreadcsv, pdwritecsv, create_dir


__all__ = [
    "MalignantState",
    "load_bayesprism_export",
    "prepare_malignant_matrix",
    "consensus_cluster",
    "choose_k",
    "cluster_signatures",
    "load_tahoe_de",
    "score_clusters_vs_tahoe",
    "scgpt_embed_pseudobulk",
]


# ===========================================================================
# 0. Container
# ===========================================================================
@dataclass
class MalignantState:
    """Everything pulled out of a BayesPrism `res` for one compartment."""

    Z: pd.DataFrame                       # samples x genes, compartment counts
    df_theta: pd.DataFrame                   # samples x cell types
    cell_name: str = "malignant"
    Z_all: Optional[Dict[str, pd.DataFrame]] = None   # ct -> samples x genes
    share: Optional[pd.DataFrame] = None  # samples x genes, malignant fraction
    meta: Dict = field(default_factory=dict)

    @property
    def theta_mal(self) -> pd.Series:
        return self.df_theta[self.cell_name]

    def compute_share(self) -> pd.DataFrame:
        """Per-sample, per-gene fraction of total expression assigned to the
        malignant compartment. Requires Z_all."""
        if self.Z_all is None:
            raise ValueError(
                "Z_all is required. Re-run export_bayesprism.R so that "
                "Z_all/ is populated, or supply `share=` directly (e.g. from "
                "PRISM.gene_compartment_share_nnls)."
            )
        genes = self.Z.columns
        total = np.zeros(self.Z.shape, dtype=float)
        for ct, Zct in self.Z_all.items():
            total += Zct.reindex(index=self.Z.index, columns=genes).values
        with np.errstate(invalid="ignore", divide="ignore"):
            sh = np.where(total > 0, self.Z.values / total, np.nan)
        self.share = pd.DataFrame(sh, index=self.Z.index, columns=genes)
        return self.share



@dataclass
class MalignantCluster:
    """Everything pulled out of a BayesPrism `res` for one compartment."""

    def __init__(self, prism, res, df_bulk, ref, root_mprog_cluster: Path, cell_name: str = "malignant"):
        self.root_mprog_cluster = root_mprog_cluster

        self.cell_name = cell_name

        self.df_theta = res.df_theta

        self.Z_all, self.genes_all = prism.full_Z(res, df_bulk, ref)
        self.Z = self.Z_all

        self.fname_theta = "df_theta.tsv"
        self.fname_z = f"Z_{self.cell_name}.tsv"
        self.dir_zall = 'z_all'

        self.program1_panel = ("FAM83A-AS1", "TFAP2A-AS2", "HOXA10-AS",
                               "HOXB-AS3", "HOXB-AS4", "LEMD1-AS1", "MIR7-3HG")
    

        self._TAHOE_REPO = "tahoebio/Tahoe-100M"


    # ===========================================================================
    # 1. run --> end-to-end convenience
    # ===========================================================================
    def run(
        self, 
        ks: Iterable[int] = range(2, 9),
        tahoe: bool = False,
        **kw,
    ):
        ms = self.load_bayesprism_export()
        X, diag = self.prepare_malignant_matrix(ms, keep_genes=self.program1_panel, **kw)
        cc = self.consensus_cluster(X, ks=ks)
        k = self.choose_k(cc)
        labels = cc[k]["labels"]
        sig = self.cluster_signatures(X, labels)

        res = {"X": X, "diag": diag, "consensus": cc, "k": k,
            "labels": labels, "signatures": sig}

        if tahoe:
            S, cond = self.load_tahoe_de(organs=("Pancreas",))
            res["tahoe"] = self.score_clusters_vs_tahoe(sig, S, cond)
        return res


    def load_bayesprism_export(
        self,
        load_all: bool = True,
    ) -> MalignantState:
        """Read the TSVs written by export_bayesprism.R."""
        
        df_theta = pdreadcsv(self.fname_theta, self.root_mprog_cluster, sep="\t", index_col=0)

        zpath = self.root_mprog_cluster / self.fname_z
        if not os.path.exists(zpath):
            cands = glob.glob(os.path.join(self.root_mprog_cluster, "Z_*.tsv"))
            raise FileNotFoundError(f"{zpath} missing; found {cands}")
        Z = pdreadcsv(self.fname_z, self.root_mprog_cluster, index_col=0)

        Z_all = None

        zdir = self.root_mprog_cluster / self.dir_zall
        if load_all and os.path.isdir(zdir):
            Z_all = {}
            for f in sorted(glob.glob(os.path.join(zdir, "*.tsv"))):
                ct = os.path.splitext(os.path.basename(f))[0]
                Z_all[ct] = pdreadcsv(f, zdir, index_col=0)

        self.Z = Z
        self.Z_all = Z_all

        if self.cell_name not in df_theta.columns:
            raise KeyError(
                f"'{self.cell_name}' not in df_theta columns {list(df_theta.columns)}. "
                "With a Peng-2019 (CRA001160) reference the malignant label is "
                "often 'Ductal cell type 2'."
            )

        ms = MalignantState(Z=Z, df_theta=df_theta, cell_name=self.cell_name, Z_all=Z_all)
        if Z_all is not None:
            ms.compute_share()
        return ms


    def from_full_Z(
        Z_full,
        genes_full: Sequence[str],
        cell_types: Sequence[str],
        cell_name: str = "Ductal cell type 2",
        samples: Optional[Sequence[str]] = None,
        df_theta: Optional[pd.DataFrame] = None,
        Z_mal: Optional[pd.DataFrame] = None,
        recovered_genes: Optional[Sequence[str]] = None,
    ) -> MalignantState:
        """Build a MalignantState straight from `prism.full_Z(res, df_bulk, ref)`.

        Preferred over `load_bayesprism_export` -- no R round-trip, and the
        per-gene compartment share is exact rather than reconstructed.

            Z_full, genes_full = prism.full_Z(res, df_bulk, ref)
            Zmal = prism.state_expression(Z_full, genes_full, res, "Ductal cell type 2")
            ms = from_full_Z(Z_full, genes_full,
                            cell_types=list(s2t.values()),      # or ref cell types
                            cell_name="Ductal cell type 2",
                            samples=df_bulk.index, Z_mal=Zmal)

        Axis order of `Z_full` is sniffed from the lengths of `genes_full`,
        `cell_types` and `samples`; pass `samples` whenever two axes could be
        confused.

        Parameters
        ----------
        df_theta
            If omitted, fractions are derived from Z_full by marginalising over
            genes. That is the *Z-implied* (df_theta.cs / first-pass) fraction, which
            is not identical to BayesPrism's updated `theta_f`. Pass the final
            thetas explicitly if you want the decoupling covariate to match what
            you report elsewhere.
        recovered_genes
            Loci that were absent from the marker-based fit and reconstructed
            post-hoc (your HOX-antisense / MIR7-3HG panel). These are recorded in
            `meta` and treated separately downstream -- see the note below.
        """
        A = np.asarray(Z_full, dtype=float)
        if A.ndim != 3:
            raise ValueError(
                f"expected a 3-D (sample, gene, cell-type) array from full_Z, "
                f"got shape {A.shape}. If your full_Z returns a dict of 2-D "
                f"matrices keyed by cell type, pass it as MalignantState(Z_all=...)"
            )

        genes, cts = list(genes_full), list(cell_types)
        ng, nct = len(genes), len(cts)
        ns = len(samples) if samples is not None else None

        valid = []
        for ga in range(3):
            for ca in range(3):
                if ga == ca or A.shape[ga] != ng or A.shape[ca] != nct:
                    continue
                sa = ({0, 1, 2} - {ga, ca}).pop()
                if ns is not None and A.shape[sa] != ns:
                    continue
                valid.append((sa, ga, ca))
        if len(valid) != 1:
            raise ValueError(
                f"cannot resolve axes of shape {A.shape} against "
                f"n_genes={ng}, n_cell_types={nct}, n_samples={ns}. "
                "Pass `samples=` to disambiguate."
            )
        sa, ga, ca = valid[0]
        A = np.moveaxis(A, (sa, ga, ca), (0, 1, 2))          # samples x genes x ct

        idx = pd.Index(samples) if samples is not None else pd.RangeIndex(A.shape[0])
        if cell_name not in cts:
            raise KeyError(f"'{cell_name}' not in cell_types {cts}")
        ci = cts.index(cell_name)

        Zc = pd.DataFrame(A[:, :, ci], index=idx, columns=genes)
        total = A.sum(axis=2)
        with np.errstate(invalid="ignore", divide="ignore"):
            share = pd.DataFrame(np.where(total > 0, A[:, :, ci] / total, np.nan),
                                index=idx, columns=genes)

        if df_theta is None:
            m = A.sum(axis=1)
            df_theta = pd.DataFrame(m / m.sum(axis=1, keepdims=True),
                                index=idx, columns=cts)
            theta_source = "derived from Z_full (Z-implied, not theta_f)"
        else:
            df_theta = df_theta.reindex(index=idx)
            if cell_name not in df_theta.columns:
                raise KeyError(f"'{cell_name}' not in supplied df_theta columns")
            theta_source = "user-supplied"

        meta = {"theta_source": theta_source,
                "axes": {"sample": sa, "gene": ga, "cell_type": ca},
                "recovered_genes": list(recovered_genes or [])}

        # Cross-check state_expression output against the Z_full slice.
        if Z_mal is not None:
            Zm = pd.DataFrame(Z_mal).reindex(index=idx)
            common = [g for g in Zm.columns if g in Zc.columns]
            if common:
                r = np.corrcoef(Zm[common].values.ravel(),
                                Zc[common].values.ravel())[0, 1]
                meta["state_expression_concordance"] = float(r)
                if r < 0.99:
                    warnings.warn(
                        f"state_expression('{cell_name}') and the Z_full slice "
                        f"agree at r={r:.3f}. If state_expression applies its own "
                        "normalisation this is expected; otherwise check that the "
                        "cell-type label maps to the same column."
                    )
            extra = [g for g in Zm.columns if g not in Zc.columns]
            if extra:
                meta["genes_only_in_state_expression"] = extra
            Zc = Zm.reindex(columns=Zm.columns)
            share = share.reindex(columns=Zc.columns)

        return MalignantState(Z=Zc, df_theta=df_theta, cell_name=cell_name,
                            share=share, meta=meta)


    # ===========================================================================
    # 1. Matrix preparation
    # ===========================================================================
    def prepare_malignant_matrix(
        self,
        ms: MalignantState,
        min_share: float = 0.50,
        share_frac_samples: float = 0.50,
        min_counts: float = 10.0,
        min_theta: float = 0.05,
        n_hvg: int = 2000,
        decouple_purity: bool = True,
        extra_covariates: Optional[pd.DataFrame] = None,
        keep_genes: Optional[Sequence[str]] = None,
        return_diagnostics: bool = True,
    ) -> Tuple[pd.DataFrame, Dict]:
        """
        Turn compartment counts into a matrix it is safe to cluster.

        Parameters
        ----------
        min_share, share_frac_samples
            Keep a gene only if >= `share_frac_samples` of samples assign at least
            `min_share` of its expression to the malignant compartment. This is the
            shrinkage guard described in the module docstring.
        min_counts
            Minimum median compartment-attributable count.
        min_theta
            Drop samples whose malignant fraction is below this -- their Z is
            almost entirely prior.
        decouple_purity
            Regress log-CPM on theta_mal (+ extra_covariates) per gene and keep the
            residuals.
        keep_genes
            Force-retain genes regardless of filters (e.g. your Program-1 panel:
            FAM83A-AS1, HOXA10-AS, HOXB-AS3, HOXB-AS4, MIR7-3HG, TFAP2A-AS2,
            LEMD1-AS1). They are still reported in diagnostics so you can see
            whether they *would* have survived.
        """
        diag: Dict = {}

        # --- sample filter ------------------------------------------------------
        keep_s = ms.theta_mal >= min_theta
        if (~keep_s).any():
            warnings.warn(
                f"dropping {int((~keep_s).sum())} samples with theta_{ms.cell_name}"
                f" < {min_theta}: {list(ms.theta_mal.index[~keep_s])}"
            )
        Z = ms.Z.loc[keep_s]
        theta_mal = ms.theta_mal.loc[keep_s]
        diag["samples_dropped"] = list(ms.theta_mal.index[~keep_s])

        # --- compartment CPM ----------------------------------------------------
        lib = Z.sum(axis=1)
        cpm = Z.div(lib, axis=0) * 1e6
        logx = np.log2(cpm + 1.0)

        # --- gene filters -------------------------------------------------------
        expressed = (Z.median(axis=0) >= min_counts)
        diag["n_genes_expressed"] = int(expressed.sum())

        if ms.share is not None:
            sh = ms.share.loc[Z.index, Z.columns]
            share_ok = (sh >= min_share).mean(axis=0) >= share_frac_samples
        else:
            warnings.warn(
                "No per-gene malignant share available -- shrinkage guard is OFF. "
                "Interpret clusters with corresponding caution."
            )
            share_ok = pd.Series(True, index=Z.columns)
        diag["n_genes_share_ok"] = int(share_ok.sum())

        keep_g = expressed & share_ok
        forced = []
        if keep_genes:
            forced = [g for g in keep_genes if g in logx.columns]
            diag["forced_genes_passing_filters"] = {
                g: bool(keep_g.get(g, False)) for g in forced
            }
            keep_g.loc[forced] = True
        logx = logx.loc[:, keep_g]
        diag["n_genes_kept"] = int(logx.shape[1])

        # --- purity decoupling --------------------------------------------------
        if decouple_purity:
            C = pd.DataFrame({"df_theta": theta_mal, "log_lib": np.log(lib)})
            if extra_covariates is not None:
                C = C.join(extra_covariates.loc[C.index], how="left")
            C = C.fillna(C.mean())
            D = np.column_stack([np.ones(len(C)), C.values.astype(float)])
            beta, *_ = np.linalg.lstsq(D, logx.values, rcond=None)
            resid = logx.values - D @ beta
            # keep the grand mean so downstream LFCs stay interpretable
            Xc = pd.DataFrame(resid + logx.values.mean(axis=0),
                            index=logx.index, columns=logx.columns)
            diag["covariates"] = list(C.columns)
        else:
            Xc = logx

        # --- HVG ----------------------------------------------------------------
        if n_hvg and n_hvg < Xc.shape[1]:
            v = Xc.var(axis=0)
            forced_set = set(forced)
            order = v.sort_values(ascending=False).index
            sel = [g for g in order if g not in forced_set][: max(0, n_hvg - len(forced_set))]
            sel = list(forced_set) + sel
            Xc = Xc.loc[:, [g for g in Xc.columns if g in set(sel)]]
        diag["n_hvg"] = int(Xc.shape[1])

        # --- purity leakage check ----------------------------------------------
        if Xc.shape[1] >= 2:
            pcs = PCA(n_components=min(5, min(Xc.shape) - 1)).fit_transform(
                (Xc - Xc.mean()).values
            )
            diag["pc_theta_pearson"] = [
                float(np.corrcoef(pcs[:, i], theta_mal.values)[0, 1])
                for i in range(pcs.shape[1])
            ]
        diag["theta_mal"] = theta_mal

        return (Xc, diag) if return_diagnostics else (Xc, {})


    # ===========================================================================
    # 2. Consensus clustering (Monti et al.)
    # ===========================================================================
    def _pac(consensus: np.ndarray, lo: float = 0.1, hi: float = 0.9) -> float:
        iu = np.triu_indices_from(consensus, k=1)
        v = consensus[iu]
        return float(((v > lo) & (v < hi)).mean())


    def consensus_cluster(
        self,
        X: pd.DataFrame,
        ks: Iterable[int] = range(2, 9),
        n_reps: int = 500,
        sample_frac: float = 0.8,
        feature_frac: float = 0.8,
        metric: str = "correlation",
        method: str = "average",
        random_state: int = 0,
    ) -> Dict[int, Dict]:
        """Resampling consensus clustering on samples x features `X`.

        `metric="correlation"` + `method="average"` mirrors the usual
        transcriptomic-subtype convention; use `metric="euclidean", method="ward"`
        if you want a spherical prior.
        """
        rng = np.random.default_rng(random_state)
        n = X.shape[0]
        Xv = X.values
        out: Dict[int, Dict] = {}

        for k in ks:
            M = np.zeros((n, n))       # co-cluster counts
            I = np.zeros((n, n))       # co-sampled counts
            for _ in range(n_reps):
                si = rng.choice(n, size=max(k + 1, int(round(sample_frac * n))),
                                replace=False)
                fi = rng.choice(Xv.shape[1],
                                size=max(2, int(round(feature_frac * Xv.shape[1]))),
                                replace=False)
                sub = Xv[np.ix_(si, fi)]
                d = pdist(sub, metric=metric)
                if not np.all(np.isfinite(d)):
                    continue
                lab = fcluster(linkage(d, method=method), k, criterion="maxclust")
                I[np.ix_(si, si)] += 1
                same = (lab[:, None] == lab[None, :]).astype(float)
                M[np.ix_(si, si)] += same
            with np.errstate(invalid="ignore", divide="ignore"):
                C = np.where(I > 0, M / I, 0.0)
            np.fill_diagonal(C, 1.0)

            d_cons = squareform(1.0 - C, checks=False)
            Zl = linkage(d_cons, method=method)
            labels = fcluster(Zl, k, criterion="maxclust")
            coph, _ = cophenet(Zl, d_cons)

            try:
                sil = float(silhouette_score(1.0 - C, labels, metric="precomputed"))
            except Exception:
                sil = np.nan

            out[k] = {
                "labels": pd.Series(labels, index=X.index, name=f"k{k}"),
                "consensus": pd.DataFrame(C, index=X.index, columns=X.index),
                "pac": _pac(C),
                "cophenetic": float(coph),
                "silhouette": sil,
                "linkage": Zl,
            }
        return out


    def choose_k(self, cc: Dict[int, Dict], pac_tol: float = 0.02) -> int:
        """Smallest k whose PAC is within `pac_tol` of the minimum PAC."""
        pacs = {k: v["pac"] for k, v in cc.items()}
        best = min(pacs.values())
        return min(k for k, p in pacs.items() if p <= best + pac_tol)


    # ===========================================================================
    # 3. Cluster signatures
    # ===========================================================================
    def cluster_signatures(
        self,
        X: pd.DataFrame,
        labels: pd.Series,
        *,
        n_top: int = 150,
        min_abs_lfc: float = 0.0,
    ) -> Dict[int, Dict[str, object]]:
        """One-vs-rest signed signatures with a simple moderated t.

        Returns {cluster: {"stat": Series (all genes, signed),
                        "up": [genes], "down": [genes]}}.
        """
        labels = labels.reindex(X.index)
        sig: Dict[int, Dict[str, object]] = {}
        n = X.shape[0]

        for c in sorted(labels.unique()):
            a = X.loc[labels == c]
            b = X.loc[labels != c]
            if len(a) < 2 or len(b) < 2:
                continue
            d = a.mean(axis=0) - b.mean(axis=0)
            va, vb = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
            se = np.sqrt(va / len(a) + vb / len(b))
            s0 = np.nanpercentile(se, 10)          # SAM-style offset
            t = d / (se + s0)
            t = t.replace([np.inf, -np.inf], np.nan).fillna(0.0)
            if min_abs_lfc:
                t = t.where(d.abs() >= min_abs_lfc, 0.0)
            up = t.sort_values(ascending=False).head(n_top)
            dn = t.sort_values().head(n_top)
            sig[int(c)] = {
                "stat": t,
                "lfc": d,
                "up": [g for g, v in up.items() if v > 0],
                "down": [g for g, v in dn.items() if v < 0],
                "n": int((labels == c).sum()),
            }
        return sig


    # ===========================================================================
    # 4. Tahoe-100M perturbation reference
    # ===========================================================================


    def load_tahoe_de(
        self,
        organs: Optional[Sequence[str]] = ("Pancreas",),
        cell_lines: Optional[Sequence[str]] = None,
        drugs: Optional[Sequence[str]] = None,
        cache_dir: Optional[str] = None,
        stat_col: Optional[str] = None,
        gene_col: Optional[str] = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Pull the Tahoe-100M pseudobulk DE table as a signature matrix.

        Returns
        -------
        S : DataFrame, genes x conditions
            Signed DE statistic per (cell_line, drug[, dose]) vs plate-matched DMSO.
        cond : DataFrame
            Condition metadata (cell line, organ, drug, MOA, targets).

        Notes
        -----
        The DE table lives under `metadata/pseudobulk_differential_expression/`
        (~4.1e9 rows across all conditions), so we push the cell-line and drug
        filters into the parquet scan rather than materialising it. Column names
        have shifted between revisions of the dataset card, so the schema is
        sniffed at runtime; override with `stat_col` / `gene_col` if the sniffing
        picks the wrong one.
        """
        import pyarrow.dataset as ds  # noqa: F401  (kept for the local-path branch)
        from huggingface_hub import snapshot_download
        import duckdb

        local = snapshot_download(
            repo_id=self._TAHOE_REPO,
            repo_type="dataset",
            allow_patterns=["metadata/*"],
            cache_dir=cache_dir,
        )
        de_glob = os.path.join(local, "metadata", "pseudobulk_differential_expression", "*.parquet")
        con = duckdb.connect()

        cl_meta = con.execute(
            f"SELECT * FROM read_parquet('{os.path.join(local,'metadata','cell_line_metadata.parquet')}')"
        ).df()
        drug_meta = con.execute(
            f"SELECT * FROM read_parquet('{os.path.join(local,'metadata','drug_metadata.parquet')}')"
        ).df()

        # --- resolve which cell lines we want ----------------------------------
        sel = cl_meta.copy()
        if organs:
            sel = sel[sel["Organ"].isin(organs)]
        if cell_lines:
            sel = sel[sel["cell_name"].isin(cell_lines)]
        if sel.empty:
            raise ValueError(f"no cell lines match organs={organs} lines={cell_lines}")
        cvcl = sorted(sel["Cell_ID_Cellosaur"].dropna().unique().tolist())

        # --- sniff schema -------------------------------------------------------
        cols = con.execute(
            f"SELECT * FROM read_parquet('{de_glob}') LIMIT 1"
        ).df().columns.tolist()

        def _pick(cands, override):
            if override:
                return override
            for c in cands:
                for col in cols:
                    if col.lower() == c:
                        return col
            for c in cands:
                for col in cols:
                    if c in col.lower():
                        return col
            raise KeyError(f"none of {cands} in DE schema {cols}")

        gcol = _pick(["gene_symbol", "gene", "genes", "symbol"], gene_col)
        scol = _pick(["log2fc", "logfoldchange", "lfc", "score", "stat",
                    "log2_fold_change"], stat_col)
        ccol = _pick(["cell_line_id", "cell_id_cellosaur", "cell_line"], None)
        dcol = _pick(["drug", "drugname_drugconc", "treatment"], None)

        where = [f"{ccol} IN ({','.join(repr(x) for x in cvcl)})"]
        if drugs:
            where.append(f"{dcol} IN ({','.join(repr(x) for x in drugs)})")

        q = f"""
            SELECT {ccol} AS cell_line_id, {dcol} AS drug,
                {gcol} AS gene, {scol} AS stat
            FROM read_parquet('{de_glob}')
            WHERE {' AND '.join(where)}
        """
        df = con.execute(q).df()
        if df.empty:
            raise ValueError("Tahoe DE query returned nothing -- check filters.")

        df["condition"] = df["cell_line_id"] + "|" + df["drug"].astype(str)
        S = df.pivot_table(index="gene", columns="condition",
                        values="stat", aggfunc="mean")

        cond = (df[["condition", "cell_line_id", "drug"]]
                .drop_duplicates()
                .merge(cl_meta[["Cell_ID_Cellosaur", "cell_name", "Organ"]],
                    left_on="cell_line_id", right_on="Cell_ID_Cellosaur",
                    how="left")
                .merge(drug_meta[["drug", "moa-fine", "moa-broad", "targets",
                                "human-approved"]],
                    on="drug", how="left")
                .set_index("condition"))
        return S, cond.loc[S.columns]


    # --- CMap-style connectivity ------------------------------------------------
    def _weighted_es(metric_desc: np.ndarray, hit_mask: np.ndarray,
                    p: float = 1.0) -> float:
        """GSEA weighted running-sum enrichment score."""
        n_hit = int(hit_mask.sum())
        if n_hit == 0 or n_hit == len(hit_mask):
            return 0.0
        w = np.abs(metric_desc) ** p
        nr = w[hit_mask].sum()
        if nr <= 0:
            return 0.0
        inc = np.where(hit_mask, w / nr, 0.0)
        dec = np.where(hit_mask, 0.0, 1.0 / (len(hit_mask) - n_hit))
        run = np.cumsum(inc - dec)
        return float(run[np.argmax(np.abs(run))])


    def _wtcs(ref: pd.Series, up: Sequence[str], down: Sequence[str]) -> float:
        """Weighted connectivity score of a query (up/down) against one reference
        perturbation signature. Sign convention: positive = the drug *recapitulates*
        the cluster program; negative = the drug *reverses* it."""
        r = ref.dropna().sort_values(ascending=False)
        idx = r.index
        m = r.values
        up_m = np.isin(idx, list(up))
        dn_m = np.isin(idx, list(down))
        if up_m.sum() < 5 or dn_m.sum() < 5:
            return np.nan
        es_u = _weighted_es(m, up_m)
        es_d = _weighted_es(m, dn_m)
        return (es_u - es_d) / 2.0 if np.sign(es_u) != np.sign(es_d) else 0.0


    def score_clusters_vs_tahoe(
        sig: Dict[int, Dict],
        S: pd.DataFrame,
        cond: pd.DataFrame,
        *,
        normalize_within: str = "cell_line_id",
    ) -> pd.DataFrame:
        """WTCS + within-cell-line normalised connectivity (NCS) for every
        (cluster, condition) pair.

        Interpretation: rank by ascending `ncs` to get compounds predicted to
        *reverse* a cluster's malignant program; descending to find compounds whose
        transcriptional footprint *is* that program (mechanistic hypothesis for
        what the cluster's cells are doing).
        """
        rows = []
        for c, s in sig.items():
            up, dn = s["up"], s["down"]
            for condition in S.columns:
                rows.append({
                    "cluster": c,
                    "condition": condition,
                    "wtcs": _wtcs(S[condition], up, dn),
                })
        R = pd.DataFrame(rows)
        R = R.join(cond, on="condition")

        # NCS: divide by the mean of the same-signed WTCS within each
        # (cluster, cell line) stratum, so magnitudes are comparable across lines.
        R["ncs"] = np.nan
        for _, idx in R.groupby(["cluster", normalize_within]).groups.items():
            v = R.loc[idx, "wtcs"]
            pos, neg = v[v > 0], v[v < 0]
            mu_p = pos.mean() if len(pos) else np.nan
            mu_n = abs(neg.mean()) if len(neg) else np.nan
            ncs = pd.Series(np.nan, index=v.index)
            if np.isfinite(mu_p) and mu_p > 0:
                ncs[v > 0] = v[v > 0] / mu_p
            if np.isfinite(mu_n) and mu_n > 0:
                ncs[v < 0] = v[v < 0] / mu_n
            ncs[v == 0] = 0.0
            R.loc[idx, "ncs"] = ncs.values

        return R.sort_values(["cluster", "ncs"]).reset_index(drop=True)


    # ===========================================================================
    # 5. scGPT embedding space
    # ===========================================================================
    def scgpt_embed_pseudobulk(
        X_counts: pd.DataFrame,
        model_dir: str,
        *,
        obs: Optional[pd.DataFrame] = None,
        gene_col: str = "index",
        batch_size: int = 16,
        device: str = "cuda",
    ) -> pd.DataFrame:
        """Embed pseudobulk profiles (samples x genes, *raw-ish counts*) with a
        pretrained scGPT and return an n x d embedding frame.

        Use this to place (a) malignant-compartment profiles and (b) Tahoe
        perturbed/DMSO pseudobulks in one space, then measure cosine distance
        between a cluster centroid and each drug's shift vector.

        Caveats worth stating in a methods section
        ------------------------------------------
        * scGPT was pretrained on single cells; pseudobulk is out of distribution.
        It works in practice but the embedding is not calibrated for it.
        * The scGPT vocabulary is protein-coding-dominant. Most of your Program-1
        HOX-antisense panel will not be tokenised -- check the "match n/m genes"
        log line. Anything you conclude in this space is about the coding
        correlate of the program, not the lncRNAs themselves.
        """
        import anndata as ad
        import scgpt as scg

        obs = obs if obs is not None else pd.DataFrame(index=X_counts.index)
        adata = ad.AnnData(
            X=np.asarray(X_counts.values, dtype=np.float32),
            obs=obs.loc[X_counts.index].copy(),
            var=pd.DataFrame(index=X_counts.columns),
        )
        emb = scg.tasks.embed_data(
            adata,
            model_dir,
            gene_col=gene_col,
            obs_to_save=list(obs.columns) or None,
            batch_size=batch_size,
            device=device,
            return_new_adata=True,
        )
        return pd.DataFrame(np.asarray(emb.X), index=X_counts.index)


    def cosine_to_perturbation_shifts(
        emb: pd.DataFrame,
        cluster_labels: pd.Series,
        pert_emb: pd.DataFrame,
        ctrl_emb: pd.DataFrame,
    ) -> pd.DataFrame:
        """Cosine between each cluster centroid (in scGPT space, centred on the
        cohort mean) and each drug's (treated - DMSO) shift vector."""
        Z = emb.sub(emb.mean(axis=0), axis=1)
        cents = Z.groupby(cluster_labels.reindex(Z.index)).mean()
        shifts = pert_emb.values - ctrl_emb.reindex(pert_emb.index).values

        def _cos(A, B):
            A = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
            B = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)
            return A @ B.T

        return pd.DataFrame(_cos(cents.values, shifts),
                            index=cents.index, columns=pert_emb.index)


