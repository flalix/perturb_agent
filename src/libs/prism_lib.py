#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/07/31
# Udated  on 2026/07/31
# @author: Flavio Lichtenstein
# @local: Home sweet home


"""
paad_deconv.py
==============
Bulk RNA-seq deconvolution for PAAD (CPTAC3 + TCGA), producing

  (1) theta : cell-type / cell-state fractions per sample
  (2) Z     : cell-type-specific *expression* per sample (gene x sample, per cell type)

Core engine is a pure-NumPy/PyTorch re-implementation of the BayesPrism
generative model solved by its analytic fixed point (this is exactly what
InstaPrism does: the posterior mean of the BayesPrism Gibbs sampler with a
flat Dirichlet prior coincides with the EM fixed point of the multinomial
mixture, so you get BayesPrism results ~100-1000x faster and deterministically).

Model
-----
    X[s, g]  ~  Multinomial( N_s , sum_k theta[s,k] * phi[k,g] )
    Z[s,k,g] =  E[ counts of gene g in sample s coming from cell type k ]

E-step:  Z[s,k,g] = X[s,g] * theta[s,k] phi[k,g] / sum_k' theta[s,k'] phi[k',g]
M-step:  theta[s,k] propto sum_g Z[s,k,g]   (+ Dirichlet alpha - 1 pseudocount)

Stage 2 ("updated reference", the part that matters for tumors): the malignant
compartment gets a *sample-specific* profile psi_mal[s,:] estimated from Z,
non-malignant types are merged into a single psi_env profile, and the model is
re-fit. This removes the bias caused by using one fixed malignant reference for
tumors that are transcriptionally heterogeneous (classical vs basal-like PDAC).

Author: written for a PAAD CPTAC3/TCGA project
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

import re
import warnings

from dataclasses import dataclass, field
from typing import Iterable, Sequence


try:
    import torch

    _HAS_TORCH = True
except Exception:  # pragma: no cover
    _HAS_TORCH = False


# =============================================================================
# 1. BULK MATRIX ASSEMBLY  (df_tumor + df_normal + metadata -> counts matrix)
# =============================================================================

# genes that must be removed before deconvolution (BayesPrism cleanup.genes)
_DROP_REGEX = re.compile(
    r"^(MT-|MTRNR|MTATP|MTCO|MTCYB|MTND|MTRF|"      # mitochondrial
    r"RPL\d|RPS\d|RPLP\d|RPSA$|"                     # cytosolic ribosomal
    r"MRPL\d|MRPS\d|"                                # mito-ribosomal
    r"HB[ABDEGQZ]\d?$|"                              # hemoglobin (blood contam.)
    r"MALAT1$|NEAT1$)"                               # nuclear-retained, huge, driven by protocol
)

DEFAULT_BIOTYPES = ("protein_coding",)

class PRISM(object):
    def __init__(self, root0: Path, root0_data: Path):

        self.ANNOT_COLS = ["geneid", "symbol", "biotype"]


    def build_bulk_matrix(self,
        df_tumor: pd.DataFrame,
        df_normal: pd.DataFrame,
        df_metadata: pd.DataFrame | None = None,
        keep_biotypes: Sequence[str] | None = DEFAULT_BIOTYPES,
        gene_key: str = "symbol",
        drop_regex: re.Pattern = _DROP_REGEX,
        min_count: int = 10,
        min_samples_frac: float = 0.05,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Merge df_tumor / df_normal into one raw-count matrix (genes x samples).

        IMPORTANT: BayesPrism/InstaPrism take **raw counts**, not TPM/FPKM/CPM and
        not log-transformed values. Do NOT normalise here.

        Returns
        -------
        df_bulk : DataFrame (genes x samples), integer counts
        df_meta : DataFrame (samples x metadata), aligned & ordered to df_bulk.columns
        """
        id_cols = [c for c in self.ANNOT_COLS if c in df_tumor.columns and c in df_normal.columns]
        if gene_key not in id_cols:
            raise ValueError(f"`{gene_key}` must be an identifier column in both frames.")

        dft = df_tumor.set_index(id_cols)
        dfn = df_normal.set_index(id_cols)

        # align on the full gene identifier (geneid+symbol+biotype)
        common = dft.index.intersection(dfn.index)
        if len(common) < 0.8 * min(len(dft), len(dfn)):
            warnings.warn(
                f"Only {len(common)} genes shared between tumor and normal "
                f"({len(dft)} / {len(dfn)}). Check that both came from the same annotation."
            )
        df_bulk = pd.concat([dft.loc[common], dfn.loc[common]], axis=1)
        df_bulk = df_bulk.reset_index()

        # --- biotype filter -------------------------------------------------
        if keep_biotypes is not None and "biotype" in df_bulk.columns:
            df_bulk = df_bulk[df_bulk["biotype"].isin(keep_biotypes)]

        # --- collapse duplicated symbols by summing counts ------------------
        genes = df_bulk[gene_key].astype(str)
        df_mat = df_bulk.drop(columns=[c for c in id_cols if c in df_bulk.columns])
        df_mat = df_mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        df_mat.index = genes.values
        df_mat = df_mat.groupby(level=0).sum()

        # --- drop MT / ribosomal / hemoglobin -------------------------------
        df_mat = df_mat[~df_mat.index.to_series().str.match(drop_regex)]

        # --- expression filter ----------------------------------------------
        keep = (df_mat >= min_count).mean(axis=1) >= min_samples_frac
        df_mat = df_mat.loc[keep]

        df_mat = df_mat.round().astype(np.int64)
        df_mat.index.name = gene_key

        # --- metadata alignment ---------------------------------------------
        if df_metadata is not None:
            df_meta = df_metadata.copy()
            df_meta.index = df_meta.index.astype(str)
            missing = [c for c in df_mat.columns if c not in df_meta.index]
            if missing:
                warnings.warn(f"{len(missing)} df_bulk samples absent from metadata: {missing[:5]}")
            shared = [c for c in df_mat.columns if c in df_meta.index]
            df_mat = df_mat[shared]
            df_meta = df_meta.loc[shared]
        else:
            df_meta = pd.DataFrame(index=df_mat.columns)

        return df_mat, df_meta


    def qc_report(self, df_bulk: pd.DataFrame, df_meta: pd.DataFrame | None = None) -> pd.DataFrame:
        """
        Per-sample QC. Use this to confirm/deny the batch split visible in the
        log10(counts) boxplots (CPTAC3 vs TCGA is the obvious suspect).
        """
        lib = df_bulk.sum(axis=0)
        qc = pd.DataFrame(
            {
                "library_size": lib,
                "log10_library": np.log10(lib.replace(0, np.nan)),
                "n_detected": (df_bulk > 0).sum(axis=0),
                "frac_zero": (df_bulk == 0).mean(axis=0),
                "median_nonzero": df_bulk.replace(0, np.nan).median(axis=0),
                "q99": df_bulk.quantile(0.99, axis=0),
                # complexity: how much of the library is eaten by the top 50 genes
                "top50_frac": df_bulk.apply(lambda c: c.nlargest(50).sum() / max(c.sum(), 1), axis=0),
            }
        )
        if df_meta is not None and len(df_meta.columns):
            qc = qc.join(df_meta, how="left")
        return qc


    # =============================================================================
    # 2. SINGLE-CELL REFERENCE  (scRNA-seq -> pseudobulk profiles per cell state)
    # =============================================================================

    # Recommended PAAD scRNA references (download as .h5ad):
    #   Peng et al. 2019, Cell Res  - CRA001160 (24 PDAC + 11 normal pancreas, ~57k cells)
    #   Steele et al. 2020, Nat Cancer - GSE155698
    #   Werba et al. 2023, Nat Commun  - GSE205013 (treatment-naive + treated)
    #   Lin et al. 2020 - GSE154778 (metastatic)
    #   Pre-annotated / harmonised: TISCH2 (PAAD_CRA001160), or the CellxGene PDAC collections.
    #
    # Required .obs columns:
    #   cell_type  : coarse compartment  (malignant, acinar, ductal, endocrine,
    #                fibroblast, stellate, endothelial, T, NK, B, plasma,
    #                myeloid, mast, ...)
    #   cell_state : fine label used for deconvolution
    #                (malignant_classical, malignant_basal, iCAF, myCAF, apCAF,
    #                 CD8_exhausted, CD8_effector, Treg, TAM_C1QC, TAM_SPP1, ...)


    def pseudobulk_reference(self,
        adata,
        state_key: str = "cell_state",
        type_key: str = "cell_type",
        layer: str | None = None,
        min_cells: int = 30,
        drop_regex: re.Pattern = _DROP_REGEX,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """
        Sum raw single-cell counts within each cell state -> reference profiles.

        Parameters
        ----------
        adata : AnnData with RAW COUNTS in .X (or in `layer`).

        Returns
        -------
        ref   : DataFrame (states x genes) of summed raw counts
        s2t   : Series mapping cell_state -> cell_type
        """
        import scipy.sparse as sp

        X = adata.layers[layer] if layer else adata.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        if X.data.size and not np.allclose(X.data[:1000], np.round(X.data[:1000])):
            warnings.warn("Reference .X does not look like raw counts (non-integer values).")

        states = pd.Series(adata.obs[state_key].astype(str).values, name=state_key)
        keep_states = states.value_counts()
        keep_states = keep_states[keep_states >= min_cells].index
        mask = states.isin(keep_states).values

        codes, uniq = pd.factorize(states[mask])
        # one-hot (states x cells) @ (cells x genes) -> (states x genes)
        onehot = sp.csr_matrix(
            (np.ones(codes.size), (codes, np.arange(codes.size))),
            shape=(len(uniq), codes.size),
        )
        prof = np.asarray((onehot @ X[mask]).todense())

        ref = pd.DataFrame(prof, index=pd.Index(uniq, name=state_key), columns=adata.var_names)
        ref = ref.loc[:, ~ref.columns.to_series().str.match(drop_regex)]
        if ref.columns.duplicated().any():           # collapse duplicated symbols
            ref = ref.T.groupby(level=0).sum().T

        s2t = (
            adata.obs.loc[mask, [state_key, type_key]]
            .astype(str)
            .drop_duplicates()
            .set_index(state_key)[type_key]
            .reindex(ref.index)
        )
        return ref, s2t


    def select_genes(self, 
        ref: pd.DataFrame,
        n_per_state: int = 200,
        min_logfc: float = 0.5,
        pseudocount: float = 1.0,
    ) -> list[str]:
        """
        Informative-gene selection (analogue of BayesPrism's get.exp.stat +
        select.marker). Keeps genes with a strong one-vs-rest log2 fold change
        in the CPM-normalised pseudobulk reference.

        Restricting to ~2-5k informative genes is what makes the fixed point both
        fast and stable; it barely changes theta but massively reduces the
        influence of protocol-driven genes.
        """
        cpm = ref.div(ref.sum(axis=1), axis=0) * 1e6
        logc = np.log2(cpm + pseudocount)

        markers: set[str] = set()
        for state in logc.index:
            rest = logc.drop(index=state)
            lfc = logc.loc[state] - rest.max(axis=0)          # vs the best competitor
            lfc = lfc[lfc > min_logfc]
            markers.update(lfc.nlargest(n_per_state).index)
        return sorted(markers)


    # =============================================================================
    # 3. THE ENGINE  (BayesPrism / InstaPrism fixed point)
    # =============================================================================


    def _to_backend(self, *arrays, device: str = "cpu", dtype=np.float64):
        if device != "cpu" and _HAS_TORCH:
            td = torch.float32 if dtype is np.float32 else torch.float64
            return [torch.as_tensor(np.asarray(a), dtype=td, device=device) for a in arrays]
        return [np.asarray(a, dtype=dtype) for a in arrays]


    def prism_em(self,
        X: np.ndarray,
        phi: np.ndarray,
        *,
        alpha: float = 1.0,
        max_iter: int = 500,
        tol: float = 1e-7,
        device: str = "cpu",
        verbose: bool = False,
    ) -> np.ndarray:
        """
        Solve for theta (S x K) given df_bulk counts X (S x G) and reference
        profiles phi (K x G, rows sum to 1).

        Memory is O(S*G) -- Z is never materialised during iteration.

        alpha : Dirichlet concentration on theta (alpha=1 -> flat prior = MLE,
                which is BayesPrism's default).
        """
        S, G = X.shape
        K = phi.shape[0]
        use_torch = device != "cpu" and _HAS_TORCH
        xp = torch if use_torch else np
        X_, phi_ = _to_backend(X, phi, device=device)

        theta = xp.ones((S, K), dtype=X_.dtype) / K
        if use_torch:
            theta = theta.to(device)

        prev = theta
        for it in range(max_iter):
            den = theta @ phi_                       # (S, G) expected mixture
            den = xp.clip(den, 1e-300, None)
            w = X_ / den                             # (S, G)
            theta = theta * (w @ phi_.T)             # (S, K)  <- E+M combined
            if alpha != 1.0:
                theta = theta + (alpha - 1.0)
                theta = xp.clip(theta, 1e-300, None)
            theta = theta / theta.sum(axis=1, keepdims=True)

            delta = float(xp.abs(theta - prev).max())
            if verbose and it % 25 == 0:
                print(f"  iter {it:4d}  max|dtheta| = {delta:.3e}")
            if delta < tol:
                break
            prev = theta

        return theta.cpu().numpy() if use_torch else np.asarray(theta)


    def posterior_Z(
        X: np.ndarray,
        theta: np.ndarray,
        phi: np.ndarray,
        *,
        chunk: int = 32,
    ) -> np.ndarray:
        """
        Cell-type-specific expression: Z (S x K x G), in count units.
        Chunked over samples to bound memory.
        """
        S, G = X.shape
        K = phi.shape[0]
        Z = np.empty((S, K, G), dtype=np.float32)
        for i in range(0, S, chunk):
            sl = slice(i, min(i + chunk, S))
            num = theta[sl, :, None] * phi[None, :, :]          # (s, K, G)
            den = np.clip(num.sum(axis=1, keepdims=True), 1e-300, None)
            Z[sl] = (X[sl, None, :] * num / den).astype(np.float32)
        return Z


    def update_reference(self,
        Z: np.ndarray,
        malignant_idx: Sequence[int],
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        BayesPrism stage 2. From the stage-1 posterior Z, build

        psi_mal[s, g] : SAMPLE-SPECIFIC malignant profile (rows sum to 1)
        psi_env[g]    : one merged non-malignant profile   (sums to 1)

        This is the step that makes BayesPrism work on tumors: the malignant
        reference is no longer forced to a fixed scRNA profile, so classical vs
        basal-like PDAC heterogeneity does not leak into the stromal/immune
        fractions.
        """
        mal = np.asarray(malignant_idx, dtype=int)
        env = np.setdiff1d(np.arange(Z.shape[1]), mal)

        psi_mal = Z[:, mal, :].sum(axis=1)                        # (S, G)
        psi_mal = psi_mal / np.clip(psi_mal.sum(axis=1, keepdims=True), 1e-300, None)

        psi_env = Z[:, env, :].sum(axis=(0, 1))                   # (G,)
        psi_env = psi_env / np.clip(psi_env.sum(), 1e-300, None)

        return psi_mal.astype(np.float64), psi_env.astype(np.float64)


    def _refit_two_component(self,
        X: np.ndarray,
        psi_mal: np.ndarray,
        psi_env: np.ndarray,
        max_iter: int = 500,
        tol: float = 1e-8,
    ) -> np.ndarray:
        """Per-sample 2-component EM with a sample-specific malignant profile."""
        S = X.shape[0]
        theta2 = np.full((S, 2), 0.5)
        for s in range(S):
            phi_s = np.vstack([psi_mal[s], psi_env])
            theta2[s] = self.prism_em(X[s : s + 1], phi_s, max_iter=max_iter, tol=tol)[0]
        return theta2


    # =============================================================================
    # 4. ORCHESTRATOR
    # =============================================================================


    @dataclass
    class DeconvResult:
        theta: pd.DataFrame                     # samples x cell states (final)
        theta_stage1: pd.DataFrame              # samples x cell states (initial)
        theta_type: pd.DataFrame                # samples x coarse cell types
        tumor_purity: pd.Series                 # malignant fraction
        genes: list[str] = field(repr=False, default_factory=list)
        Z: np.ndarray | None = field(repr=False, default=None)      # S x K x G
        states: list[str] = field(repr=False, default_factory=list)

        def cell_type_expression(self, state: str, *, cpm: bool = True) -> pd.DataFrame:
            """Return the deconvolved expression of one cell state: genes x samples."""
            if self.Z is None:
                raise ValueError("Re-run with return_Z=True.")
            k = self.states.index(state)
            Zk = pd.DataFrame(self.Z[:, k, :].T, index=self.genes, columns=self.theta.index)
            if cpm:
                Zk = Zk.div(Zk.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
            return Zk


    def run_bayesprism(self, 
        df_bulk: pd.DataFrame,
        ref: pd.DataFrame,
        state_to_type: pd.Series,
        malignant_types: Iterable[str] = ("malignant", "tumor", "epithelial_malignant"),
        gene_subset: Sequence[str] | None = None,
        update_malignant_reference: bool = True,
        return_Z: bool = True,
        device: str = "cpu",
        verbose: bool = True,
    ) -> DeconvResult:
        """
        Full two-stage BayesPrism on `df_bulk` (genes x samples) with reference
        `ref` (states x genes).
        """
        # ---- align genes ----------------------------------------------------
        genes = df_bulk.index.intersection(ref.columns)
        if gene_subset is not None:
            genes = genes.intersection(pd.Index(gene_subset))
        genes = genes.sort_values()
        if len(genes) < 500:
            warnings.warn(f"Only {len(genes)} genes shared -- deconvolution will be unstable.")
        if verbose:
            print(f"[prism] {df_bulk.shape[1]} samples x {len(genes)} genes x {ref.shape[0]} states")

        X = df_bulk.loc[genes].T.to_numpy(dtype=np.float64)           # (S, G)
        R = ref[genes].to_numpy(dtype=np.float64)                  # (K, G)
        R = R + 1e-8
        phi = R / R.sum(axis=1, keepdims=True)

        states = list(ref.index)
        types = state_to_type.reindex(states).astype(str)
        mal_idx = [i for i, s in enumerate(states) if types.iloc[i].lower() in
                {m.lower() for m in malignant_types}]
        if not mal_idx:
            warnings.warn("No malignant state found; stage 2 reference update disabled.")
            update_malignant_reference = False

        # ---- stage 1 --------------------------------------------------------
        theta1 = self.prism_em(X, phi, device=device, verbose=verbose)
        theta_final = theta1
        Z = self.posterior_Z(X, theta1, phi) if (return_Z or update_malignant_reference) else None

        # ---- stage 2 --------------------------------------------------------
        if update_malignant_reference:
            if verbose:
                print("[prism] stage 2: sample-specific malignant reference")
            psi_mal, psi_env = self.update_reference(Z, mal_idx)
            theta2 = self._refit_two_component(X, psi_mal, psi_env)      # (S, 2)

            env_idx = np.setdiff1d(np.arange(len(states)), mal_idx)
            theta_final = theta1.copy()
            # rescale malignant states to the updated purity, env states to 1 - purity
            mal_share = theta1[:, mal_idx].sum(axis=1, keepdims=True)
            env_share = theta1[:, env_idx].sum(axis=1, keepdims=True)
            theta_final[:, mal_idx] = theta1[:, mal_idx] / np.clip(mal_share, 1e-12, None) * theta2[:, [0]]
            theta_final[:, env_idx] = theta1[:, env_idx] / np.clip(env_share, 1e-12, None) * theta2[:, [1]]
            theta_final /= theta_final.sum(axis=1, keepdims=True)

            if return_Z:  # recompute Z under the updated mixture
                Z = posterior_Z(X, theta_final, phi)

        # ---- package --------------------------------------------------------
        idx = df_bulk.columns
        th = pd.DataFrame(theta_final, index=idx, columns=states)
        th1 = pd.DataFrame(theta1, index=idx, columns=states)
        th_type = th.T.groupby(types.values).sum().T
        purity = th.iloc[:, mal_idx].sum(axis=1) if mal_idx else pd.Series(np.nan, index=idx)
        purity.name = "tumor_purity"

        return DeconvResult(
            theta=th,
            theta_stage1=th1,
            theta_type=th_type,
            tumor_purity=purity,
            genes=list(genes),
            Z=Z if return_Z else None,
            states=states,
        )


    # =============================================================================
    # 5. BASELINES / CROSS-CHECKS
    # =============================================================================


    def nnls_deconvolve(df_bulk: pd.DataFrame, ref: pd.DataFrame, genes=None) -> pd.DataFrame:
        """Fast NNLS baseline on CPM signature space (dtangle/CIBERSORT-like)."""
        from scipy.optimize import nnls

        g = df_bulk.index.intersection(ref.columns)
        if genes is not None:
            g = g.intersection(pd.Index(genes))
        B = (df_bulk.loc[g] / df_bulk.loc[g].sum() * 1e6).to_numpy()
        S = (ref[g].div(ref[g].sum(axis=1), axis=0) * 1e6).T.to_numpy()
        out = np.vstack([nnls(S, B[:, j])[0] for j in range(B.shape[1])])
        out = out / np.clip(out.sum(axis=1, keepdims=True), 1e-12, None)
        return pd.DataFrame(out, index=df_bulk.columns, columns=ref.index)


    def run_instaprism_r(df_bulk: pd.DataFrame, ref: pd.DataFrame, state_to_type: pd.Series):
        """
        Canonical implementation via rpy2, for validation of the Python engine.

            R: remotes::install_github("humengying0907/InstaPrism")
            R: devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
        """
        from rpy2.robjects import numpy2ri, pandas2ri, r
        from rpy2.robjects.packages import importr

        numpy2ri.activate(); pandas2ri.activate()
        ip = importr("InstaPrism")
        g = df_bulk.index.intersection(ref.columns)
        r.assign("df_bulk", df_bulk.loc[g].to_numpy())
        r.assign("ref", ref[g].T.to_numpy())
        r.assign("ct", np.array(state_to_type.reindex(ref.index).astype(str)))
        r.assign("cs", np.array(ref.index.astype(str)))
        r.assign("gn", np.array(g.astype(str)))
        r("""
        rownames(df_bulk) <- gn; colnames(df_bulk) <- paste0("S", 1:ncol(df_bulk))
        rownames(ref)  <- gn; colnames(ref)  <- cs
        res <- InstaPrism(bulk_Expr = df_bulk, refPhi_cs = new("refPhi_cs",
                            phi.cs = ref, cell.type.labels = ct, cell.state.labels = cs),
                            dfn.iter = 100, update = TRUE)
        """)
        theta = np.asarray(r("dft(res@Post.updated.ct@theta)"))
        return pd.DataFrame(theta, index=df_bulk.columns)


    # =============================================================================
    # 6. DOWNSTREAM: purity-corrected PDAC subtyping on the malignant compartment
    # =============================================================================

    # Representative marker sets (extend with the full published lists).
    MOFFITT_BASAL = ["KRT17", "KRT5", "KRT6A", "KRT6C", "KRT14", "KRT15", "S100A2",
                    "SPRR1B", "SPRR3", "TNS4", "DHRS9", "LY6D", "SCEL", "SERPINB3",
                    "SERPINB4", "CST6", "AREG", "FAM83A", "GPR87", "VGLL1"]
    MOFFITT_CLASSICAL = ["GATA6", "TFF1", "TFF2", "TFF3", "AGR2", "LGALS4", "CTSE",
                        "BTNL8", "ANXA10", "CEACAM6", "LYZ", "REG4", "TSPAN8",
                        "ST6GALNAC1", "MYO1A", "CLRN3", "VSIG2", "SPINK4"]
    MOFFITT_ACTIVATED_STROMA = ["SPARC", "COL10A1", "COL11A1", "POSTN", "FN1",
                                "INHBA", "WNT2", "SFRP2", "THBS2", "CTHRC1", "FAP"]
    MOFFITT_NORMAL_STROMA = ["ACTA2", "MYH11", "DES", "VIM", "IGF1", "RSPO3", "SPOCK1"]


    def score_signature(expr: pd.DataFrame, gene_set: Sequence[str], *, log: bool = True) -> pd.Series:
        """Mean z-score of a gene set across samples (expr = genes x samples, CPM)."""
        g = [x for x in gene_set if x in expr.index]
        if not g:
            return pd.Series(np.nan, index=expr.columns)
        sub = np.log2(expr.loc[g] + 1) if log else expr.loc[g]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
        return z.mean(axis=0)


    def subtype_malignant(res: DeconvResult, malignant_state: str) -> pd.DataFrame:
        """
        Classical vs basal-like call made on the *deconvolved malignant* expression
        instead of df_bulk -- which is the whole point: df_bulk subtype calls are
        confounded by stromal content in PDAC.
        """
        Zmal = res.cell_type_expression(malignant_state, cpm=True)
        out = pd.DataFrame(
            {
                "classical_score": score_signature(Zmal, MOFFITT_CLASSICAL),
                "basal_score": score_signature(Zmal, MOFFITT_BASAL),
                "tumor_purity": res.tumor_purity,
            }
        )
        out["basal_minus_classical"] = out["basal_score"] - out["classical_score"]
        out["subtype"] = np.where(out["basal_minus_classical"] > 0, "basal-like", "classical")
        return out

