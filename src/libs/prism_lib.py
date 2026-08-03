#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/07/31
# Udated  on 2026/07/31
# @author: Flavio Lichtenstein
# @local: Home sweet home


"""
first study
===========
written for a PAAD CPTAC3/TCGA
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


"""

from __future__ import annotations

from pathlib import Path
import datetime as _dt
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import h5py

import re
import json
import os
from dataclasses import asdict, fields, dataclass, is_dataclass, field
import warnings

import anndata as ad

import scipy.sparse as sp
from scipy.optimize import nnls

from libs.cBioPortal_lib import cBioPortal
from libs.Basic import pdreadcsv, pdwritecsv, create_dir

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

DEFAULT_BIOTYPES = ("protein_coding", "lncRNA", "miRNA",)

class PRISM(object):
    def __init__(self, root0: Path, root0_data: Path):

        self.cbio = cBioPortal(root0=root0, root0_data=root0_data, memory_restriction=False)

        self.ANNOT_COLS = ["geneid", "symbol", "biotype"]

        self.root_colab = Path()
        self.root_gtex  = Path()

        self.root_mprog = Path()
        self.root_singc = Path()

        self.fname_bulk = "bulk_matrix.tsv"
        self.fname_meta = "bulk_metadata.tsv"

        # =============================================================================
        # 6. DOWNSTREAM: purity-corrected PDAC subtyping on the malignant compartment
        # =============================================================================

        # Representative marker sets (extend with the full published lists).
        self.MOFFITT_BASAL = ["KRT17", "KRT5", "KRT6A", "KRT6C", "KRT14", "KRT15", "S100A2",
                        "SPRR1B", "SPRR3", "TNS4", "DHRS9", "LY6D", "SCEL", "SERPINB3",
                        "SERPINB4", "CST6", "AREG", "FAM83A", "GPR87", "VGLL1"]
        self.MOFFITT_CLASSICAL = ["GATA6", "TFF1", "TFF2", "TFF3", "AGR2", "LGALS4", "CTSE",
                            "BTNL8", "ANXA10", "CEACAM6", "LYZ", "REG4", "TSPAN8",
                            "ST6GALNAC1", "MYO1A", "CLRN3", "VSIG2", "SPINK4"]
        self.MOFFITT_ACTIVATED_STROMA = ["SPARC", "COL10A1", "COL11A1", "POSTN", "FN1",
                                    "INHBA", "WNT2", "SFRP2", "THBS2", "CTHRC1", "FAP"]
        self.MOFFITT_NORMAL_STROMA = ["ACTA2", "MYH11", "DES", "VIM", "IGF1", "RSPO3", "SPOCK1"]


    def build_bulk_matrix(self,
        df_tumor: pd.DataFrame,
        df_normal: pd.DataFrame,
        df_metadata: pd.DataFrame | None = None,
        keep_biotypes: Sequence[str] | None = DEFAULT_BIOTYPES,
        gene_key: str = "symbol",
        drop_regex: re.Pattern = _DROP_REGEX,
        min_count: int = 10,
        min_samples_frac: float = 0.05,
        force: bool = False,
        verbose: bool = False,
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

        filename_bulk = self.root_singc / self.fname_bulk
        filename_meta = self.root_singc / self.fname_meta

        if filename_bulk.exists() and filename_meta.exists() and not force:
            df_bulk = pdreadcsv(self.fname_bulk, self.root_singc, index_col=0, verbose=verbose)
            df_meta = pdreadcsv(self.fname_meta, self.root_singc, index_col=0, verbose=verbose)

            return df_bulk, df_meta

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
        df_bulk = df_bulk.drop(columns=[c for c in id_cols if c in df_bulk.columns])
        df_bulk = df_bulk.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        df_bulk.index = genes.values
        df_bulk = df_bulk.groupby(level=0).sum()

        # --- drop MT / ribosomal / hemoglobin -------------------------------
        df_bulk = df_bulk[~df_bulk.index.to_series().str.match(drop_regex)]

        # --- expression filter ----------------------------------------------
        keep = (df_bulk >= min_count).mean(axis=1) >= min_samples_frac
        df_bulk = df_bulk.loc[keep]

        df_bulk = df_bulk.round().astype(np.int64)
        df_bulk.index.name = gene_key

        # --- metadata alignment ---------------------------------------------
        if df_metadata is not None:
            df_meta = df_metadata.copy()
            df_meta.index = df_meta.index.astype(str)
            missing = [c for c in df_bulk.columns if c not in df_meta.index]
            if missing:
                warnings.warn(f"{len(missing)} df_bulk samples absent from metadata: {missing[:5]}")
            shared = [c for c in df_bulk.columns if c in df_meta.index]
            df_bulk = df_bulk[shared]
            df_meta = df_meta.loc[shared]
        else:
            df_meta = pd.DataFrame(index=df_bulk.columns)

        _ = pdwritecsv(df_bulk, self.fname_bulk, self.root_singc, index=True, verbose=verbose)
        _ = pdwritecsv(df_meta, self.fname_meta, self.root_singc, index=True, verbose=verbose)

        return df_bulk, df_meta


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


    """
    =============================================================================
    2. SINGLE-CELL REFERENCE  (scRNA-seq -> pseudobulk profiles per cell state)
    =============================================================================

    Recommended PAAD scRNA references (download as .h5ad):
      Peng et al. 2019, Cell Res  - CRA001160 (24 PDAC + 11 normal pancreas, ~57k cells)
      Steele et al. 2020, Nat Cancer - GSE155698
      Werba et al. 2023, Nat Commun  - GSE205013 (treatment-naive + treated)
      Lin et al. 2020 - GSE154778 (metastatic)
      Pre-annotated / harmonised: TISCH2 (PAAD_CRA001160), or the CellxGene PDAC collections.
    #
    Required .obs columns:
      cell_type  : coarse compartment  (malignant, acinar, ductal, endocrine,
                   fibroblast, stellate, endothelial, T, NK, B, plasma,
                   myeloid, mast, ...)
      cell_state : fine label used for deconvolution
                   (malignant_classical, malignant_basal, iCAF, myCAF, apCAF,
                   CD8_exhausted, CD8_effector, Treg, TAM_C1QC, TAM_SPP1, ...)
    """
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


    '''
    =============================================================================
    3. THE ENGINE  (BayesPrism / InstaPrism fixed point)
    =============================================================================
    '''

    def _to_backend(self, *arrays, device: str = "cpu", dtype=np.float64):
        if device != "cpu" and _HAS_TORCH:
            td = torch.float32 if dtype is np.float32 else torch.float64
            return [torch.as_tensor(np.asarray(a), dtype=td, device=device) for a in arrays]
        return [np.asarray(a, dtype=dtype) for a in arrays]


    def prism_em(self,
        X: np.ndarray,
        phi: np.ndarray,
        alpha: float = 1.0,
        max_iter: int = 1000,
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

        try:
            S, G = X.shape
        except Exception as e:
            print(f"Error: {e}")
            print(type(X))
            raise ValueError("\n\n-------- stop ------------\n\n")


        K = phi.shape[0]
        use_torch = device != "cpu" and _HAS_TORCH
        xp = torch if use_torch else np
        X_, phi_ = self._to_backend(X, phi, device=device)

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
            self,
            X: np.ndarray,
            theta: np.ndarray,
            phi: np.ndarray,
            chunk: int = 32, ) -> np.ndarray:
        """
        Cell-type-specific expression: Z (S x K x G), in count units.
        Chunked over samples to bound memory.
        """
        try:
            S, G = X.shape
        except Exception as e:
            print(f"Error: {e}")
            print(type(X))
            print("\n------------\n")
            print(X)
            print("\n------------\n")
            raise ValueError("\n\n-------- stop ------------\n\n")


        K = phi.shape[0]
        Z = np.empty((S, K, G), dtype=np.float32)
        for i in range(0, S, chunk):
            sl = slice(i, min(i + chunk, S))
            num = theta[sl, :, None] * phi[None, :, :]          # (s, K, G)
            den = np.clip(num.sum(axis=1, keepdims=True), 1e-300, None)
            Z[sl] = (X[sl, None, :] * num / den).astype(np.float32)
        return Z


    def full_Z(self, res: DeconvResult, 
               df_bulk: pd.DataFrame, 
               ref: pd.DataFrame, 
               chunk: int = 32) -> tuple[np.ndarray, list]:
        """ 
        Recompute Z over ALL shared genes with theta held fixed.
        """

        genes = df_bulk.index.intersection(ref.columns).sort_values()

        X = df_bulk.loc[genes].T.to_numpy(dtype=np.float64)
        R = ref[genes].to_numpy(dtype=np.float64) + 1e-8
        phi = R / R.sum(axis=1, keepdims=True)
        theta = res.theta.to_numpy()

        return self.posterior_Z(X, theta, phi, chunk=chunk), list(genes)


    def state_expression(self, Z: np.ndarray, genes: list[str], res,
                         state: str, cpm: bool = True) -> pd.DataFrame:
        """
        genes x samples for one cell state, from a full-gene Z.
        """
        if state not in res.states:
            raise KeyError(f"{state!r} not among states: {res.states}")
        k  = res.states.index(state)
        Zk = pd.DataFrame(Z[:, k, :].T, index=genes, columns=res.theta.index)
        return Zk.div(Zk.sum(axis=0).replace(0, np.nan), axis=1) * 1e6 if cpm else Zk

    def gene_compartment_share(self, Z: np.ndarray, genes: list[str], res,
                               gene: str) -> pd.Series:
        """
        Fraction of a gene's counts attributed to each cell state.
        """
        if gene not in genes:
            raise KeyError(f"{gene!r} not in the matrix (check _DROP_REGEX and biotype filter).")
        share = Z[:, :, genes.index(gene)].sum(axis=0)
        return (pd.Series(share / share.sum(), index=res.states)
                  .sort_values(ascending=False))

    
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

    '''
    =============================================================================
    4. ORCHESTRATOR
    =============================================================================
    '''

    def open_bayesprism(self, 
        verbose: bool = False,
    ) -> DeconvResult | None:
        self.fname_dec = "deconv.h5ad"
        filename_ad = self.root_singc / self.fname_dec

        if filename_ad.exists():
            res = DeconvResult.load(filename_ad, verbose=verbose)
            return res

        print(f"Could not find deconvolution results: {filename_ad}")
        return None
        

    def run_bayesprism(self, 
        df_bulk: pd.DataFrame,
        meta_desc: dict, 
        ref: pd.DataFrame,
        state_to_type: pd.Series,
        malignant_types: Iterable[str] = ("malignant", "tumor", "epithelial_malignant"),
        gene_subset: Sequence[str] | None = None,
        update_malignant_reference: bool = True,
        return_Z: bool = True,
        device: str = "cpu",
        force: bool = False,
        verbose: bool = False,
    ) -> DeconvResult:
        """
        Full two-stage BayesPrism on `df_bulk` (genes x samples) with reference
        `ref` (states x genes).
        """

        self.fname_dec = "deconv.h5ad"
        filename_ad = self.root_singc / self.fname_dec

        if filename_ad.exists() and not force:
            res = DeconvResult.load(filename_ad, verbose=verbose)
            return res


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
            theta_final = theta1
            # rescale malignant states to the updated purity, env states to 1 - purity
            mal_share = theta1[:, mal_idx].sum(axis=1, keepdims=True)
            env_share = theta1[:, env_idx].sum(axis=1, keepdims=True)
            theta_final[:, mal_idx] = theta1[:, mal_idx] / np.clip(mal_share, 1e-12, None) * theta2[:, [0]]
            theta_final[:, env_idx] = theta1[:, env_idx] / np.clip(env_share, 1e-12, None) * theta2[:, [1]]
            theta_final /= theta_final.sum(axis=1, keepdims=True)

            if return_Z:  # recompute Z under the updated mixture
                Z = self.posterior_Z(X, theta_final, phi)

        # ---- package --------------------------------------------------------
        idx = df_bulk.columns
        th = pd.DataFrame(theta_final, index=idx, columns=states)
        th1 = pd.DataFrame(theta1, index=idx, columns=states)
        th_type = th.T.groupby(types.values).sum().T
        purity = th.iloc[:, mal_idx].sum(axis=1) if mal_idx else pd.Series(np.nan, index=idx)
        purity.name = "tumor_purity"

        res = DeconvResult(
            theta=th,
            theta_stage1=th1,
            theta_type=th_type,
            tumor_purity=purity,
            genes=list(genes),
            Z=Z if return_Z else None,
            states=states,
        )

        res.save(
            filename_ad,
            meta=meta_desc,
            verbose=verbose
        )

        return res

    '''
    =============================================================================
    5. BASELINES / CROSS-CHECKS
    =============================================================================
    '''
    def nnls_deconvolve(self, df_bulk: pd.DataFrame, 
                        ref: pd.DataFrame, genes=None) -> pd.DataFrame:

        """
        Fast NNLS baseline on CPM signature space (dtangle/CIBERSORT-like).
        """

        g = df_bulk.index.intersection(ref.columns)
        if genes is not None:
            g = g.intersection(pd.Index(genes))
        B = (df_bulk.loc[g] / df_bulk.loc[g].sum() * 1e6).to_numpy()
        S = (ref[g].div(ref[g].sum(axis=1), axis=0) * 1e6).T.to_numpy()
        out = np.vstack([nnls(S, B[:, j])[0] for j in range(B.shape[1])])
        out = out / np.clip(out.sum(axis=1, keepdims=True), 1e-12, None)
        return pd.DataFrame(out, index=df_bulk.columns, columns=ref.index)


    def run_instaprism_r(self, df_bulk: pd.DataFrame, ref: pd.DataFrame, state_to_type: pd.Series):
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
        r.assign("df_ct", np.array(state_to_type.reindex(ref.index).astype(str)))
        r.assign("cs", np.array(ref.index.astype(str)))
        r.assign("gn", np.array(g.astype(str)))
        r("""
        rownames(df_bulk) <- gn; colnames(df_bulk) <- paste0("S", 1:ncol(df_bulk))
        rownames(ref)  <- gn; colnames(ref)  <- cs
        res <- InstaPrism(bulk_Expr = df_bulk, refPhi_cs = new("refPhi_cs",
                            phi.cs = ref, cell.type.labels = df_ct, cell.state.labels = cs),
                            dfn.iter = 100, update = TRUE)
        """)
        theta = np.asarray(r("dft(res@Post.updated.df_ct@theta)"))
        return pd.DataFrame(theta, index=df_bulk.columns)


    def score_signature(self, expr: pd.DataFrame, gene_set: Sequence[str], log: bool = True) -> pd.Series:
        """Mean z-score of a gene set across samples (expr = genes x samples, CPM)."""
        g = [x for x in gene_set if x in expr.index]
        if not g:
            return pd.Series(np.nan, index=expr.columns)
        sub = np.log2(expr.loc[g] + 1) if log else expr.loc[g]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
        return z.mean(axis=0)


    def subtype_malignant(self, res: DeconvResult, malignant_state: str) -> pd.DataFrame:
        """
        Classical vs basal-like call made on the *deconvolved malignant* expression
        instead of df_bulk -- which is the whole point: df_bulk subtype calls are
        confounded by stromal content in PDAC.
        """
        Zmal = res.cell_type_expression(malignant_state, cpm=True)
        out = pd.DataFrame(
            {
                "classical_score": self.score_signature(Zmal, self.MOFFITT_CLASSICAL),
                "basal_score": self.score_signature(Zmal, self.MOFFITT_BASAL),
                "tumor_purity": res.tumor_purity,
            }
        )
        out["basal_minus_classical"] = out["basal_score"] - out["classical_score"]
        out["subtype"] = np.where(out["basal_minus_classical"] > 0, "basal-like", "classical")
        return out


    # =============================================================================
    # ROUTE A -- TISCH2  (recommended: Peng 2019 data, already annotated)
    # =============================================================================
    # Portal:  http://tisch.comp-genomics.org/
    # Dataset: PAAD_CRA001160   (Peng et al. 2019, ~57k cells, 13 cell types)
    #          PAAD_GSE165399   (second PDAC cohort, good robustness check)
    #
    # Download from the dataset page:
    #     PAAD_CRA001160_expression.h5          10x-style HDF5 matrix
    #     PAAD_CRA001160_CellMetainfo_table.tsv per-cell annotation
    #
    # The metadata table carries several annotation depths; typical columns are
    # 'Celltype (malignancy)', 'Celltype (major-lineage)', 'Celltype (minor-lineage)',
    # 'Patient', 'Sample', 'Tissue'. Print them before mapping -- do not assume.
    
    
    def load_tisch2(self, h5_filename: str, meta_filename: str, verbose: bool = True):
        """
        Build an AnnData from a TISCH2 .h5 + CellMetainfo table.
        """
        import scanpy as sc
    
        adata = sc.read_10x_h5(h5_filename)
        adata.var_names_make_unique()
    
        df_meta = pd.read_csv(meta_filename, sep="\t", index_col=0)
        if verbose:
            print("metadata columns:", list(df_meta.columns))
    
        common = adata.obs_names.intersection(df_meta.index)
        adata = adata[common].copy()
        adata.obs = adata.obs.join(df_meta.loc[common])
    
        major = "Celltype (major-lineage)"
        minor = "Celltype (minor-lineage)"
        adata.obs["cell_type"] = adata.obs[major].astype(str)
        adata.obs["cell_state"] = adata.obs.get(minor, adata.obs[major]).astype(str)
    
        # BayesPrism needs one compartment flagged as malignant
        adata.obs["cell_type"] = adata.obs["cell_type"].replace({"Malignant": "malignant"})
        return adata
    
    
    '''    
    =============================================================================
    ROUTE B -- CZ CELLxGENE Census (programmatic; run the query, then choose)
    =============================================================================
    # pip install cellxgene-census
    '''    
    
    def census_find_pdac(self, census_version: str = "stable") -> pd.DataFrame:
        """
        List candidate pancreas / PDAC datasets currently in the Census.
        """
        import cellxgene_census
    
        with cellxgene_census.open_soma(census_version=census_version) as census:
            ds = census["census_info"]["datasets"].read().concat().to_pandas()
        hit = ds[
            ds["dataset_title"].str.contains("pancrea|PDAC|ductal", case=False, na=False)
            | ds["collection_name"].str.contains("pancrea|PDAC|ductal", case=False, na=False)
        ]
        return hit[["dataset_id", "collection_name", "dataset_title", "dataset_total_cell_count"]]
    
    
    def census_fetch(self, dataset_id: str, out_h5ad: str | None = None, census_version: str = "stable"):
        """
        Pull one Census dataset as AnnData, or download the untouched source H5AD.
    
        IMPORTANT: the Census `raw` layer holds the counts. `get_anndata` may return
        normalised values in .X depending on the build -- always verify with
        `check_reference` below before using it as a BayesPrism reference.
        """

        import cellxgene_census

        if out_h5ad:
            cellxgene_census.download_source_h5ad(dataset_id, to_path=out_h5ad)
            return ad.read_h5ad(out_h5ad)
    
        with cellxgene_census.open_soma(census_version=census_version) as census:
            return cellxgene_census.get_anndata(
                census,
                organism="Homo sapiens",
                obs_value_filter=f"dataset_id == '{dataset_id}'",
            )
    

    '''    
    =============================================================================
    ROUTE C -- GEO (processed matrices, need your own annotation)
    =============================================================================
       GSE155698  Steele et al. 2020, Nat Cancer      16 PDAC + 3 adjacent
       GSE154778  Lin et al. 2020                     primary + metastatic
       GSE205013  Werba et al. 2023, Nat Commun       treatment-naive + treated
       GSE165399  second TISCH2 PDAC cohort
    
    These ship as 10x barcodes/features/matrix triplets; read with
    scanpy.read_10x_mtx per sample and concatenate, then annotate yourself.
        
    =============================================================================
    VALIDATION -- run this before handing anything to run_bayesprism
    =============================================================================
    '''    
    def check_reference(self, adata,
                        state_key: str = "cell_state", 
                        type_key: str = "cell_type") -> dict:
        """
        Verify the reference satisfies BayesPrism's assumptions. The single most
        common failure is a normalised/log-transformed .X: the model is multinomial
        over counts, and log values silently produce plausible-looking but wrong
        fractions.
        """
        X = adata.X
        sub = X[:200].data if sp.issparse(X) else np.asarray(X[:200]).ravel()
        sub = sub[np.isfinite(sub)]
    
        is_int = bool(np.allclose(sub, np.round(sub))) if sub.size else False
        report = {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "X_looks_like_raw_counts": is_int,
            "X_max": float(sub.max()) if sub.size else np.nan,
            "has_negative_values": bool((sub < 0).any()) if sub.size else False,
            "layers_available": list(adata.layers.keys()),
            "raw_available": adata.raw is not None,
            f"{type_key}_present": type_key in adata.obs,
            f"{state_key}_present": state_key in adata.obs,
        }
        if type_key in adata.obs:
            report["cell_types"] = adata.obs[type_key].value_counts().to_dict()
            report["has_malignant"] = any(
                "malign" in str(t).lower() for t in adata.obs[type_key].unique()
            )
        if state_key in adata.obs:
            vc = adata.obs[state_key].value_counts()
            report["n_states"] = int(vc.size)
            report["states_under_30_cells"] = vc[vc < 30].to_dict()
    
        problems = []
        if not is_int:
            problems.append(
                "X is not integer counts. Try adata.layers['counts'] or adata.raw.X; "
                "if only normalised data exists, this reference is unusable for BayesPrism."
            )
        if report.get("has_negative_values"):
            problems.append("X contains negatives -- this is scaled/z-scored data.")
        if not report.get("has_malignant", True):
            problems.append(
                "No malignant compartment found. Stage-2 reference update will be "
                "disabled and tumor purity cannot be estimated."
            )
        report["problems"] = problems
        return report


    def set_program_and_primary_site(self, prog_id:str, psi_id:str, verbose:bool=False) -> pd.DataFrame:
        df_psi = self.cbio.set_program_and_primary_site(prog_id=prog_id, psi_id=psi_id, verbose=verbose)

        self.root_colab = self.cbio.root_colab
        self.root_gtex  = self.cbio.root_gtex

        self.root_mprog = self.cbio.root_mprog
        self.root_singc = self.cbio.root_singc

        self.df_psi = df_psi
        return df_psi

    def load_and_view_matrix_txt(self, 
                                 fname: str | Path, 
                                 sep: str | None = None, 
                                 nrows: int = 4) -> bool:
        """
        Print the corner of the file. ALWAYS run this first -- confirm the
        orientation (genes as rows vs cells as rows) and the separator.
        """

        filename = self.root_singc / fname

        if not filename.exists():
            print(f"File not found: {filename}")
            return False

        if sep is None:
            with open(filename) as f:
                f.readline()
                line = f.readline()

            counts = {"\t": line.count("\t"),
                      ",": line.count(","),
                      " ": line.count(" "), 
                      ";": line.count(";")}

            sep = max(counts, key=counts.get)

            print("field counts:", {repr(k): v for k, v in counts.items()},
                "-> sep =", repr(sep))

            if counts[sep] == 0:
                print("No delimiter found; file may be compressed or fixed-width.")
                return False
        
        df = pd.read_csv(filename, sep="\t", nrows=nrows, index_col=0)
        print(f"first {nrows} rows x 5 cols:\n{df.iloc[:, :5]}\n")
        print("row labels look like:", list(df.index[:3]))
        print("col labels look like:", list(df.columns[:3]))
        print("dtypes:", df.dtypes.unique())

        return True

    """
    Convert the Peng 2019 GSA deposit into an AnnData ready for
    `paad_deconv.pseudobulk_reference()`.
    
    Input (from ftp://download.big.ac.cn/gsa/CRA001160/):
        count-matrix.txt   2.77 GB dense TSV, genes x cells
        all_celltype.txt   2.1 MB, per-cell annotation
    
    The matrix is dense text: ~20k genes x ~57k cells is ~1.1e9 values, which is
    10-13 GB as a dense float array but well under 1 GB as CSR, since scRNA counts
    are >90% zeros. So it is parsed in row chunks and sparsified incrementally --
    never materialised dense.
    """
    def load_matrix(self, 
        fname: str | Path,
        chunksize: int = 1000,
        dtype=np.float32,
        genes_are_rows: bool = True,
        sep: str = ' ',
        compression:str = "gzip",
        force: bool = False,
        verbose: bool = False, ) -> ad.AnnData:
        """
        Stream the dense TSV into a sparse AnnData (cells x genes).
        head -2 count-matrix.txt | cut -c1-200
        awk 'NR==2{print "tabs:",gsub(/\t/,""); print "spaces:",gsub(/ /,""); print "commas:",gsub(/,/,"")}' count-matrix.txt
        """

        fname_ad = str(fname).replace('.txt', '.h5ad')
        if not fname_ad.endswith('.h5ad'):
            fname_ad += '.h5ad'
        filename_ad = self.root_singc / fname_ad

        if filename_ad.exists() and not force:
            adata = ad.read_h5ad(filename_ad)
            if verbose: 
                print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes | obs: {list(adata.obs.columns)}")
            return adata

        
        filename = self.root_singc / fname

        if not filename.exists():
            print(f"File not found: {filename}")
            return ad.AnnData()

        blocks, labels, cols = [], [], None
        for i, chunk in enumerate(
            pd.read_csv(filename, sep=sep, index_col=0, chunksize=chunksize)
        ):
            if cols is None:
                cols = chunk.columns
            labels.extend(chunk.index.astype(str))
            blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=dtype)))
            if verbose and i % 5 == 0:
                print(f"  chunk {i}: {len(labels):,} rows parsed", flush=True)

        X = sp.vstack(blocks, format="csr")
        del blocks

        if genes_are_rows:
            X = X.T.tocsr()
            obs_names, var_names = list(cols.astype(str)), labels
        else:
            obs_names, var_names = labels, list(cols.astype(str))

        adata = ad.AnnData(
            X,
            obs=pd.DataFrame(index=pd.Index(obs_names, name="cell")),
            var=pd.DataFrame(index=pd.Index(var_names, name="gene")),
        )
        adata.var_names_make_unique()
        if verbose:
            try:
                dens = X.nnz / (X.shape[0] * X.shape[1])
                print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes | "
                    f"density {dens:.3f} | {X.data.nbytes/1e9:.2f} GB in memory")
            except Exception as e:
                print(f"Error occurred while calculating density: {e}")
                print(X.shape)

        print("\n\n------------------- end --------------------\n")

        adata.write_h5ad(filename_ad, compression=compression)
        if verbose:
            print(f"AData saved as {filename_ad},  ({filename_ad.stat().st_size/1e6:.0f} MB), compressed with {compression}")

        return adata


    def attach_celltypes(
        self,
        adata: ad.AnnData,
        fname_celltype: str,
        type_col: str | None = None,
        state_col: str | None = None,
        verbose: bool = True,
    ):
        """
        Join all_celltype.txt onto the AnnData and create `cell_type` / `cell_state`.

        The file's column names are not assumed -- they are printed so you can pass
        the right ones explicitly.
        """

        filename = self.root_singc / fname_celltype

        df_ct = pd.read_csv(filename, sep="\t", index_col=0)
        if verbose:
            print("all_celltype.txt columns:", list(df_ct.columns))
            print(df_ct.head(3))

        if type_col is None:
            cands = [c for c in df_ct.columns if "type" in c.lower() or "cluster" in c.lower()]
            if not cands:
                raise ValueError(f"Could not guess the label column from {list(df_ct.columns)}")
            type_col = cands[0]
            if verbose:
                print(f"using type_col='{type_col}'")

        shared = adata.obs_names.intersection(df_ct.index)
        if verbose:
            print(f"barcode overlap: {len(shared):,} / {adata.n_obs:,}")
        if len(shared) < 0.5 * adata.n_obs:
            raise ValueError(
                "Barcode overlap under 50%. The two files likely use different "
                "barcode conventions (sample prefixes, -1 suffixes). Inspect both "
                "index formats and reconcile before joining."
            )

        adata = adata[shared].copy()
        adata.obs = adata.obs.join(df_ct.loc[shared])
        adata.obs["cell_type"] = adata.obs[type_col].astype(str)
        adata.obs["cell_state"] = adata.obs[state_col or type_col].astype(str)

        # BayesPrism needs one compartment flagged malignant for the stage-2 update
        adata.obs["cell_type"] = adata.obs["cell_type"].replace(
            {"Malignant": "malignant", "Ductal cell type 2": "malignant"}
        )
        if verbose:
            print(adata.obs["cell_type"].value_counts())
        return adata


@dataclass
class DeconvResult:
    theta: pd.DataFrame                     # samples x cell states (final)
    theta_stage1: pd.DataFrame              # samples x cell states (initial)
    theta_type: pd.DataFrame                # samples x coarse cell types
    tumor_purity: pd.Series                 # malignant fraction
    genes: list[str] = field(repr=False, default_factory=list)
    Z: np.ndarray | None = field(repr=False, default=None)      # S x K x G
    states: list[str] = field(repr=False, default_factory=list)
 
    def cell_type_expression(self, state: str, cpm: bool = True) -> pd.DataFrame:
        """Return the deconvolved expression of one cell state: genes x samples."""
        if self.Z is None:
            raise ValueError("Re-run with return_Z=True.")
        k = self.states.index(state)
        Zk = pd.DataFrame(self.Z[:, k, :].T, index=self.genes, columns=self.theta.index)
        if cpm:
            Zk = Zk.div(Zk.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
        return Zk
 
    # ------------------------------------------------------------------ #
    # persistence — thin wrappers over the module-level functions below   #
    # ------------------------------------------------------------------ #
    def save(self, filename, **kw) -> str:
        return save_deconv_result(self, filename, **kw)
 
    @classmethod
    def load(cls, filename, **kw) -> "DeconvResult":
        return load_deconv_result(filename, cls=cls, **kw)
 
    def summary(self) -> pd.DataFrame:
        th = _as_frame(self.theta)
        return pd.DataFrame({
            "mean": th.mean(), "median": th.median(),
            "min": th.min(), "max": th.max(),
            "frac_below_2pct": (th < 0.02).mean(),
        }).round(4)
 
 
# =============================================================================
# 7. PERSISTENCE FOR DeconvResult
# =============================================================================
# fmt="h5ad" : Z flattened to long format (obs = sample x state, X = genes),
#              theta/purity in .uns. Reads back via the "long" layout.
# fmt="h5"   : native 3-D tensor /Z (S, K, G), chunked along the sample axis.
#
# These are MODULE-LEVEL functions on purpose. They take no `self`; putting them
# in a class body is what produced the _as_frame redeclaration.
 
 
def _as_frame(obj, index=None, columns=None, name=None) -> pd.DataFrame | None:
    """theta may be a DataFrame, ndarray or Series depending on the backend."""
    if obj is None:
        return None
    if isinstance(obj, pd.DataFrame):
        return obj
    if isinstance(obj, pd.Series):
        return obj.to_frame(name or "value")
    arr = np.asarray(obj)
    if arr.ndim == 1:
        return pd.DataFrame({name or "value": arr}, index=index)
    return pd.DataFrame(arr, index=index, columns=columns)
 
 
def _samples_of(res) -> list[str]:
    th = res.theta
    if isinstance(th, (pd.DataFrame, pd.Series)):
        return th.index.astype(str).tolist()
    return [f"S{i:04d}" for i in range(np.asarray(th).shape[0])]
 
 
def _vlen_str(g, key: str, values):
    values = list(values)
    if not values:                                   # h5py rejects empty object arrays
        g.create_dataset(key, shape=(0,), dtype=h5py.string_dtype("utf-8"))
        return
    g.create_dataset(key, data=np.asarray(values, dtype=object),
                     dtype=h5py.string_dtype("utf-8"))
 
 
def _default_meta(extra: dict | None = None) -> dict:
    meta = {"created": _dt.datetime.now().isoformat(timespec="seconds"),
            "format_version": "1.0"}
    if extra:
        meta.update(extra)
    return meta
 
 
def save_deconv_result(res,
                       filename: str | Path,
                       fmt: str = "h5ad",
                       dtype=np.float32,
                       compression: str = "gzip",
                       compression_opts: int = 4,
                       meta: dict | None = None,
                       verbose: bool = True) -> str:
    """
    dtype : float32 is right. EM fixed-point values carry nowhere near float64
            precision and it halves the file.
    meta  : provenance — reference dataset, strand mode, cohort, git hash.
    """
    filename = str(filename)                          # accept Path
 
    samples = _samples_of(res)
    states = list(res.states) if res.states else None
    genes = list(res.genes)
 
    th = _as_frame(res.theta, index=samples, columns=states)
    if states is None and th is not None:
        states = th.columns.astype(str).tolist()
 
    th1 = _as_frame(res.theta_stage1, index=samples, columns=states)
    tht = _as_frame(res.theta_type, index=samples)
    pur = _as_frame(res.tumor_purity, index=samples, name="tumor_purity")
 
    meta = _default_meta(meta)
    meta.update(n_samples=len(samples), n_states=len(states or []),
                n_genes=len(genes), has_Z=res.Z is not None)
 
    if fmt == "h5ad":
        filename = _save_h5ad(res, filename, samples, states, genes,
                              th, th1, tht, pur, dtype, compression, meta)
    elif fmt == "h5":
        filename = _save_h5(res, filename, samples, states, genes,
                            th, th1, tht, pur, dtype, compression,
                            compression_opts, meta)
    else:
        raise ValueError("fmt must be 'h5ad' or 'h5'")
 
    if verbose:
        mb = os.path.getsize(filename) / 1e6
        print(f"saved {filename} ({mb:,.1f} MB) | {len(samples)} samples x "
              f"{len(states or [])} states x {len(genes)} genes"
              f"{' | Z included' if res.Z is not None else ' | theta only'}")
    return filename
 
 
def _orient_Z(Z, n_samples: int, n_states: int) -> np.ndarray:
    """Z must be (S, K, G). Fail loudly rather than guess."""
    Z = np.asarray(Z)
    if Z.ndim != 3:
        raise ValueError(f"Z must be 3-D (S,K,G), got shape {Z.shape}")
    if Z.shape[:2] == (n_samples, n_states):
        return Z
    if Z.shape[:2] == (n_states, n_samples):
        warnings.warn("Z looked like (K,S,G); transposing to (S,K,G).")
        return np.transpose(Z, (1, 0, 2))
    raise ValueError(f"Z shape {Z.shape} matches neither (S,K,G)="
                     f"({n_samples},{n_states},G) nor (K,S,G)")
 
 
def _save_h5ad(res, filename, samples, states, genes, th, th1, tht, pur,
               dtype, compression, meta):
    if not filename.endswith(".h5ad"):
        filename += ".h5ad"
 
    if res.Z is not None:
        Z = _orient_Z(res.Z, len(samples), len(states))
        ns, nk, ng = Z.shape
        X = Z.reshape(ns * nk, ng).astype(dtype)
        obs = pd.DataFrame({
            "sample": np.repeat(np.asarray(samples, dtype=object), nk),
            "state": np.tile(np.asarray(states, dtype=object), ns),
        })
        obs.index = pd.Index(obs["sample"].astype(str) + "|" + obs["state"].astype(str))
        obs["sample"] = obs["sample"].astype("category")
        obs["state"] = obs["state"].astype("category")
        adata = ad.AnnData(X=X, obs=obs,
                           var=pd.DataFrame(index=pd.Index(genes, name="geneid")))
    else:
        adata = ad.AnnData(X=np.asarray(th.values, dtype=dtype),
                           obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
                           var=pd.DataFrame(index=pd.Index(states, name="state")))
 
    adata.uns["theta"] = th
    if th1 is not None:
        adata.uns["theta_stage1"] = th1
    if tht is not None:
        adata.uns["theta_type"] = tht
    if pur is not None:
        adata.uns["tumor_purity"] = pur
    adata.uns["states"] = np.asarray(states, dtype=object)
    adata.uns["sample_order"] = np.asarray(samples, dtype=object)
    adata.uns["meta"] = json.dumps(meta)
 
    adata.write_h5ad(filename, compression=compression)
    return filename
 
 
def _save_h5(res, filename, samples, states, genes, th, th1, tht, pur,
             dtype, compression, compression_opts, meta):
    if not filename.endswith((".h5", ".hdf5")):
        filename += ".h5"
 
    with h5py.File(filename, "w") as f:
        f.attrs["meta"] = json.dumps(meta)
        _vlen_str(f, "samples", samples)
        _vlen_str(f, "states", states or [])
        _vlen_str(f, "genes", genes)
 
        if res.Z is not None:
            Z = _orient_Z(res.Z, len(samples), len(states))
            f.create_dataset("Z", data=Z.astype(dtype),
                             chunks=(1, Z.shape[1], Z.shape[2]),
                             compression=compression,
                             compression_opts=compression_opts, shuffle=True)
 
        for key, df in (("theta", th), ("theta_stage1", th1),
                        ("theta_type", tht), ("tumor_purity", pur)):
            if df is None:
                continue
            g = f.create_group(key)
            g.create_dataset("values", data=np.asarray(df.values, dtype=np.float64),
                             compression=compression)
            _vlen_str(g, "index", df.index.astype(str))
            _vlen_str(g, "columns", df.columns.astype(str))
    return filename
 
 
def load_deconv_result(filename: str | Path, cls=None, mmap: bool = False, verbose: bool = False):
    """
    Returns a DeconvResult if `cls` is given, else a plain dict.
 
    mmap=True (fmt='h5' only) leaves Z as an open h5py dataset — the file handle
    stays open under out['_h5file'], so slice what you need and do not let the
    object outlive the process.
    """
    filename = str(filename)
    out = _load_h5ad(filename) if filename.endswith(".h5ad") else _load_h5(filename, mmap=mmap)
 
    # tumor_purity is declared as a Series on the dataclass; it round-trips as a
    # 1-column DataFrame, so squeeze it back.
    pur = out.get("tumor_purity")
    if isinstance(pur, pd.DataFrame) and pur.shape[1] == 1:
        out["tumor_purity"] = pur.iloc[:, 0].rename("tumor_purity")
 
    if cls is None or not is_dataclass(cls):
        return out

    valid = {f.name for f in fields(cls)}

    if verbose:
        mb = os.path.getsize(filename) / 1e6
        print(f"Loaded {filename} ({mb:,.1f} MB)")
            
    return cls(**{k: v for k, v in out.items() if k in valid})
 
 
def _load_h5ad(filename):
    adata = ad.read_h5ad(filename)
    uns = adata.uns
    meta = json.loads(str(uns.get("meta", "{}")))
 
    Z = None
    if {"sample", "state"}.issubset(adata.obs.columns):
        states = [str(s) for s in np.asarray(uns["states"]).ravel()]
        if "sample_order" in uns:
            s_order = [str(s) for s in np.asarray(uns["sample_order"]).ravel()]
        else:
            s_order = adata.obs["sample"].astype(str).drop_duplicates().tolist()
        nk, ng = len(states), adata.n_vars
        if adata.n_obs != len(s_order) * nk:
            raise ValueError(f"obs rows ({adata.n_obs}) != samples x states "
                             f"({len(s_order)} x {nk}); file is inconsistent.")
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        Z = X.reshape(len(s_order), nk, ng)
    else:
        states = [str(s) for s in adata.var_names]
 
    return dict(theta=uns.get("theta"), theta_stage1=uns.get("theta_stage1"),
                theta_type=uns.get("theta_type"), tumor_purity=uns.get("tumor_purity"),
                genes=adata.var_names.astype(str).tolist(), Z=Z,
                states=states, meta=meta)
 
 
def _load_h5(filename, mmap=False):
    f = h5py.File(filename, "r")
    meta = json.loads(f.attrs.get("meta", "{}"))
 
    def _dec(v):
        return v.decode() if isinstance(v, bytes) else str(v)
 
    def _read_df(key):
        if key not in f:
            return None
        g = f[key]
        return pd.DataFrame(g["values"][:],
                            index=[_dec(s) for s in g["index"][:]],
                            columns=[_dec(s) for s in g["columns"][:]])
 
    def _read_str(key):
        return [_dec(s) for s in f[key][:]]
 
    Z = f["Z"] if (mmap and "Z" in f) else (f["Z"][:] if "Z" in f else None)
 
    out = dict(theta=_read_df("theta"), theta_stage1=_read_df("theta_stage1"),
               theta_type=_read_df("theta_type"), tumor_purity=_read_df("tumor_purity"),
               genes=_read_str("genes"), Z=Z, states=_read_str("states"), meta=meta)
    if mmap:
        out["_h5file"] = f
    else:
        f.close()
    return out