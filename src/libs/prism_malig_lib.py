#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/07/31
# Udated  on 2026/07/31
# @author: Flavio Lichtenstein
# @local: Home sweet home

"""
Claude - Opus 5 High code

Sample-level clustering of the *malignant compartment* recovered by BayesPrism,
with perturbation anchoring against Tahoe-100M and (optionally) scGPT.

Pipeline
--------
    init()                                  # -> CellNameState
        |
    prepare_malignant_matrix()              # compartment CPM, share filter,
        |                                   # purity decoupling, HVG - high-variable genes + keep_genes
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

import hashlib
import json
import os
import warnings
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from pathlib import Path

from pymupdf import canon
from scipy.cluster.hierarchy import fcluster, linkage, cophenet
from scipy.spatial.distance import pdist, squareform
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from libs.Basic import pdreadcsv, pdwritecsv, create_dir, title_replace


__version__ = "0.35.0"          # bump on every edit; check with pml.__version__

__all__ = ["CellNameState", "MalignantCluster", "__version__"]


# ===========================================================================
# 0. Container
# ===========================================================================
@dataclass
class CellNameState:
    """Everything pulled out of a BayesPrism `res` for one compartment."""

    Z: pd.DataFrame                       # samples x genes, compartment counts
    df_theta: pd.DataFrame                   # samples x cell types
    cell_name: str = "malignant"
    Z_all: Optional[Dict[str, pd.DataFrame]] = None   # ct -> samples x genes
    share: Optional[pd.DataFrame] = None  # samples x genes, malignant fraction
    meta: Dict = field(default_factory=dict)

    @property
    def filt_cell_name_theta(self) -> pd.Series:
        return self.df_theta[self.cell_name]

    def compute_share(self) -> pd.DataFrame:
        """Per-sample, per-gene fraction of total expression assigned to the
        malignant compartment. Requires Z_all."""
        if self.Z_all is None:
            raise ValueError(
                "Z_all is required. Re-run prism.open_bayesprism() so that "
                "Z_all/ is populated, or supply `share=` directly (e.g. from "
                "PRISM.gene_compartment_share_nnls)."
            )

        genes = self.Z.columns

        if isinstance(self.Z_all, dict):
            total = np.zeros(self.Z.shape, dtype=float)
            for ct, Zct in self.Z_all.items():
                total += Zct.reindex(index=self.Z.index, columns=genes).values
        else:
            raise TypeError(
                "Z_all must be a dict {cell_type: samples x genes}. A raw 3-D "
                "array from full_Z should go through MalignantCluster.build_CellNameState_from_full_Z, "
                "which resolves the axis order and computes `share` itself."
            )

        with np.errstate(invalid="ignore", divide="ignore"):
            sh = np.where(total > 0, self.Z.values / total, np.nan)

        self.share = pd.DataFrame(sh, index=self.Z.index, columns=genes)

        return self.share


class MalignantCluster:
    """Driver: BayesPrism `res` -> malignant-compartment sample clusters.

    Holds the deconvolution inputs so the pipeline can be re-run with
    different filters without re-deconvolving.
    """

    def __init__(self, prism, res, df_bulk, ref,
                 root_mprog_disease: Path,
                 organ: str = "Pancreas",
                 mal_cell_name: str = "Ductal cell type 2",
                 cell_types: Optional[Sequence[str]] = None):

        self.prism, self.res, self.df_bulk, self.ref = prism, res, df_bulk, ref

        self.root_mprog_disease = root_mprog_disease
        self.root_cluster = create_dir(root_mprog_disease / "cluster")
        self.root_tahoe   = create_dir(root_mprog_disease / "tahoe")
        self.root_prism   = create_dir(root_mprog_disease / "deconv")
        self.root_lfc     = create_dir(root_mprog_disease / "lfc")

        self.__lib_version__ = __version__
        self.mal_cell_name = mal_cell_name
        self.organ = organ
        self.df_theta = res.theta

        if self.mal_cell_name not in self.df_theta.columns:
            raise KeyError(
                f"'{self.mal_cell_name}' not in df_theta columns {list(self.df_theta.columns)}. "
            )
        '''
            Pancreas: particular case, Peng-2019 (CRA001160) reference the malignant label is often 'Ductal cell type 2'
        '''

        self.program1_panel = ("FAM83A-AS1", "TFAP2A-AS2", "HOXA10-AS",
                               "HOXB-AS3", "HOXB-AS4", "LEMD1-AS1", "MIR7-3HG")
        self._TAHOE_REPO = "tahoebio/Tahoe-100M"

        self.program_gene_list = []

        self.Z_full, self.genes_full = prism.full_Z(res, df_bulk, ref)
        self.Z_mal = prism.state_expression(self.Z_full, self.genes_full,
                                            res, self.mal_cell_name)
        self.cell_types = (list(cell_types) if cell_types is not None
                           else list(self.df_theta.columns))

        # Route through build_CellNameState_from_full_Z so axis resolution, the state_expression
        # concordance check and recovered-gene bookkeeping actually run.
        self.cns = self.build_CellNameState_from_full_Z(
            Z_full = self.Z_full, 
            genes_full = self.genes_full,
            cell_types=self.cell_types,
            cell_name=self.mal_cell_name,
            samples=self.df_theta.index,
            df_theta=self.df_theta,
            Z_mal=self.Z_mal,
            recovered_genes=self.program1_panel,
        )
        self.Z = self.cns.Z

    # ===========================================================================
    # 1. run --> end-to-end convenience
    # ===========================================================================
    def run(
        self, 
        ks: Iterable[int] = range(2, 9),
        tahoe: bool = False,
        tahoe_kw: Optional[Dict] = None,
        **kw,
    ):
        """Cluster the malignant compartment; optionally anchor against Tahoe.

        `**kw` goes to prepare_malignant_matrix; `tahoe_kw` to load_tahoe_de
        (e.g. dict(mode="download", hf_token=..., threads=2)).
        """
        X, diag = self.prepare_malignant_matrix(keep_genes=self.program1_panel, **kw)
        cc = self.consensus_cluster(X, ks=ks)
        k = self.choose_k(cc)
        labels = cc[k]["labels"]
        sig = self.cluster_signatures(X, labels)

        out = {"X": X, "diag": diag, "consensus": cc, "k": k,
               "labels": labels, "signatures": sig}

        if tahoe:
            # genes=X.columns is not optional: without it the pancreas subset
            # is ~5e8 rows. Callers may override, but never silently drop it.
            tkw = {"organs": (self.organ,), "genes": X.columns,
                   "save_to": self.tahoe_derived_dir()}
            tkw.update(tahoe_kw or {})
            if tkw.get("genes") is None:
                warnings.warn(
                    "run(tahoe=True) with genes=None will pull the full gene "
                    "universe and is likely to exhaust memory.")
            S, cond = self.load_tahoe_de(**tkw)
            out["tahoe"] = self.score_clusters_vs_tahoe(sig, S, cond)
            out["tahoe_S"], out["tahoe_cond"] = S, cond
        
        return out


    def build_CellNameState_from_full_Z(
        self,
        Z_full,
        genes_full: Sequence[str],
        cell_types: Sequence[str],
        cell_name: str = "Ductal cell type 2",
        samples: Optional[Sequence[str]] = None,
        df_theta: Optional[pd.DataFrame] = None,
        Z_mal: Optional[pd.DataFrame] = None,
        recovered_genes: Optional[Sequence[str]] = None,
    ) -> CellNameState:
        """
        Build a CellNameState straight from `prism.full_Z(res, df_bulk, ref)`.

        Preferred over `load_bayesprism_export` -- no R round-trip, and the
        per-gene compartment share is exact rather than reconstructed.

            Z_full, genes_full = prism.full_Z(res, df_bulk, ref)
            Zmal = prism.state_expression(Z_full, genes_full, res, "Ductal cell type 2")
            cns = build_CellNameState_from_full_Z(Z_full, genes_full,
                            cell_types=list(s2t.values()),      # or ref cell types
                            cell_name=cell_name,
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
        if isinstance(Z_full, dict):
            cell_types = list(Z_full.keys())
            first = next(iter(Z_full.values()))
            if samples is None:
                samples = list(first.index)
            A = np.stack([np.asarray(pd.DataFrame(Z_full[c])
                                     .reindex(index=samples, columns=genes_full)
                                     .values, dtype=float)
                          for c in cell_types], axis=2)
        else:
            A = np.asarray(Z_full)
            if A.dtype.kind != "f":            # keep f32 as f32
                A = A.astype(np.float32)
        if A.ndim != 3:
            raise ValueError(
                f"expected a 3-D (sample, gene, cell-type) array or a dict of "
                f"2-D matrices from full_Z, got shape {getattr(A, 'shape', None)}."
            )

        genes, cell_type_list = list(genes_full), list(cell_types)
        ng, nct = len(genes), len(cell_type_list)
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
        # samples x genes x ct
        A = np.moveaxis(A, (sa, ga, ca), (0, 1, 2))          

        idx = pd.Index(samples) if samples is not None else pd.RangeIndex(A.shape[0])
        if cell_name not in cell_type_list:
            raise KeyError(f"'{cell_name}' not in cell_types {cell_type_list}")

        ci = cell_type_list.index(cell_name)

        Zc = pd.DataFrame(np.ascontiguousarray(A[:, :, ci]),
                          index=idx, columns=genes)
        # Accumulate the denominator in f64, then divide in place and store the
        # share as f32: at 300 x 60k x 10 the naive version holds three extra
        # f64 copies of a samples x genes plane.
        total = A.sum(axis=2, dtype=np.float64)
        with np.errstate(invalid="ignore", divide="ignore"):
            zero = total == 0
            np.divide(A[:, :, ci], total, out=total, where=~zero)
            total[zero] = np.nan
        share = pd.DataFrame(total.astype(np.float32), index=idx, columns=genes)
        del total, zero

        if df_theta is None:
            m = A.sum(axis=1, dtype=np.float64)
            df_theta = pd.DataFrame(m / m.sum(axis=1, keepdims=True),
                                index=idx, columns=cell_type_list)
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
            Zm = pd.DataFrame(Z_mal)
            # state_expression may hand back genes x samples; detect and fix,
            # otherwise the reindex below yields an all-NaN frame and every
            # gene silently fails the share filter.
            hit_rows = Zm.index.isin(idx).mean()
            hit_cols = Zm.columns.isin(idx).mean()
            if hit_cols > hit_rows:
                warnings.warn(
                    f"Z_mal looks transposed (index matches samples at "
                    f"{hit_rows:.0%}, columns at {hit_cols:.0%}); transposing.")
                Zm = Zm.T
            # An empty frame defeats the transpose heuristic above: with zero
            # rows both hit rates are 0.0, so nothing is transposed and the
            # error below reports an empty index instead of the real cause.
            if 0 in Zm.shape:
                raise ValueError(
                    f"Z_mal is empty (shape {Zm.shape}); state_expression() "
                    "returned nothing, so the failure is upstream of this "
                    "method. Most likely the gene vocabularies disagree: after "
                    "the ensembl switch df_bulk is ENSG-indexed, so full_Z() "
                    "must be called with an ENSG-indexed reference too. "
                    f"Check len(genes_full)={len(genes)}, "
                    f"Z_full.shape={getattr(A, 'shape', None)}, and that "
                    "res.theta.index is non-empty.")
            if not Zm.index.isin(idx).any():
                raise ValueError(
                    "Z_mal shares no index labels with the sample index. "
                    f"Z_mal index[:3]={list(Zm.index[:3])}, "
                    f"expected sample labels like {list(idx[:3])}.")
            Zm = Zm.reindex(index=idx)
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

        return CellNameState(Z=Zc, df_theta=df_theta, cell_name=cell_name, share=share, meta=meta)


    # ===========================================================================
    # 2. Matrix preparation
    # ===========================================================================
    def prepare_malignant_matrix(
        self,
        min_share: float = 0.50,
        share_frac_samples: float = 0.50,
        min_counts: float = 10.0,
        min_theta: float = 0.05,
        keep_samples: Optional[Sequence[str]] = None,
        drop_pattern: Optional[str] = None,
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
            Drop samples whose malignant fraction is below this -- their Z is almost entirely prior.
        decouple_purity
            Regress log-CPM on theta_mal (+ extra_covariates) per gene and keep the residuals.
        keep_genes
            Force-retain genes regardless of filters (e.g. your Program-1 panel:
            FAM83A-AS1, HOXA10-AS, HOXB-AS3, HOXB-AS4, MIR7-3HG, TFAP2A-AS2,
            LEMD1-AS1). They are still reported in diagnostics so you can see
            whether they *would* have survived.

        Explanation: Purity is the default answer, not the interesting one. get.exp() returns counts 
        attributable to the malignant compartment, so row sums scale with theta_mal × depth. 
        Clustering that matrix recovers tumour cellularity — the exact Program-1/Program-2 axis 
        you already suspect is confounded. The module CPMs within the compartment, 
        then optionally regresses each gene on theta_mal. diag["pc_theta_pearson"] is the check: 
        run with --no-decouple too, and if structure only survives there, it's purity.

        """
        diag: Dict = {}

        # --- sample filter ------------------------------------------------------
        # Normal / adjacent-tissue samples have no malignant compartment, so
        # their Z_mal is close to pure BayesPrism prior -- the reference
        # profile, not the sample. Left in, they dominate the first split and
        # the "subtypes" are really tumour-vs-normal.
        _idx = self.cns.Z.index
        _sel = pd.Series(True, index=_idx)

        if keep_samples is not None:
            _want = set(map(str, keep_samples))
            _sel &= pd.Series(_idx.astype(str).isin(_want), index=_idx)
            _unk = _want - set(_idx.astype(str))
            if _unk:
                warnings.warn(f"keep_samples not in Z index: {sorted(_unk)[:5]}")

        if drop_pattern:
            _hit = pd.Series(_idx.astype(str).str.contains(drop_pattern, regex=True),
                             index=_idx)
            _sel &= ~_hit

        _n_excl = int((~_sel).sum())
        if _n_excl:
            print(f"excluded {_n_excl}/{len(_idx)} samples by keep_samples/drop_pattern")

        diag["samples_excluded_by_filter"] = list(_idx[~_sel])

        keep_s = _sel & (self.cns.filt_cell_name_theta.reindex(_idx) >= min_theta)
        if (~keep_s).any():
            warnings.warn(
                f"dropping {int((~keep_s).sum())} samples with theta_{self.cns.cell_name}"
                f" < {min_theta}: {list(self.cns.filt_cell_name_theta.index[~keep_s])}"
            )
        Z = self.cns.Z.loc[keep_s]
        theta_mal = self.cns.filt_cell_name_theta.loc[keep_s]
        diag["samples_dropped"] = list(self.cns.filt_cell_name_theta.index[~keep_s])

        # --- compartment CPM ----------------------------------------------------
        lib = Z.sum(axis=1)
        cpm = Z.div(lib, axis=0) * 1e6
        logx = np.log2(cpm + 1.0)

        # --- gene filters -------------------------------------------------------
        expressed = (Z.median(axis=0) >= min_counts)
        diag["n_genes_expressed"] = int(expressed.sum())

        if self.cns.share is not None:
            sh = self.cns.share.reindex(index=Z.index, columns=Z.columns)
            computable = sh.notna().any(axis=0)
            share_ok = (sh >= min_share).mean(axis=0) >= share_frac_samples
            diag["n_genes_share_not_computable"] = int((~computable).sum())
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
            # "no_share" != "failed": recovered loci absent from Z_full have
            # no computable share and would be dropped for lack of data, not
            # for lack of compartment specificity.
            _sh_ok = locals().get("computable")
            diag["forced_genes_status"] = {
                g: ("passed" if bool(keep_g.get(g, False))
                    else ("no_share" if (_sh_ok is not None
                                         and not bool(_sh_ok.get(g, False)))
                          else "failed"))
                for g in forced
            }
            keep_g.loc[forced] = True
        logx = logx.loc[:, keep_g]
        diag["n_genes_kept"] = int(logx.shape[1])

        if logx.shape[1] < 2:
            raise ValueError(
                "gene filters left {n} gene(s). Cascade:\n"
                "  total genes in Z            : {tot}\n"
                "  median count >= {mc:<10g}  : {ex}\n"
                "  share computable            : {sc}\n"
                "  >= {cns:g} share in >= {sf:.0%} samples : {so}\n"
                "  intersection (kept)         : {n}\n"
                "Most common causes: (a) Z is not on a count scale, so "
                "min_counts={mc:g} removes everything -- check "
                "Z.median().median()={zmed:.3g}; (b) min_share={cns:g} is too "
                "strict for a gene_subset fit, where the share denominator "
                "covers only the fitted compartments. Call "
                "diagnose_filters() to sweep thresholds."
                .format(n=int(logx.shape[1]), tot=int(len(keep_g)),
                        mc=min_counts, ex=int(expressed.sum()),
                        sc=diag.get("n_genes_share_not_computable", "n/a") if
                           isinstance(diag.get("n_genes_share_not_computable"), str)
                           else int(len(keep_g)) - int(diag.get("n_genes_share_not_computable", 0)),
                        ms=min_share, sf=share_frac_samples,
                        so=int(share_ok.sum()),
                        zmed=float(np.nanmedian(Z.values))))

        # --- purity decoupling --------------------------------------------------
        if decouple_purity:
            C = pd.DataFrame({"theta_mal": theta_mal, "log_lib": np.log(lib)})
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
            sel_set = forced_set | set(sel)
            Xc = Xc.loc[:, [g for g in Xc.columns if g in sel_set]]
        diag["n_hvg"] = int(Xc.shape[1])

        # --- purity leakage check ----------------------------------------------
        # Correlate PCs with theta on BOTH matrices. On the decoupled matrix
        # this is ~0 by construction (residuals are orthogonal to the
        # covariates), so it proves nothing -- the informative number is the
        # pre-decoupling one.
        def _pc_theta(M: pd.DataFrame, name: str):
            # PCA needs >=2 samples AND >=2 features; theta needs to vary or
            # the correlation is undefined (0/0 -> nan).
            if M.shape[0] < 2 or M.shape[1] < 2:
                return f"unavailable - {name}.shape={M.shape}, need >=2 x >=2"
            if float(theta_mal.std()) == 0:
                return "unavailable - theta_mal is constant across samples"
            n_comp = min(5, min(M.shape) - 1)
            if n_comp < 1:
                return f"unavailable - {name} supports 0 components"
            pcs = PCA(n_components=n_comp).fit_transform((M - M.mean()).values)
            return [float(np.corrcoef(pcs[:, i], theta_mal.values)[0, 1])
                    for i in range(pcs.shape[1])]

        diag["pc_theta_pearson_raw"] = _pc_theta(logx, "logx")

        diag["pc_theta_pearson"] = _pc_theta(Xc, "Xc")
        diag["decouple_purity"] = bool(decouple_purity)
        if decouple_purity:
            diag["pc_theta_note"] = (
                "pc_theta_pearson is ~0 by construction under "
                "decouple_purity=True (residuals are orthogonal to the "
                "covariates); read pc_theta_pearson_raw instead.")
        else:
            diag["pc_theta_note"] = (
                "decouple_purity=False, so pc_theta_pearson and "
                "pc_theta_pearson_raw are the same matrix and both are "
                "informative: a large |r| on an early PC means the clustering "
                "is tracking tumour purity.")

        # Global signal level per sample -- catches samples whose malignant
        # compartment is near-empty or degraded. These cluster out on their own
        # and look like a subtype until you notice every program is down.
        diag["sample_mean_expr"] = Xc.mean(axis=1).round(3)
        diag["sample_total_Z"] = self.cns.Z.loc[Xc.index].sum(axis=1)

        diag["theta_mal"] = theta_mal
        diag["n_samples_used"] = int(len(theta_mal))
        self.samples_used = list(theta_mal.index)
        if _n_excl:
            diag["theta_excluded"] = self.cns.filt_cell_name_theta[~_sel].describe().round(4)
            diag["theta_kept"] = self.cns.filt_cell_name_theta[_sel].describe().round(4)

        return (Xc, diag) if return_diagnostics else (Xc, {})


    def diagnose_filters(
        self,
        min_shares: Sequence[float] = (0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
        min_counts_grid: Sequence[float] = (0.0, 1.0, 5.0, 10.0, 50.0),
        share_frac_samples: float = 0.50,
        min_theta: float = 0.05,
    ) -> Dict:
        """Sweep the gene filters and report how many genes survive each.

        Run this before touching the thresholds. It answers the two questions
        that matter: is Z on a count scale (does min_counts bite at all?), and
        is min_share calibrated for this fit (a `gene_subset` deconvolution has
        a restricted share denominator, so shares run high and 0.5 may be
        either trivially permissive or, if compartments are collinear, fatal).
        """
        keep_s = self.cns.filt_cell_name_theta >= min_theta
        Z = self.cns.Z.loc[keep_s]

        out = {
            "n_samples_kept": int(keep_s.sum()),
            "n_samples_dropped": int((~keep_s).sum()),
            "n_genes_total": int(Z.shape[1]),
            "Z_median_of_medians": float(np.nanmedian(Z.median(axis=0))),
            "Z_looks_like_counts": bool(np.nanmedian(Z.values) > 1.0),
        }

        med = Z.median(axis=0)
        out["by_min_counts"] = pd.Series(
            {c: int((med >= c).sum()) for c in min_counts_grid},
            name="n_genes").rename_axis("min_counts")

        if self.cns.share is None:
            out["share"] = "unavailable -- shrinkage guard would be OFF"
            return out

        sh = self.cns.share.reindex(index=Z.index, columns=Z.columns)
        out["n_genes_share_not_computable"] = int((~sh.notna().any(axis=0)).sum())
        out["share_quantiles"] = sh.stack().quantile(
            [0.05, 0.25, 0.5, 0.75, 0.95]).round(3)
        out["by_min_share"] = pd.Series(
            {m: int((((sh >= m).mean(axis=0) >= share_frac_samples)).sum())
             for m in min_shares},
            name="n_genes").rename_axis("min_share")

        grid = pd.DataFrame(
            {m: {c: int((((sh >= m).mean(axis=0) >= share_frac_samples)
                         & (med >= c)).sum())
                 for c in min_counts_grid} for m in min_shares})
        grid.index.name = "min_counts"; grid.columns.name = "min_share"
        out["joint_grid"] = grid
        return out


    # ===========================================================================
    # 3. Consensus clustering (Monti et al.)
    # ===========================================================================
    @staticmethod
    def _pac(consensus: np.ndarray, lo: float = 0.1, hi: float = 0.9) -> float:
        # Accept a DataFrame too: cc[k]["consensus"] is one, and label-based
        # __getitem__ would misread the positional triu indices.
        consensus = np.asarray(consensus)
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
        if X.shape[1] < 2:
            raise ValueError(
                f"X has {X.shape[1]} feature(s); nothing to cluster on. This "
                "almost always means prepare_malignant_matrix() filtered out "
                "every gene -- run diagnose_filters() to see which threshold.")
        
        if X.shape[0] < 3:
            raise ValueError(f"X has {X.shape[0]} sample(s); need >= 3.")

        rng = np.random.default_rng(random_state)
        n = X.shape[0]
        Xv = X.values
        out: Dict[int, Dict] = {}

        ks = [k for k in ks if k < n]
        if not ks:
            raise ValueError(f"no k in the requested range is < n_samples={n}.")

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
                "pac": self._pac(C),
                "cophenetic": float(coph),
                "silhouette": sil,
                "linkage": Zl,
            }
        return out


    def choose_k(self, cc: Dict[int, Dict], pac_tol: float = 0.02,
                 min_cluster_frac: float = 0.0) -> int:
        """Smallest k whose PAC is within `pac_tol` of the minimum.

        `min_cluster_frac` rejects degenerate solutions first. PAC rewards
        reproducibility, and peeling 3 outliers off 108 samples is maximally
        reproducible -- every resample isolates them identically -- so PAC
        approaches 0 while the split carries no biology. Setting
        min_cluster_frac=0.10 requires every cluster to hold >=10% of samples
        before a k is eligible.
        """
        pacs = {k: v["pac"] for k, v in cc.items()}

        if min_cluster_frac > 0:
            ok = {}
            for k, v in cc.items():
                counts = v["labels"].value_counts()
                if counts.min() / counts.sum() >= min_cluster_frac:
                    ok[k] = pacs[k]
            if not ok:
                warnings.warn(
                    f"no k has all clusters >= {min_cluster_frac:.0%} of "
                    f"samples; smallest-cluster fractions were "
                    f"{ {k: round(v['labels'].value_counts().min()/len(v['labels']),3) for k, v in cc.items()} }. "
                    "This usually means the data has a quality gradient rather "
                    "than discrete groups -- fix that before choosing k.")
            else:
                pacs = ok

        best = min(pacs.values())
        return min(k for k, p in pacs.items() if p <= best + pac_tol)

    def cluster_summary(self, cc: Dict[int, Dict]) -> pd.DataFrame:
        """PAC / cophenetic / silhouette AND cluster sizes for every k.

        Always read sizes alongside PAC: a low PAC with a 3-vs-108 split is an
        outlier detector, not a subtyping.
        """
        rows = {}
        for k, v in cc.items():
            n = v["labels"].value_counts()
            rows[k] = {"pac": round(v["pac"], 4),
                       "cophenetic": round(v["cophenetic"], 4),
                       "silhouette": round(v["silhouette"], 4)
                       if v["silhouette"] == v["silhouette"] else np.nan,
                       "sizes": sorted(n.tolist()),
                       "min_frac": round(n.min() / n.sum(), 3)}
        return pd.DataFrame(rows).T

    @staticmethod
    def _bh(p: pd.Series) -> pd.Series:
        """Benjamini-Hochberg FDR."""
        p = p.fillna(1.0)
        n = len(p)
        order = np.argsort(p.values)
        ranked = p.values[order]
        q = ranked * n / (np.arange(n) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        out = np.empty(n)
        out[order] = np.clip(q, 0, 1)
        return pd.Series(out, index=p.index)

    def signature_table(
        self,
        X: pd.DataFrame,
        labels: pd.Series,
        heldout: Optional[Dict[int, pd.DataFrame]] = None,
        min_abs_lfc: float = 0.0,
        sort_by: str = "stat",
    ) -> pd.DataFrame:
        """Tidy per-cluster DE table: lfc, moderated stat, p, FDR.

        One row per (cluster, gene). Columns:
          lfc            difference of means on the clustering scale (log2 CPM,
                         purity-residualised if decouple_purity was on -- so it
                         is a log2 fold-change ONLY if that scale was log2 CPM)
          stat           SAM-style moderated t used for the WTCS ranking
          p_nominal      Welch t p-value, circular (see below)
          fdr_nominal    BH across genes within each cluster
          frac_splits_sig / median_p / fdr_of_median  if `heldout` is supplied

        On k=2: one-vs-rest is a single contrast, so cluster 1 and cluster 2
        rows are mirror images -- identical p, sign-flipped lfc. Report one.

        p_nominal/fdr_nominal are anticonservative because the clusters were
        chosen to separate these samples. Measured on null data: ~1-8% of genes
        pass fdr_nominal<0.05 with no structure present, worsening as n falls.
        Use them to rank; use `heldout` columns to claim.
        """
        sig = self.cluster_signatures(X, labels, n_top=X.shape[1],
                                      min_abs_lfc=min_abs_lfc, pvalues=True)
        rows = []
        for c, v in sig.items():
            df = pd.DataFrame({
                "cluster": c,
                "gene": v["stat"].index,
                "lfc": v["lfc"].values,
                "stat": v["stat"].values,
                "p_nominal": v["p_nominal"].values,
                "fdr_nominal": v["fdr_nominal"].values,
                "sd_cluster": v["sd_cluster"].values,
                "sd_rest": v["sd_rest"].values,
                "welch_df": v["welch_df"].values,
                "n_cluster": v["n"],
            })
            if heldout is not None and c in heldout:
                h = heldout[c]
                df = df.merge(
                    h[["median_p", "fdr_of_median", "frac_splits_sig"]],
                    left_on="gene", right_index=True, how="left")
            rows.append(df)

        out = pd.concat(rows, ignore_index=True)
        out["direction"] = np.where(out["lfc"] > 0, "up", "down")

        # Smaller-is-better for p-values/FDR, larger-is-better for effect
        # sizes; sorting everything descending puts the best gene last and
        # makes `rank` read backwards.
        _ascending = sort_by in {"p_nominal", "fdr_nominal",
                                 "median_p", "fdr_of_median"}
        _key = out[sort_by].abs() if sort_by in {"lfc", "stat", "mean_lfc"} \
            else out[sort_by]
        out = (out.assign(_k=_key)
                  .sort_values(["cluster", "_k"], ascending=[True, _ascending])
                  .drop(columns="_k").reset_index(drop=True))
        out["rank"] = out.groupby("cluster").cumcount() + 1
        return out

    def heldout_signatures(
        self,
        X: pd.DataFrame,
        k: int,
        n_splits: int = 50,
        frac: float = 0.5,
        reference_labels: Optional[pd.Series] = None,
        alpha: float = 0.05,
        random_state: int = 0,
        **cc_kw,
    ) -> Dict[int, pd.DataFrame]:
        """Valid p-values by sample splitting: cluster on A, test on B.

        Each split clusters a random `frac` of samples (A), assigns the held-out
        samples (B) to those clusters by nearest centroid, then runs a Welch
        t-test **within B only**. Because B played no part in defining the
        clusters, its p-values are not circular.

        Clusters are matched to `reference_labels` (your full-data solution) by
        majority overlap on A, so results are comparable across splits.

        Returns {cluster: DataFrame(median_p, fdr_of_median, frac_splits_sig,
        mean_lfc, n_splits_tested)}. `frac_splits_sig` is the useful column: a
        gene significant in 45/50 splits is far more trustworthy than one with a
        small p in a single split.

        Caveat: each split halves n, so power drops sharply. Genes that fail
        here are not shown to be null -- they are underpowered. Absence of
        evidence, in a design that deliberately spends evidence on validity.
        """
        rng = np.random.default_rng(random_state)
        samples = np.array(X.index)
        n = len(samples)
        acc: Dict[int, Dict[str, list]] = {}

        for rep in range(n_splits):
            idx = rng.permutation(n)
            nA = max(k + 1, int(round(frac * n)))
            A, B = samples[idx[:nA]], samples[idx[nA:]]
            if len(B) < 2 * k:
                continue

            cc = self.consensus_cluster(X.loc[A], ks=[k],
                                        random_state=int(rng.integers(1e6)),
                                        **cc_kw)
            labA = cc[k]["labels"]

            cent = X.loc[A].groupby(labA).mean()
            # assign B by max correlation to a cluster centroid
            Bc = X.loc[B].sub(X.loc[B].mean(axis=1), axis=0)
            Cc = cent.sub(cent.mean(axis=1), axis=0)
            corr = (Bc.values @ Cc.values.T) / (
                np.linalg.norm(Bc.values, axis=1)[:, None]
                * np.linalg.norm(Cc.values, axis=1)[None, :] + 1e-12)
            labB = pd.Series(cent.index[np.argmax(corr, axis=1)], index=B)

            # relabel A-clusters to reference ids by majority overlap
            if reference_labels is not None:
                ct = pd.crosstab(labA, reference_labels.reindex(A))
                mapping = ct.idxmax(axis=1).to_dict()
                labB = labB.map(mapping)

            for c in sorted(labB.dropna().unique()):
                g1, g2 = X.loc[labB[labB == c].index], X.loc[labB[labB != c].index]
                if len(g1) < 2 or len(g2) < 2:
                    continue
                d = g1.mean(axis=0) - g2.mean(axis=0)
                v1, v2 = g1.var(axis=0, ddof=1), g2.var(axis=0, ddof=1)
                n1, n2 = len(g1), len(g2)
                se = np.sqrt(v1 / n1 + v2 / n2)
                with np.errstate(invalid="ignore", divide="ignore"):
                    tw = d / se
                    dfree = ((v1 / n1 + v2 / n2) ** 2 /
                             ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)))
                pv = pd.Series(2 * stats.t.sf(np.abs(tw.values), dfree.values),
                               index=X.columns).fillna(1.0)
                a_ = acc.setdefault(int(c), {"p": [], "lfc": []})
                a_["p"].append(pv)
                a_["lfc"].append(d)

        out = {}
        for c, v in acc.items():
            P = pd.concat(v["p"], axis=1)
            L = pd.concat(v["lfc"], axis=1)
            med = P.median(axis=1)
            sig_frac = P.apply(lambda col: self._bh(col) < alpha).mean(axis=1)
            out[c] = pd.DataFrame({
                "median_p": med,
                "fdr_of_median": self._bh(med),
                "frac_splits_sig": sig_frac,
                "mean_lfc": L.mean(axis=1),
                "n_splits_tested": P.shape[1],
            }).sort_values("frac_splits_sig", ascending=False)
        return out

    def cluster_signatures(
        self,
        X: pd.DataFrame,
        labels: pd.Series,
        n_top: int = 150,
        min_abs_lfc: float = 0.0,
        pvalues: bool = False,
    ) -> Dict[int, Dict[str, object]]:
        """One-vs-rest signed signatures with a simple moderated t.

        `pvalues=True` adds `p_nominal` / `fdr_nominal`. THESE ARE NOT VALID
        INFERENCE. The clusters were chosen to maximise separation in this same
        X, so testing that separation here is circular ("double dipping"): under
        a null of one homogeneous group, any clustering still yields tiny
        p-values. Use them to rank genes within a cluster, never to claim a
        cluster is real or that a gene is "significantly" associated.

        For p-values you can defend, use heldout_signatures(), which clusters on
        one half of the samples and tests on the other.

        Returns {cluster: {"stat": Series (all genes, signed),
                        "up": [genes], "down": [genes]}}.
        """
        labels = labels.reindex(X.index)
        sig: Dict[int, Dict[str, object]] = {}

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

            entry_p = {}
            if pvalues:
                # Welch t WITHOUT the s0 offset: with s0 > 0 the statistic is
                # not t-distributed, so p would be wrong even ignoring the
                # circularity below.
                na, nb = len(a), len(b)
                # Genes with zero variance in the SMALLER group are the trap:
                # Welch df collapses from ~na-1 to ~nb-1 and se shrinks to
                # sqrt(vb/nb), so t explodes and p reaches 1e-40 from a group
                # of 6. Floor each group's variance at a small positive
                # quantile so a single degenerate gene cannot dominate.
                _floor_a = float(np.nanpercentile(va[va > 0], 5)) if (va > 0).any() else 0.0
                _floor_b = float(np.nanpercentile(vb[vb > 0], 5)) if (vb > 0).any() else 0.0
                va_f, vb_f = va.clip(lower=_floor_a), vb.clip(lower=_floor_b)
                se_f = np.sqrt(va_f / na + vb_f / nb)
                with np.errstate(invalid="ignore", divide="ignore"):
                    t_w = d / se_f
                    dfree = ((va_f / na + vb_f / nb) ** 2 /
                             ((va_f / na) ** 2 / (na - 1)
                              + (vb_f / nb) ** 2 / (nb - 1)))
                pv = pd.Series(2 * stats.t.sf(np.abs(t_w.values), dfree.values),
                               index=t.index).fillna(1.0)
                entry_p = {
                    "p_nominal": pv,
                    "fdr_nominal": self._bh(pv),
                    "sd_cluster": np.sqrt(va),
                    "sd_rest": np.sqrt(vb),
                    "welch_df": pd.Series(dfree.values, index=t.index),
                    "n_zero_var_in_cluster": int((va == 0).sum()),
                }

            sig[int(c)] = {
                "stat": t,
                "lfc": d,
                **entry_p,
                "up": [g for g, v in up.items() if v > 0],
                "down": [g for g, v in dn.items() if v < 0],
                "n": int((labels == c).sum()),
            }
        return sig



    # ===========================================================================
    # 4b. The other axis: TME composition (see compare_axes)
    # ===========================================================================
    def tme_axis(self, exclude: Optional[Sequence[str]] = None,
                 ks: Iterable[int] = range(2, 7), pseudo: float = 1e-4,
                 **cc_kw) -> Dict:
        """Cluster samples on TME *composition*, orthogonalised against purity.

        Fractions are renormalised among non-malignant compartments only, so the
        axis is "what is the stroma made of", not "how much stroma is there".
        CLR before distance: the sum-to-one constraint on compositional data
        induces structure that Euclidean/correlation distance will happily
        cluster on.
        """
        drop = {self.mal_cell_name} | set(exclude or [])
        keep = [c for c in self.df_theta.columns if c not in drop]
        if len(keep) < 2:
            raise ValueError(f"need >=2 non-malignant compartments, have {keep}")

        T = self.df_theta[keep]
        if getattr(self, "samples_used", None):
            T = T.reindex(self.samples_used)   # align with clustered set
        T = T.div(T.sum(axis=1), axis=0)
        L = np.log(T.values + pseudo)
        clr = pd.DataFrame(L - L.mean(axis=1, keepdims=True),
                           index=T.index, columns=keep)

        cc = self.consensus_cluster(clr, ks=ks, **cc_kw)
        k = self.choose_k(cc)
        return {"clr": clr, "consensus": cc, "k": k,
                "labels": cc[k]["labels"], "compartments": keep}

    @staticmethod
    def compare_axes(a: pd.Series, b: pd.Series,
                     names: Tuple[str, str] = ("malignant", "TME")) -> Dict:
        """Cross-tabulate two label sets and quantify their dependence.

        Low ARI => the axes carry different information, report the joint
        (tumour x stroma) label. High ARI => check whether the malignant
        clustering simply recovered purity.
        """
        from scipy.stats import chi2_contingency
        from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

        idx = a.index.intersection(b.index)
        a, b = a.loc[idx], b.loc[idx]
        ct = pd.crosstab(a.rename(names[0]), b.rename(names[1]))
        chi2, pval, dof, _ = chi2_contingency(ct)
        n = ct.values.sum()
        v = (float(np.sqrt(chi2 / (n * (min(ct.shape) - 1))))
             if min(ct.shape) > 1 else np.nan)
        return {"crosstab": ct, "chi2": float(chi2), "p": float(pval),
                "cramers_v": v, "ari": float(adjusted_rand_score(a, b)),
                "ami": float(adjusted_mutual_info_score(a, b))}

    # ===========================================================================
    # 5. Tahoe-100M perturbation reference
    # ===========================================================================

    # --- where Tahoe lives on disk -------------------------------------------
    _TAHOE_SMALL = ["metadata/cell_line_metadata.parquet",
                    "metadata/drug_metadata.parquet",
                    "metadata/gene_metadata.parquet",
                    "metadata/sample_metadata.parquet"]
    _TAHOE_DE_GLOB = "metadata/pseudobulk_differential_expression/*.parquet"


    def tahoe_derived_dir(self, save_to=None) -> Path:
        """Where derived S/cond parquet is cached. Defaults to
        `<root_tahoe>/derived` so every call site shares it."""
        d = Path(save_to) if save_to not in (None, True) \
            else self.root_tahoe / "derived"
        create_dir(d)
        return d

    def download_tahoe(self, de: bool = False,
                       force: bool = False):
        """
        Download Tahoe-100M metadata into a plain directory of your choice.
        How to load a Parquet file into a table
        https://duckdb.org/docs/current/guides/overview

        Parameters
        ----------
        de
            Also pull the pseudobulk DE table. This is the expensive part: 21
            parquet shards covering 4.09e9 rows. Leave it False and use
            `mode="stream"` unless you need repeated offline access.

        Note
        ----
        `allow_patterns=["metadata/*"]` would also drag down obs_metadata
        (101M rows) and the whole DE table. The small tables below are a few MB
        total, which is all `load_tahoe_de` needs for cell-line/drug lookup.
        """
        from huggingface_hub import snapshot_download

        patterns = list(self._TAHOE_SMALL)
        if de:
            patterns.append(self._TAHOE_DE_GLOB)

        # huggingface_hub.snapshot_download
        snapshot_download(
            repo_id=self._TAHOE_REPO,
            repo_type="dataset",
            allow_patterns=patterns,
            local_dir=self.root_tahoe,
            force_download=force,
        )
   

    # Column-name candidates for the DE table, most-specific first. The real
    # schema uses Cell_ID_Cellosaur / Cell_Name_Vevo / gene_name -- there is no
    # `cell_line_id` column, and `gene_id` (Ensembl) sits alongside `gene_name`
    # (HGNC), so ordering here decides whether a query returns rows at all.
    _DE_GENE_CANDS = ["gene_name", "gene_symbol", "symbol", "gene", "genes"]
    # DESeq2-style DE output uses log2FoldChange + a signed Wald `stat`.
    # `stat` is preferred for connectivity scoring: it is variance-stabilised,
    # so it is the closer analogue of CMap's moderated z-score, whereas raw
    # log2FC is dominated by low-count genes unless shrunk. Both are signed;
    # anything unsigned breaks the WTCS sign convention entirely.
    _DE_STAT_CANDS = ["stat", "log2foldchange", "log2fc", "log2_fold_change",
                      "logfoldchange", "logfc", "lfc", "shrunk_lfc",
                      "vision_score", "score"]
    _DE_CELL_CANDS = ["cell_id_cellosaur", "cell_line_id", "cell_name_vevo",
                      "cell_line"]
    _DE_DRUG_CANDS = ["drug", "drugname_drugconc", "treatment"]

    @staticmethod
    def _pick_col(cols, cands, override=None, label=""):
        """Resolve a column name against a schema, exact match first."""
        if override:
            if override not in cols:
                raise KeyError(f"{override!r} not in DE schema {list(cols)}")
            return override
        for c in cands:
            for col in cols:
                if col.lower() == c:
                    return col
        for c in cands:
            hits = [col for col in cols if c in col.lower()]
            if hits:
                if len(hits) > 1:
                    warnings.warn(
                        f"ambiguous {label or 'column'} resolution: {hits} all "
                        f"contain {c!r}; using {hits[0]!r}. Pass an explicit "
                        f"override if that is wrong.")
                return hits[0]
        raise KeyError(f"none of {cands} in DE schema {list(cols)}")

    def resolve_de_columns(self, cols, gene_col=None, stat_col=None) -> Dict[str, str]:
        """Which DE columns the query will actually use. Check this before a
        long scan -- a wrong `stat` column returns rows but meaningless scores,
        which is worse than returning none."""
        return {
            "gene": self._pick_col(cols, self._DE_GENE_CANDS, gene_col, "gene"),
            "stat": self._pick_col(cols, self._DE_STAT_CANDS, stat_col, "statistic"),
            "cell_line": self._pick_col(cols, self._DE_CELL_CANDS, None, "cell line"),
            "drug": self._pick_col(cols, self._DE_DRUG_CANDS, None, "drug"),
        }

    def fetch_de_shard(self, i: int = 0) -> Path:
        """Download a single DE shard to disk (cached). One file, not 1026."""
        from huggingface_hub import hf_hub_download
        root_mprog_probe = create_dir(self.root_tahoe / "_probe")
        f = hf_hub_download(repo_id=self._TAHOE_REPO, repo_type="dataset",
                            filename=self.list_de_shards()[i],
                            local_dir=root_mprog_probe)
        return Path(f)

    def inspect_de_schema(self, genes=None, shard: int = 0, n_distinct: int = 10,
                          n_rows: int = 500_000):
        """
        Report what the DE columns actually contain. Run on 0-row queries.

        Pulls ONE shard to disk and does every check in pandas. Inspecting over
        HTTP means a separate remote scan per column, which is why the previous
        version took minutes; this is one download then local work.

        Answers: does the cell-line column hold Cellosaurus ids or cell names,
        and does the gene column hold HGNC symbols, Ensembl ids, or integer
        token ids?
        """

        """
        ParquetFile
        https://github.com/ueshin/apache-arrow/blob/master/python/pyarrow/parquet.py
        """
        import pyarrow.parquet as pq

        f = self.fetch_de_shard(shard)
        pf = pq.ParquetFile(f)
        df = next(pf.iter_batches(batch_size=min(n_rows, 200_000))).to_pandas()

        out = {"shard_file": str(f),
               "file_mb": round(Path(f).stat().st_size / 1e6, 1),
               "total_rows_in_shard": pf.metadata.num_rows,
               "columns": list(df.columns),
               "dtypes": {c: str(t) for c, t in df.dtypes.items()},
               "head": df.head(3)}

        for c in df.columns:
            u = df[c].dropna().unique()
            out[f"distinct_{c}"] = list(u[:n_distinct])
            out[f"n_distinct_{c}"] = int(len(u))

        # Arrow-backed frames may give string[pyarrow], not object.
        obj_cols = [c for c in df.columns
                    if not pd.api.types.is_numeric_dtype(df[c])]

        # which of our cell-line vocabularies actually appears?
        cl_path = self.root_tahoe / "metadata" / "cell_line_metadata.parquet"
        if cl_path.exists():
            cl = pd.read_parquet(cl_path)
            out["cell_line_metadata_columns"] = list(cl.columns)
            for key in ("Cell_ID_Cellosaur", "cell_name"):
                if key not in cl.columns:
                    continue
                cand = set(cl[key].dropna().astype(str))
                for c in obj_cols:
                    n = df[c].astype(str).isin(cand).sum()
                    if n:
                        out[f"MATCH cell_line_metadata.{key} -> DE.{c}"] = int(n)
        else:
            out["cell_line_metadata"] = f"not downloaded yet at {cl_path}"

        # which gene vocabulary matches our query genes?
        if genes is not None:
            gset = {str(g) for g in genes}
            for c in obj_cols:
                n = df[c].astype(str).isin(gset).sum()
                if n:
                    out[f"MATCH query genes -> DE.{c}"] = int(n)
            gm = self.root_tahoe / "metadata" / "gene_metadata.parquet"
            if gm.exists():
                g_meta = pd.read_parquet(gm)
                out["gene_metadata_columns"] = list(g_meta.columns)
                for gk in g_meta.columns:
                    hit = g_meta[gk].astype(str).isin(gset).sum()
                    if hit:
                        out[f"query genes look like gene_metadata.{gk}"] = int(hit)

        num = df.select_dtypes(include=[np.number])
        if not num.empty:
            prof = pd.DataFrame({
                "min": num.min(), "max": num.max(), "mean": num.mean(),
                "frac_negative": (num < 0).mean(),
                "n_unique": num.nunique(),
            }).round(4)
            out["numeric_profile"] = prof
            # A statistic used for WTCS must be signed and roughly symmetric.
            out["signed_candidates"] = [
                c for c in num.columns
                if 0.2 < float((num[c] < 0).mean()) < 0.8]

        try:
            out["resolved_columns"] = self.resolve_de_columns(list(df.columns))
        except KeyError as e:
            out["resolved_columns"] = f"unresolved: {e}"

        out["matches"] = [k for k in out if k.startswith("MATCH")]
        if not out["matches"]:
            out["verdict"] = ("no filter vocabulary matched this shard -- "
                              "compare distinct_* values against what you pass")
        return out

    def list_de_shards(self) -> List[str]:
        """Repo-relative paths of the pseudobulk DE parquet shards."""
        from huggingface_hub import HfApi
        pre = "metadata/pseudobulk_differential_expression/"
        files = HfApi().list_repo_files(self._TAHOE_REPO, repo_type="dataset")
        return sorted(f for f in files if f.startswith(pre) and f.endswith(".parquet"))

    def shard_index(self, column: str = "Cell_ID_Cellosaur",
                    memory_limit: str = "4GB", verbose: bool = True,
                    stride: int = 1, force: bool = False,):
        """Accumulate a per-shard min/max index over `column` from parquet footers.

        Incremental: results merge into one canonical file keyed by shard index,
        so refining from stride=20 to stride=6 reads only the ~119 footers you
        do not already have, not all 172. Reads footer statistics only, never
        data.

        Each DE shard holds exactly one cell line (lo == hi), but the blocks are
        not in sorted id order, so this index is the only way to know which
        shards to scan.
        """
        import duckdb
        import time

        canon = self.root_tahoe / f"_shard_index_{column}.parquet"
        have = pd.DataFrame(columns=["i", "shard", "lo", "hi"])
        if canon.exists() and not force:
            have = pd.read_parquet(canon)

        all_shards = self.list_de_shards()
        want = list(range(0, len(all_shards), max(1, stride)))
        if want[-1] != len(all_shards) - 1:
            want.append(len(all_shards) - 1)
        known = set(have["i"].astype(int)) if len(have) else set()
        todo = [i for i in want if i not in known]

        if verbose:
            print(f"shard index: {len(known)} cached, {len(todo)} to read "
                  f"(stride={stride}, {len(all_shards)} shards total)", flush=True)
        if not todo:
            return have.sort_values("i").reset_index(drop=True)

        con = duckdb.connect()
        con.execute(f"SET memory_limit='{memory_limit}'")
        con.execute("INSTALL httpfs; LOAD httpfs;")
        tok = os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")
        if tok:
            con.execute(f"CREATE OR REPLACE SECRET hf_tok "
                        f"(TYPE huggingface, TOKEN '{tok}')")

        rows, _t0 = [], time.time()
        for j, i in enumerate(todo):
            f = all_shards[i]
            url = f"hf://datasets/{self._TAHOE_REPO}/{f}"
            try:
                m = self._duck_retry(con, f"""
                    SELECT min(coalesce(stats_min_value, stats_min)) AS lo,
                           max(coalesce(stats_max_value, stats_max)) AS hi
                    FROM parquet_metadata('{url}')
                    WHERE path_in_schema = '{column}'""")
                rows.append({"i": i, "shard": f, "lo": m["lo"][0], "hi": m["hi"][0]})
            except Exception as e:
                rows.append({"i": i, "shard": f, "lo": None, "hi": None,
                             "error": str(e)[:80]})
            if verbose and ((j + 1) % 5 == 0 or j == 0 or j == len(todo) - 1):
                el = time.time() - _t0
                eta = el / (j + 1) * (len(todo) - j - 1)
                print(f"  footer {j+1}/{len(todo)}  [shard {i}]  "
                      f"elapsed={el:5.0f}s  eta={eta:5.0f}s", flush=True)
            if (j + 1) % 25 == 0:                     # checkpoint mid-flight
                pd.concat([have, pd.DataFrame(rows)], ignore_index=True) \
                  .drop_duplicates("i").sort_values("i").to_parquet(canon)

        idx = (pd.concat([have, pd.DataFrame(rows)], ignore_index=True)
                 .drop_duplicates("i").sort_values("i").reset_index(drop=True))
        idx.to_parquet(canon)
        return idx

    def organ_cell_lines(self) -> List[str]:
        """Cellosaurus ids for one organ, from the cached cell_line_metadata."""
        f = self.root_tahoe / "metadata" / "cell_line_metadata.parquet"
        if not f.exists():
            self.download_tahoe(de=False)

        cl = pd.read_parquet(f)
        
        sel = cl[cl["Organ"] == self.organ]

        return sorted(sel["Cell_ID_Cellosaur"].dropna().astype(str).unique())

    def shard_sizes(self, shards: Optional[Sequence[str]] = None) -> pd.DataFrame:
        """Byte size of each DE shard, from the repo file listing."""
        from huggingface_hub import HfApi
        info = HfApi().repo_info(self._TAHOE_REPO, repo_type="dataset",
                                 files_metadata=True)
        sz = {s.rfilename: (s.size or 0) for s in info.siblings}
        keep = list(shards) if shards is not None else self.list_de_shards()
        return pd.DataFrame({"shard": keep,
                             "bytes": [sz.get(f, 0) for f in keep]})

    def download_shards(self, shards: Sequence[str],
                        max_workers: int = 8, dry_run: bool = False):
        """Download a specific set of DE shards in parallel.

        Why this exists: DuckDB's `hf://` reader fetches column chunks over
        single-threaded HTTP range requests, so a shard that actually contains
        matching rows costs minutes, not seconds -- while a shard that can be
        skipped entirely via footer statistics costs seconds. Timing a scan on
        skippable shards badly underestimates the real cost. Downloading the
        shards you need, in parallel, then querying them locally is far faster.

        Set HF_HUB_ENABLE_HF_TRANSFER=1 (pip install hf_transfer) for best speed.
        """
        from huggingface_hub import snapshot_download

        sizes = self.shard_sizes(shards)
        total = int(sizes["bytes"].sum())
        print(f"{len(shards)} shards, {total/1e9:.2f} GB "
              f"(median {sizes['bytes'].median()/1e6:.0f} MB/shard)")
        if dry_run:
            return sizes

        snapshot_download(repo_id=self._TAHOE_REPO, repo_type="dataset",
                          allow_patterns=list(shards), local_dir=self.root_tahoe,
                          max_workers=max_workers)
        return

    def index_coverage(self, column: str = "Cell_ID_Cellosaur"):
        """Distinguish "not sampled yet" from "not in the DE table at all".

        Refining the stride only helps for lines whose block is smaller than the
        gap. If the index already resolves ~every line the table contains, a
        target that is still missing is absent from the DE data (the atlas ships
        metadata for more lines than survive QC), and no stride will find it.
        """
        idx = self.shard_index(column=column, stride=10 ** 6, verbose=False)
        seen = sorted(idx["lo"].dropna().astype(str).unique())
        blocks = (idx.dropna(subset=["lo"]).groupby("lo")["i"]
                     .agg(["min", "max", "count"])
                     .sort_values("min"))
        out = {"n_probes": int(len(idx)),
               "n_lines_seen": len(seen),
               "lines_seen": seen,
               "block_probe_counts": blocks["count"].describe().round(1)}

        f = self.root_tahoe / "metadata" / "cell_line_metadata.parquet"
        if f.exists():
            cl = pd.read_parquet(f)
            allc = set(cl["Cell_ID_Cellosaur"].dropna().astype(str))
            out["n_lines_in_metadata"] = len(allc)
            out["in_metadata_not_in_index"] = sorted(allc - set(seen))
            gap = int(idx["i"].diff().median()) if len(idx) > 1 else 1
            out["smallest_block_probes"] = int(blocks["count"].min())
            out["verdict"] = (
                f"index resolves {len(seen)} lines from {len(idx)} probes at "
                f"gap~{gap}. Every block seen spans >= "
                f"{int(blocks['count'].min())} probe(s). Lines still absent are "
                f"most likely NOT in the DE table (metadata lists "
                f"{len(allc)}), rather than under-sampled.")
        return out

    def find_shards_for(self, values: Sequence[str],
                        column: str = "Cell_ID_Cellosaur",
                        report: bool = True, **kw) -> List[str]:
        """Shards that can contain any of `values`, via footer statistics.

        The DE table is grouped by cell line (lo == hi per shard) but the blocks
        are NOT in sorted id order, so a single min-to-max band would cover
        almost everything. We instead expand +/- one stride around each
        individual hit and union the results.

        With a coarse index some target lines may fall entirely between probes;
        those are reported as `missing` and need a finer stride.
        """
        kw.setdefault("stride", 10 ** 6)      # read nothing new by default
        kw.setdefault("verbose", False)
        idx = self.shard_index(column=column, **kw)
        ok = idx.dropna(subset=["lo", "hi"]).copy()
        all_shards = self.list_de_shards()
        if ok.empty:
            warnings.warn("no usable footer statistics; scanning all shards.")
            return list(all_shards)

        vals = sorted({str(v) for v in values})
        hit_mask = [any(str(r.lo) <= v <= str(r.hi) for v in vals)
                    for r in ok.itertuples()]
        hits = ok[hit_mask]

        diffs = ok["i"].diff().dropna()
        step = int(diffs.median()) if len(diffs) else 1

        keep_i = set()
        for i in hits["i"].astype(int):
            keep_i.update(range(max(0, i - step), min(len(all_shards), i + step + 1)))
        shards = [all_shards[i] for i in sorted(keep_i)]

        found = sorted({str(r.lo) for r in hits.itertuples()} |
                       {str(r.hi) for r in hits.itertuples()}) if len(hits) else []
        found = [f for f in found if f in vals]
        missing = [v for v in vals if v not in found]

        if report:
            print(f"targets: {len(vals)} | located: {len(found)} | "
                  f"missing: {len(missing)}")
            print(f"shards to scan: {len(shards)}/{len(all_shards)} "
                  f"({len(shards)/len(all_shards):.0%})  "
                  f"~{len(shards)*4.3/60:.0f} min at 4.3 s/shard")
            if missing:
                print(f"  NOT located: {missing}")
                n_seen = ok["lo"].dropna().astype(str).nunique()
                print(f"  index resolves {n_seen} distinct {column} values "
                      f"from {len(ok)} probes (gap~{step}).")
                print("  -> run index_coverage(): if these are absent from the "
                      "DE table, no stride will find them and you should "
                      "proceed without them.")

        if missing:
            warnings.warn(
                f"{len(missing)} target(s) not located at stride~{step}; their "
                "conditions will be silently absent from the result. Rebuild "
                "the index with a finer stride before the full scan.")
        if len(shards) / len(all_shards) > 0.9:
            warnings.warn("pruning saves nothing; consider download_tahoe(de=True).")
        return shards

    def probe_tahoe(self, genes=None, organs=None, cell_lines=None,
                    n_shards: int = 1, **kw):
        """Time a filtered scan of `n_shards` and extrapolate the full run.

        Run this BEFORE launching a full scan. If one shard takes 20 s, 1026
        shards take ~6 h at threads=1 -- which is what a silent multi-hour
        `load_tahoe_de` actually means. The answer that number should drive is
        stream-vs-download, not patience.
        """
        import time
        shards = self.list_de_shards()
        # Sample evenly across the range: the first n shards are unrepresentative
        # if the table is clustered by cell line, and would time an empty scan.
        idx = np.unique(np.linspace(0, len(shards) - 1, n_shards).astype(int))
        picked = [shards[i] for i in idx]

        t0 = time.time()
        S, cond = self.load_tahoe_de(
            genes=genes, organs=organs, cell_lines=cell_lines,
            save_to=False, _shard_subset=picked, **kw)
        dt = time.time() - t0
        est = dt / max(len(picked), 1) * len(shards)
        return {"n_shards_total": len(shards), "n_probed": len(picked),
                "shard_indices": idx.tolist(),
                "seconds_per_shard": round(dt / max(len(picked), 1), 1),
                "estimated_full_scan_hours": round(est / 3600, 2),
                "conditions_found": int(S.shape[1]) if S.size else 0,
                "note": ("no rows in the sampled shards -- this timing "
                         "measures fetching nothing"
                         if not S.size else "timing reflects real work"),
                "recommendation": (
                    "INCONCLUSIVE: sampled shards held no matching rows. Run "
                    "shard_index()/find_shards_for() to locate the shards that "
                    "do, then probe those." if not S.size
                    else "stream is fine" if est < 1800
                    else "too slow to stream -- use download_tahoe(de=True) "
                         "once, then mode='download'")}

    def _scan_de_sharded(self, con, shards, select_sql, where_sql,
                         checkpoint_dir, retries, wait_s, verbose=True):
        """Query shards one at a time, checkpointing each to parquet.

        Turns an opaque multi-hour scan into a resumable job: killed halfway,
        the next call skips completed shards. Also caps peak memory at one
        shard's filtered output rather than the whole result.
        """
        import time
        ck = Path(checkpoint_dir); create_dir(ck)
        t0, parts, n = time.time(), [], len(shards)

        for i, f in enumerate(shards):
            # Name by shard, never by loop position: probes and full scans pass
            # different subsets, and a positional name silently serves shard 1's
            # cached result when you asked for shard 1025.
            out = ck / f"part_{Path(f).stem}.parquet"
            if out.exists():
                parts.append(out)
                continue
            url = f"hf://datasets/{self._TAHOE_REPO}/{f}"
            q = f"{select_sql} FROM read_parquet('{url}') {where_sql}"
            df = self._duck_retry(con, q, retries=retries, wait_s=wait_s)
            df.to_parquet(out)
            parts.append(out)
            if verbose:
                done = i + 1
                el = time.time() - t0
                eta = el / done * (n - done)
                print(f"  shard {done}/{n}  rows={len(df):>7,}  "
                      f"elapsed={el/60:5.1f}m  eta={eta/60:5.1f}m", flush=True)

        if not parts:
            return pd.DataFrame()
        return pd.concat((pd.read_parquet(f) for f in parts), ignore_index=True)

    @staticmethod
    def _tahoe_cache_key(organs, cell_lines, drugs, genes) -> str:
        """Short hash of the query so different filters get different files.

        Keying the cache on the directory alone means changing `genes` (which
        run() derives from X.columns, so it moves whenever n_hvg/min_share
        change) silently returns a stale signature matrix.
        """
        payload = json.dumps({
            "organs": sorted(map(str, organs)) if organs else None,
            "cell_lines": sorted(map(str, cell_lines)) if cell_lines else None,
            "drugs": sorted(map(str, drugs)) if drugs else None,
            "genes": sorted({str(g) for g in genes}) if genes is not None else None,
        }, sort_keys=True)
        return hashlib.sha1(payload.encode()).hexdigest()[:12]

    @staticmethod
    def _duck_retry(con, query: str, retries: int = 8, wait_s: float = 5.0):
        """Run a DuckDB query, backing off on HTTP 429/5xx.

        HF returns 429 when a scan fans out across many shards at once. Backoff
        is exponential and capped; non-HTTP errors are re-raised immediately so
        real bugs are not retried.
        """
        import time
        last = None
        for attempt in range(retries):
            try:
                return con.execute(query).df()
            except Exception as e:                       # duckdb.HTTPException
                msg = str(e)
                retryable = any(c in msg for c in
                                ("429", "Too Many Requests", "503", "504",
                                 "500", "connection", "timeout", "Timeout"))
                if not retryable:
                    raise
                last = e
                sleep = min(wait_s * (2 ** attempt), 120.0)
                warnings.warn(
                    f"Tahoe query attempt {attempt + 1}/{retries} failed "
                    f"({msg.splitlines()[0][:90]}); retrying in {sleep:.0f}s.")
                time.sleep(sleep)
        raise RuntimeError(
            f"Tahoe query failed after {retries} attempts. Last error: {last}\n"
            "Options: set HF_TOKEN, lower threads=, or switch to "
            "mode='download' for a one-off local copy."
        ) from last


    # ===========================================================================
    # 3c. Multi-compartment model: programs per compartment + their coupling
    # ===========================================================================

    #: PDAC program markers, per compartment. Override freely.
    """
    link: https://pathway-viewer.toolforge.org/embed/WP5078
    better than a generic Treg pathway, because WP5078 describes how pancreatic cancer cells create a microenvironment 
    by secreting cytokines and chemokines that recruit or activate stromal cells, 
    promoting desmoplasia and immune evasion, 
    with stromal cells inhibiting CD8+ T cells through factors like IL-10, TGFβ, PD-L1, and IDO. 
    That's organised by emitting cell type, which maps onto your compartments.


        1)  WikiPathways WP5078 (T cell modulation and desmoplasia in PDAC)
            Checkpoint ligands + immunosuppressive enzymes EMITTED by the tumour.
            Scoreable even though your T-cell compartment is not: theta_mal is high,
            so this measures the emitting side where the receptor side is pure prior.

            Galectins collide with immune_evasion (LGALS1/3/9, all three), and 
            IL10 collides cross-compartment with mdsc_suppressive. Shared genes make two scores partly the same measurement.

        
        2) "desmoplastic_secretome"
            What the tumour secretes to build its own desmoplasia.
            Why desmoplastic_secretome is the most valuable addition?
                It gives you a directed hypothesis, which none of your current programs do. 
                WP5078 says:
                - the tumour emits PDGF/SHH/TGFβ and 
                - stellate cells respond by depositing ECM. 
            So malignant.desmoplastic_secretome ↔ fibroblast.myCAF is a specific, mechanistically motivated pair — 
            not one of ~50 undirected correlations. Test it as a single pre-specified hypothesis rather than letting BH dilute it.

        3) A program score is a mean of z-scores — meaningful only if the genes co-vary. 
            Fibrogenic (PDGF/SHH/TGFβ), angiogenic (VEGF/PGF), and immunosuppressive (galectins/IL10) signals move independently, 
            so the mean dilutes all three. Splitting them costs nothing and gains interpretability.

        4) CCN2 (formerly CTGF) is a well-established PDAC desmoplasia driver and belongs here more than VEGF does.
            IL1A/IL1B are the tumour-derived signals that induce the iCAF state specifically — 
            which matters because it gives you a directed prediction with a sign: 
                tumour IL-1 → iCAF up. 
            Your one replicated finding is prolif ↔ iCAF negative, so this is a genuine test rather than a fishing expedition.        

        5) CAF
            The terms myCAF, iCAF, and apCAF stand for three main types of cancer-associated fibroblasts (CAF) found inside tumors. 
            They are myofibroblastic CAF (myCAF), inflammatory CAF (iCAF), and antigen-presenting CAF (apCAF)

            One caveat worth holding: apCAF is defined by MHC-II without costimulation, 
            and MHC-II is far more abundantly expressed by macrophages and B cells than by fibroblasts. 
            Given your compartments are correlated at r≈0.9, an apCAF score is at real risk of reading myeloid bleed rather than fibroblast biology. 
            SLPI and CD74 are the least myeloid-specific of the set, which is thin protection.

            So I'd treat any apCAF result with more scepticism than myCAF/iCAF, and check it against macrophage.TAM 
            — if fibroblast.apCAF correlates strongly with a myeloid program, bleed is the likelier explanation than antigen-presenting fibroblasts.


        6) Ductal
            Normal-tissue compartments. Their programs largely index residual normal parenchyma 
            -- i.e. the same cellularity gradient that derailed gene-level clustering -- 
            so treat any coupling involving them as a composition artifact until shown otherwise.

        7) Comments:
            For ANGPT2 and IL1B, keep the compartment where the gene is canonically expressed
            — ANGPT2 is an endothelial Weibel-Palade product, IL1B is predominantly myeloid. 
            Putting them in the malignant secretome would correlate a tumour score against the very compartments those genes define. 
            
            VEGFA is the reverse: it's a classic tumour-hypoxia product, so keep it in malignant and drop it from the already-thin macrophage.SPP1.       

        8) Comments2:
            angiogenic_secretome won't survive. VEGFD and ANGPT2 are both absent from your matrix (you checked), 
            and removing ANGPT2 leaves VEGFA, VEGFB, VEGFC, PGF — four, under the min_genes=5 you'd want. 
            Either add ANGPTL4, HIF1A, SLC2A1, CA9 (hypoxia-response, tumour-intrinsic) or drop the program.

            SHH and CCN2 are likely dead too. You confirmed CCN2 isn't in genes_full, though the harmonised ref_h may now recover it 
            — worth re-checking after the re-run. SHH was present but is low-abundance and may not clear min_counts=10.

            Also, malignant.emt has FN1 and SPARC, and fibroblast.myCAF has COL1A1/POSTN/THBS2 — no shared genes, 
            but they're the same ECM biology in two compartments correlated at r≈0.9. 
            An emt↔myCAF coupling would be hard to distinguish from bleed. Not a bug, but worth remembering when interpreting it.

            Finally, 28 programs means couple_compartments runs ~340 cross-compartment pairs. 
            Your prolif↔iCAF sat at FDR 0.041–0.050 with far fewer tests — it will not survive that burden. 
            Run the WP5078 additions as a separate pre-specified family rather than folding them into one BH correction.

    """

    PROGRAMS: Dict[str, Dict[str, Sequence[str]]] = {
        "malignant": {
            "basal":     ["KRT5", "KRT6A", "KRT14", "KRT17", "KRT81", "S100A2", "TP63", "SPRR3", "DHRS9", "VGLL1", "SERPINB3", "LY6D"],
            "classical": ["GATA6", "TFF1", "TFF2", "TFF3", "AGR2", "LGALS4", "CEACAM6", "CLDN18", "REG4", "ANXA10", "CTSE", "MUC13"],
            "emt":       ["VIM", "ZEB1", "SNAI2", "CDH2", "FN1", "SPARC"],
            "prolif":    ["MKI67", "TOP2A", "CCNB1", "AURKA", "BIRC5", "PLK1"],
            # galectins stay in immune_evasion only; FAP and IL10 removed entirely
            "immune_evasion": ["CD274", "PDCD1LG2", "CD276", "VTCN1", "VSIR", "IDO1", "LGALS1", "LGALS3", "LGALS9", "PVR", "NT5E", "ENTPD1"],
            # Antigen-presentation loss: LOW score = escape active. Read inverted.
            "antigen_presentation": ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "NLRC5"],
            "desmoplastic_secretome": ["PDGFA", "PDGFB", "PDGFC", "PDGFD", "SHH", "TGFB1", "TGFB2", "TGFB3", "CCN2", "IL1A"], # "IL1B"
            "angiogenic_secretome":   ["VEGFA", "VEGFB", "VEGFC", "VEGFD", "PGF"],  # , "ANGPT2"

        },
        "fibroblast": {
            "myCAF":  ["ACTA2", "TAGLN", "POSTN", "COL1A1", "COL1A2", "THBS2", "CTHRC1", "INHBA"],
            "iCAF":   ["IL6", "CXCL12", "PDGFRA", "HAS1", "HAS2", "CFD", "LMNA", "C3", "CXCL14"],
            # Elyada's human apCAF signature rests on MHC class II, and your list already has three of those:
            # "SAA3" is a maouse gene
            "apCAF":  ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DQA1", "SLPI"],
        },


        "immune": {
            "cytotoxic":  ["CD8A", "GZMB", "PRF1", "GNLY", "NKG7", "IFNG"],
            "exhaustion": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX"],
            "treg":       ["FOXP3", "IL2RA", "IKZF2", "TNFRSF4", "CCR8"], # removed: "CTLA4", it's already in exhaustion:
            "costim":     ["CD28", "ICOS", "CD226", "CD40LG", "TNFRSF9"],   # opposes `exhaustion`
        },

        "macrophage": {
            "M1":   ["IL1B", "TNF", "CXCL9", "CXCL10", "NOS2",  "CCL5"], # removed "IL6",
            "TAM":  ["CD163", "MRC1", "MSR1", "TREM2", "APOE", "C1QA", "C1QB"],
            "SPP1": ["SPP1", "MARCO",  "MMP9"],  # removed "FN1", "VEGFA",
            "mdsc_suppressive": ["ARG1", "IL10", "CCL17", "CCL22", "IL4", "IL13", "CD80", "CD86", "CD40"],
        },
        
        "endothelial": {
            "tip_angio": ["ESM1", "ANGPT2", "DLL4", "APLN", "CXCR4", "KDR"],
            "lymphatic": ["PROX1", "LYVE1", "PDPN", "CCL21", "FLT4"],
            "activated": ["VCAM1", "ICAM1", "SELE", "ACKR1"],
        },
        "humoral": {
            # Tertiary Lymphoid Structures (TLS)
            "TLS":    ["CXCL13", "CR2", "FDCSP", "MS4A1", "CD79A", "CCL19"],
            "plasma": ["MZB1", "JCHAIN", "XBP1", "DERL3", "TNFRSF17"],
        },

        "ductal": {
            "normal_duct": ["CFTR", "SCTR", "AQP1",  "MMP7"],  # removed "SPP1",
            "ADM":         ["SOX9", "KRT19", "ONECUT1", "HNF1B"],
        },
        "acinar": {
            "acinar_identity": ["PRSS1", "CPA1", "CPB1", "CELA3A", "CTRB1", "PNLIP", "CLPS"],
            "stress":          ["REG1A", "REG3A", "CLU", "MT1G"],
        },
    }

    SIGNATURES: Dict[str, Dict[str, Sequence[str]]] = {
        # --------------------------------------------------------------------------- #
        # ACINAR
        # --------------------------------------------------------------------------- #

        'ACINAR_CORE': {
            # serine proteases / zymogens -- the highest-abundance transcripts in pancreas
            "proteases": [
                "PRSS1", "PRSS2", "PRSS3", "CTRB1", "CTRB2", "CTRC", "CTRL",
                "CELA2A", "CELA2B", "CELA3A", "CELA3B", "CELA1",
                "CPA1", "CPA2", "CPB1", "CPO",
            ],
            # lipases, amylases, other digestive enzymes
            "lipase_amylase": [
                "PNLIP", "PNLIPRP1", "PNLIPRP2", "CLPS", "CEL", "PLA2G1B",
                "AMY2A", "AMY2B", "AMY1A",
            ],
            # zymogen granule / regulated secretory machinery
            "secretory_machinery": [
                "SYCN", "GP2", "CUZD1", "ZG16", "PDIA2", "ERP27", "SEC11C",
                "AQP8", "AQP12A", "AQP12B", "SERPINI2", "KLK1", "DPEP1",
            ],
            # acinar identity transcription factors
            "identity_tf": [
                "PTF1A", "RBPJL", "BHLHA15", "NR5A2", "GATA4", "MYRF",
            ],
        },

        'ACINAR_FLAGGED': {
            "regenerating_adm": ["REG1A", "REG1B", "REG3A", "REG3G", "REG4"],
            "shared_epithelial": ["SPINK1", "RNASE1", "MUC1", "TFF1", "TFF2", "LGALS4"],
            "stress_metaplasia": ["ONECUT1", "SOX9", "KRT19", "CFTR"],
        },

        # --------------------------------------------------------------------------- #
        # FIBROBLAST / STELLATE
        # --------------------------------------------------------------------------- #
        'FIBROBLAST_CORE': {
            "collagens": [
                "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL5A3",
                "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL10A1", "COL11A1",
                "COL12A1", "COL15A1", "COL16A1",
            ],
            "ecm_glycoprotein": [
                "FBN1", "FBLN1", "FBLN2", "FBLN5", "THBS2", "TNC", "POSTN",
                "VCAN", "LUM", "DCN", "FMOD", "OGN", "ASPN", "PRELP", "MGP",
                "ELN", "EMILIN1", "LTBP1", "LTBP2", "MFAP2", "MFAP5",
                "NID2", "SPON1", "ISLR",
            ],
            "remodeling_crosslink": [
                "LOX", "LOXL1", "LOXL2", "PLOD1", "PLOD2", "SERPINH1",
                "MMP2", "MMP11", "MMP14", "MMP23B", "TIMP2", "TIMP3",
                "ADAM12", "ADAMTS2", "ADAMTS12", "SULF1", "PCOLCE",
            ],
            "identity_receptor": [
                "PDGFRA", "PDGFRB", "FAP", "THY1", "CD248", "ITGA11", "DDR2",
                "LRRC15", "ANTXR1", "GLI1", "PDPN",
            ],
            "identity_tf": [
                "TWIST1", "TWIST2", "PRRX1", "SNAI2", "FOXF1", "FOXF2", "EBF1",
            ],
            "matricellular": [
                "CCN1", "CCN2", "CCN3", "CCN5", "SPARCL1",
            ],
        },

        # Schwann — as an attribution control, not a biology panel. PLP1, SOX10, MPZ, S100B, EGR2, CDH19, PMP22 are Schwann-restricted with no ECM annotation. Peng's reference has no neural cell type, so nerve-derived transcripts must land somewhere, and this panel measures where. That makes it the same class of instrument as your "acinar core in fibroblast" check: known expected answer, tests the pipeline rather than describing it.
        # Crucially it must exclude NGFR and RELN, which you've now correctly placed in quiescent_stellate. Keeping them in both would make the two panels correlate by construction — the exact thing the zero-duplicate design avoids.
        'CONTROL': {
            # No neural cell type in the Peng reference, so nerve-derived transcripts
            # have nowhere to go. Non-zero here means the compartment is absorbing them.
            # NGFR and RELN deliberately excluded -- shared with quiescent stellate.
            "schwann_attribution": ["PLP1", "SOX10", "MPZ", "S100B", "EGR2", "CDH19", "PMP22"],
        },

        'FIBROBLAST_FLAGGED': {
            "pericyte_smc_overlap": ["ACTA2", "TAGLN", "MYL9", "CNN1", "TPM1",
                                    "TPM2", "PDLIM3", "RGS5", "NOTCH3", "MYH11"],
            "broadly_expressed":   ["VIM", "SPARC", "FN1", "THBS1", "TIMP1", "S100A4"],
            "immune_shared":       ["IL6", "CXCL1", "CXCL2", "CXCL12", "CCL2",
                                    "C3", "CFD", "C1S", "C1R", "LIF", "IL11",
                                    "HAS1", "HAS2", "PTGDS"],
            "quiescent_stellate":  ["LRAT", "DES", "NGFR", "CYGB", "RELN", "PDZRN4"],

            # FLAGGED rather than CORE because your own tier definition is "real markers, but each carries an attribution or interpretation risk." 
            # This panel's risk is circularity: 
            # it was derived from this cohort's fibroblast residuals, so scoring it here isn't independent evidence. 
            # That's exactly what FLAGGED exists to isolate — score CORE alone first, then add FLAGGED and see whether conclusions move.
            "quiescent_fibroblast_derived": [
                                    "C7", "MFAP4", "COL14A1", "THBS4", "ABCA6", "ABCA8", "ABCA9",
                                    "ABCA10", "CHRDL1", "SFRP1", "TNXB", "ADAMTSL3", "GPC3"],
        },

        'RISK': {
            # ---- acinar ----
            "regenerating_adm": (
                "REG family is induced by inflammation and acinar-to-ductal metaplasia. "
                "In adjacent-normal PDAC tissue (pancreatitis, ADM) these are up for reasons "
                "unrelated to acinar abundance, so they corrupt the normal arm of the contrast."
            ),
            "shared_epithelial": (
                "Expressed by ductal and malignant cells. SPINK1 in particular is elevated in "
                "PDAC -- including it puts a tumour-up gene inside a compartment that is "
                "tumour-down, which inverts the contrast."
            ),
            "stress_metaplasia": (
                "Ductal-programme genes switched on during ADM. High values mean acinar cells "
                "are becoming ductal, i.e. the compartment identity itself is unstable -- "
                "which breaks the same-compartment-two-states assumption."
            ),
            # ---- fibroblast ----
            "pericyte_smc_overlap": (
                "Shared with pericytes and vascular smooth muscle. Direct attribution risk "
                "against the endothelial compartment; RGS5 is a pericyte marker first. "
                "Any fibroblast<->endothelial correlation built on these is suspect."
            ),
            "broadly_expressed": (
                "Expressed at meaningful levels by most cell types including malignant. "
                "Low compartment share, so BayesPrism's split is prior-dominated -- exactly "
                "the genes that manufacture cross-compartment coupling."
            ),
            "immune_shared": (
                "The iCAF programme proper, but overlapping macrophage and endothelial "
                "signatures. These are the genes most likely to be carrying the "
                "pan-compartment inflammatory axis rather than fibroblast biology."
            ),
            "quiescent_stellate": (
                "Marks quiescent PSCs, near-absent in PDAC. Useful for the NORMAL arm and for "
                "a chronic-pancreatitis comparator; near-zero and prior-dominated in tumour."
            ),
            "quiescent_fibroblast_derived": (
                "Derived from the PAAD fibroblast-vs-acinar residual in this cohort "
                "(median resid -0.566, 13/13 negative, AveExpr_fib 5.2-9.9, no abundance "
                "gradient). Scoring it on the same data is circular -- it describes rather "
                "than tests. Independent evidence requires a second cohort or reference. "
                "Distinct from quiescent_stellate, which is canonical and largely "
                "unmeasurable here (LRAT absent)."
            ),
            "schwann_attribution": (
                "Not a biology panel -- an attribution control. The Peng reference has no "
                "neural cell type, so nerve-derived transcripts must be absorbed by some "
                "compartment. Expected value is ~zero; non-zero means the compartment is "
                "taking up nerve signal. NGFR and RELN are deliberately excluded (shared "
                "with quiescent_stellate) so the two panels stay independent."
            ),            
        }
    }

    @classmethod
    def check_signatures(cls):
        from collections import Counter
        sig = {k: v for k, v in cls.SIGNATURES.items() if k != "RISK"}
        c = Counter(g for ps in sig.values() for gs in ps.values() for g in gs)
        dup = {g: n for g, n in c.items() if n > 1}
        assert not dup, f"genes in >1 panel (inflates panel correlation by k/n): {dup}"
        missing = [f"{t}.{p}" for t, ps in sig.items() if "FLAGGED" in t or "CONTROL" in t
                for p in ps if p not in cls.SIGNATURES.get("RISK", {})]
        assert not missing, f"FLAGGED/CONTROL panels without RISK text: {missing}"
        return True

    def compartment_matrix(
            self,
            cell_name: str,
            min_share: float = 0.0,
            min_counts: float = 1.0,
            samples: Optional[Sequence[str]] = None,
            min_theta: float = 0.02,
            decouple: bool = True,
            batch: Optional[pd.Series] = None,
            batch_scale: bool = False,
            cond_warn: float = 100.0,
        ) -> pd.DataFrame:
        """
        compartment_matrix receives cell_name and build
        cns = self.build_CellNameState_from_full_Z()
        therefore, build_CellNameState_from_full_Z is for any cell type (cell name)

        log2-CPM profile of one compartment, normalised WITHIN that compartment.

        Compartment-internal CPM is what makes cross-compartment comparison
        meaningful: raw Z scales with that compartment's abundance, so
        correlating two raw Z matrices mostly recovers the composition, not the
        biology.

        `batch` residualises a cohort label alongside theta and depth. Without
        it, a multi-cohort bulk passes its batch structure straight through:
        BayesPrism corrects reference-vs-bulk mismatch, not batch WITHIN the
        bulk, so every compartment inherits the same cohort axis and
        cross-compartment correlations measure cohort membership. Observed on
        TCGA+CPTAC3 PAAD: compartment PC1 explained ~38% of variance, separated
        the cohorts at point-biserial r = 0.89, and its gene loadings agreed
        across compartments at r = 0.996.

        One coefficient per gene per batch level, so the correction absorbs
        gene-specific batch effects -- including the length-dependent ones that
        come from differing library prep, which a single global shift cannot
        touch. It is still a LOCATION correction: see `batch_scale` and the
        note in `harmonize_scores` on heteroscedastic cohort effects.

        Parameters
        ----------
        batch : Series indexed by sample, or None
            Cohort label. Any dtype; converted to dummies with the first level
            absorbed into the intercept.
        batch_scale : bool
            Also equalise per-gene within-batch SD. Off by default: a cohort
            with genuinely more heterogeneous tumours would have real spread
            flattened, and with no overlapping samples the two cases cannot be
            distinguished. Report both versions when it matters.
        cond_warn : float
            Warn when cond(D) exceeds this. theta and log(lib) are near
            collinear by construction -- Z scales with theta, so
            log(lib) ~ log(theta) + log(depth) -- and the rank check below only
            catches EXACT collinearity. Simulated cond for that pair alone is
            ~430, so a rank-3-of-3 design can still be badly conditioned and
            split shared variance between the two columns arbitrarily.
        """

        cns = self.build_CellNameState_from_full_Z(
            self.Z_full, self.genes_full, cell_types=self.cell_types,
            cell_name=cell_name, samples=self.df_theta.index, df_theta=self.df_theta)
        
        Z = cns.Z
        if samples is not None:
            Z = Z.reindex([s for s in samples if s in Z.index])

        '''
        Drop samples with essentially none of this cell type. Their Z is the
        reference prior, not a measurement -- BayesPrism can return theta on
        the order of 1e-50, and such a sample will still receive a program
        score, usually an extreme one, and can single-handedly create a "subtype".
        '''
        th = self.df_theta[cell_name].reindex(Z.index)
        low = th < min_theta
        if low.any():
            warnings.warn(
                f"{cell_name}: dropping {int(low.sum())} sample(s) with theta < "
                f"{min_theta} (min observed {th.min():.2e}); their compartment "
                "profile is prior, not signal.")
            Z = Z.loc[~low]
            th = th.loc[Z.index]

        Z = Z.loc[:, Z.median(axis=0) >= min_counts]

        if min_share > 0 and cns.share is not None:
            sh = cns.share.reindex(index=Z.index, columns=Z.columns)
            keep = (sh >= min_share).mean(axis=0) >= 0.5
            '''
            The share filter selects on composition: a gene must clear
            min_share in half the samples, and in low-theta samples nothing
            can. So the surviving gene set is chosen by the high-theta half
            of the cohort. If theta differs systematically between batches,
            the gene set is batch-selected too -- report it rather than
            silently inheriting it.
            '''
            if batch is not None:
                b = batch.reindex(Z.index)
                per_b = {lv: set(Z.columns[(sh.loc[b == lv] >= min_share)
                                        .mean(axis=0) >= 0.5])
                        for lv in b.dropna().unique()}
                if len(per_b) > 1:
                    sets = list(per_b.values())
                    shared = set.intersection(*sets)
                    union = set.union(*sets)
                    if union and len(shared) / len(union) < 0.8:
                        warnings.warn(
                            f"{cell_name}: min_share={min_share} selects "
                            f"different gene sets per batch "
                            f"({ {k: len(v) for k, v in per_b.items()} }, "
                            f"Jaccard {len(shared)/len(union):.2f}). The gene "
                            "filter is partly a batch filter; consider "
                            "intersecting per-batch sets or lowering "
                            "min_share.")
            Z = Z.loc[:, keep]

        # A sample can lose all its genes to the filters above; its row sum is
        # then 0 and CPM becomes 0/0. Drop those before they poison the design
        # matrix with NaN/-inf and crash the least-squares fit.
        lib = Z.sum(axis=1)
        dead = ~np.isfinite(lib) | (lib <= 0)
        if dead.any():
            warnings.warn(
                f"{cell_name}: dropping {int(dead.sum())} sample(s) with zero "
                f"counts after filtering: {list(Z.index[dead])[:5]}")
            Z, lib, th = Z.loc[~dead], lib.loc[~dead], th.loc[Z.index[~dead]]
        if Z.shape[0] < 3 or Z.shape[1] < 2:
            raise ValueError(
                f"{cell_name}: only {Z.shape[0]} samples x {Z.shape[1]} genes "
                f"survive min_share={min_share}, min_counts={min_counts}, "
                f"min_theta={min_theta}. Relax the filters.")

        cpm = Z.div(lib, axis=0) * 1e6
        logx = np.log2(cpm + 1.0)
        logx = logx.loc[:, np.isfinite(logx.values).all(axis=0)]

        # Keep the covariates that built this matrix reachable. Recomputing
        # them outside means re-deriving every filter above, which silently
        # diverges the moment a filter changes.
        self._last_covariates = {"cell_name": cell_name, "theta": th.copy(),
                                "lib": lib.copy(),
                                "batch": (batch.reindex(logx.index).copy()
                                        if batch is not None else None)}

        if decouple and logx.shape[1] >= 2 and float(th.std()) > 0:
            # Same argument as prepare_malignant_matrix: shrinkage toward the
            # reference is graded by this compartment's own abundance, so
            # residualise on it before scoring programs.
            cols = [np.ones(len(th)), th.values, np.log(lib.values)]
            names = ["intercept", "theta", "log_lib"]

            if batch is not None:

                print("Fixing cohort batch effects...")

                b = batch.reindex(logx.index)
                if b.isna().any():
                    raise ValueError(
                        f"{cell_name}: batch is missing for "
                        f"{int(b.isna().sum())} sample(s), e.g. "
                        f"{list(b.index[b.isna()])[:5]}. Dropping them "
                        "silently would change the sample set between "
                        "compartments; supply a complete label.")
                lv = list(pd.unique(b))
                if len(lv) < 2:
                    warnings.warn(f"{cell_name}: batch has one level "
                                f"({lv[0]!r}); no batch term added.")
                else:
                    counts = b.value_counts()
                    if counts.min() < 3:
                        warnings.warn(
                            f"{cell_name}: batch level(s) with <3 samples "
                            f"{counts[counts < 3].to_dict()}; their "
                            "coefficient is fit on almost no data.")
                    # drop-first: the reference level lives in the intercept
                    for level in lv[1:]:
                        cols.append((b == level).values.astype(float))
                        names.append(f"batch[{level}]")

            D = np.column_stack(cols)

            if not np.isfinite(D).all():
                warnings.warn(f"{cell_name}: non-finite covariates; "
                            "skipping purity decoupling.")
            elif np.linalg.matrix_rank(D) < D.shape[1]:
                warnings.warn(
                    f"{cell_name}: covariates are collinear (rank "
                    f"{np.linalg.matrix_rank(D)} < {D.shape[1]}: {names}); "
                    "skipping purity decoupling.")
            else:
                cond = float(np.linalg.cond(D))
                if cond > cond_warn:
                    warnings.warn(
                        f"{cell_name}: cond(D) = {cond:.0f} over {names}. The "
                        "fitted values and residuals are fine, but individual "
                        "coefficients are not identified -- do not read beta "
                        "as 'the theta effect' or 'the batch effect' "
                        "separately.")
                try:
                    # lstsq can fail to converge on badly scaled input; the
                    # normal equations via pinv are slower but robust.
                    beta, *_ = np.linalg.lstsq(D, logx.values, rcond=None)
                except np.linalg.LinAlgError:
                    warnings.warn(f"{cell_name}: lstsq did not converge; "
                                "falling back to pinv.")
                    beta = np.linalg.pinv(D) @ logx.values
                logx = pd.DataFrame(
                    logx.values - D @ beta + logx.values.mean(axis=0),
                    index=logx.index, columns=logx.columns)
                self._last_covariates.update(beta=beta, design_names=names,
                                            cond=cond)

        # Location correction leaves per-gene within-batch dispersion intact.
        # See harmonize_scores: cohort effects on deconvolved values are often
        # heteroscedastic, and mean-centring alone leaves every extreme -- and
        # so every tail cluster -- coming from the wider cohort.
        if batch_scale and batch is not None:
            b = batch.reindex(logx.index)
            if b.nunique() > 1:
                grand = logx.std(axis=0)
                for level in b.unique():
                    m = (b == level).values
                    if m.sum() < 3:
                        continue
                    sub = logx.loc[m]
                    sd = sub.std(axis=0).replace(0.0, np.nan)
                    logx.loc[m] = (((sub - sub.mean(axis=0))
                                    .div(sd, axis=1).fillna(0.0)
                                    .mul(grand, axis=1))
                                + sub.mean(axis=0))

        return logx


    def compartment_readiness(
        self,
        compartment_map: Dict[str, str],
        min_share: float = 0.30,
        min_theta_median: float = 0.03,
        min_markers: int = 4,
    ) -> pd.DataFrame:
        """Go/no-go per compartment BEFORE scoring programs.

        A compartment with low theta has little compartment-specific evidence,
        so BayesPrism shrinks its Z toward the reference: the resulting
        "programs" are the reference profile re-scored, and any coupling you
        find is between two shrunken copies of Peng 2019, not between tumours.
        In PDAC bulk, B cells and endothelium routinely sit at 1-3% theta.

        Columns: median/min theta, genes surviving the share filter, and marker
        coverage per program. `verdict` flags compartments to drop.
        """
        programs = self.PROGRAMS
        rows = []
        for key, ct in compartment_map.items():
            if ct not in self.df_theta.columns:
                rows.append({"compartment": key, "cell_type": ct,
                            "verdict": "MISSING from df_theta"})
                continue
            th = self.df_theta[ct]
            try:
                M = self.compartment_matrix(ct, min_share=min_share, min_counts=MIN_COUNTS, batch=batch, )
                n_genes = M.shape[1]
            except Exception as e:
                rows.append({"compartment": key, "cell_type": ct,
                            "theta_median": round(float(th.median()), 4),
                            "verdict": f"matrix failed: {str(e)[:40]}"})
                continue

            progs = programs.get(key, {})
            cov = {pr: len([g for g in gs if g in M.columns])
                for pr, gs in progs.items()}
            thin = [pr for pr, n in cov.items() if n < min_markers]

            if not progs:
                verdict = "NO PROGRAMS defined -- add to mc.PROGRAMS"
            elif th.median() < min_theta_median:
                verdict = f"DROP - theta median {th.median():.3f} too low"
            elif len(thin) == len(cov):
                verdict = "DROP - no program has enough markers"
            elif thin:
                verdict = f"CAUTION - thin programs: {thin}"
            else:
                verdict = "OK"

            rows.append({"compartment": key, "cell_type": ct,
                        "theta_median": round(float(th.median()), 4),
                        "theta_min": round(float(th.min()), 4),
                        "n_genes_after_share": n_genes,
                        "marker_coverage": cov,
                        "verdict": verdict})
        return pd.DataFrame(rows)

    def program_scores(
        self,
        compartment_map: Optional[Dict[str, str]] = None,
        samples: Optional[Sequence[str]] = None,
        min_genes: int = 3,
        min_share: float = 0.1,
        min_counts: int = 10,
        batch: pd.Series | Optional[str] = None,
        **mat_kw,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Continuous program scores per sample, across compartments.

        `compartment_map` maps a PROGRAMS key to the cell-type name in your
        reference, e.g.
            {"malignant": "Ductal cell type 2",
             "fibroblast": "Fibroblast cell",
             "immune": "T cell"}

        Each program is the mean z-score (across samples) of its marker genes
        within its own compartment. Scores are contrasts within a compartment,
        so a sample whose whole profile is depressed does not automatically
        score low on everything -- which is what makes this robust to the
        quality gradient that derails hard clustering.

        Returns (scores, coverage): samples x program, and how many markers of
        each program were found.
        """

        # --- identifier vocabulary -------------------------------------------
        # PROGRAMS are curated as gene symbols, but the matrices are ensembl-
        # indexed once build_bulk_matrix(gene_key="geneid") is used. Translate
        # on the fly rather than making the caller remember.
        all_genes = [g for progs in self.PROGRAMS.values()
                     for genes in progs.values() for g in genes]
        if not all_genes:
            raise ValueError(
                "PROGRAMS contains no genes. If a previous call overwrote it "
                "with an empty translation, restore from the class: "
                "mc.PROGRAMS = MalignantCluster.PROGRAMS")
        # Majority vote over every gene, not the first one of an arbitrary
        # compartment: dict order is not meaningful and one stray entry should
        # not decide the vocabulary for all of them.
        n_ensg = sum(g.startswith("ENSG") for g in all_genes)
        translate = n_ensg < len(all_genes) / 2
        if 0 < n_ensg < len(all_genes):
            warnings.warn(
                f"PROGRAMS mixes vocabularies: {n_ensg}/{len(all_genes)} look "
                "like ensembl ids. Translating the symbols only.")

        programs = self.PROGRAMS
        sym_missing = {}
        if translate:
            gene_map = self.prism.load_gene_map("geneid")
            sym2id = {s: i for i, s in gene_map["symbol"].astype(str).items()}

            # NOT written back to self.PROGRAMS: that made the transform
            # destructive -- a translation that silently produced empty lists
            # became its own input on the next call and could not be recovered.
            programs = {
                comp: {prog: [sym2id.get(g, g) if not g.startswith("ENSG") else g
                              for g in genes
                              if g.startswith("ENSG") or g in sym2id]
                       for prog, genes in progs.items()}
                for comp, progs in self.PROGRAMS.items()
            }
            sym_missing = {f"{c}.{p}": [g for g in gs
                                        if not g.startswith("ENSG") and g not in sym2id]
                           for c, ps in self.PROGRAMS.items() for p, gs in ps.items()}
            n_lost = sum(len(v) for v in sym_missing.values())
            if n_lost:
                warnings.warn(
                    f"{n_lost} marker symbol(s) have no ensembl id and were "
                    "dropped before scoring; see the 'missing_symbol' column "
                    "of the coverage table.")
            self.PROGRAMS_ENSG = programs        # keep it, but do not clobber

        if compartment_map is None:
            compartment_map = {"malignant": self.mal_cell_name}

        out, cov = {}, []
        for key, cell_type in compartment_map.items():
            if key not in programs:
                warnings.warn(f"no PROGRAMS entry for '{key}'; skipping")
                continue

            M = self.compartment_matrix(cell_type, min_share=min_share, min_counts=min_counts, batch=batch, samples=samples, **mat_kw)

            Zs = (M - M.mean()) / M.std().replace(0, np.nan)

            for prog, genes in programs[key].items():
                found = [g for g in genes if g in Zs.columns]

                cov.append({"compartment": key, "cell_type": cell_type,
                            "program": prog, "n_found": len(found),
                            "n_total": len(genes),
                            "missing": [g for g in genes if g not in Zs.columns],
                            "missing_symbol": sym_missing.get(f"{key}.{prog}", [])})
                if len(found) < min_genes:
                    warnings.warn(
                        f"{key}/{prog}: only {len(found)}/{len(genes)} markers "
                        f"present; skipping (min_genes={min_genes})")
                    continue
                out[f"{key}.{prog}"] = Zs[found].mean(axis=1)

        scores = pd.DataFrame(out)

        # A program score that tracks its own compartment's abundance is an
        # abundance readout, not a phenotype.
        for rec in cov:
            col = f"{rec['compartment']}.{rec['program']}"
            if col not in scores:
                continue
            th = self.df_theta[rec["cell_type"]].reindex(scores.index)
            ok = scores[col].notna() & th.notna()
            rec["r_with_theta"] = (round(float(np.corrcoef(
                scores.loc[ok, col], th[ok])[0, 1]), 3) if ok.sum() > 5 else np.nan)
            
        # Signed axes are more stable than either pole alone. Name them with
        # their OWN compartment prefix so couple_compartments() still treats
        # them as belonging to it -- an "axis." prefix would make every axis
        # look like a separate compartment and produce trivial self-correlations
        # against its own poles.
        if "malignant.basal" in scores and "malignant.classical" in scores:
            scores["malignant.axis_basal_minus_classical"] = (
                scores["malignant.basal"] - scores["malignant.classical"])
        if "fibroblast.myCAF" in scores and "fibroblast.iCAF" in scores:
            scores["fibroblast.axis_myCAF_minus_iCAF"] = (
                scores["fibroblast.myCAF"] - scores["fibroblast.iCAF"])
        if "immune.cytotoxic" in scores and "immune.exhaustion" in scores:
            scores["immune.axis_cytotoxic_minus_exhaustion"] = (
                scores["immune.cytotoxic"] - scores["immune.exhaustion"])
        return scores, pd.DataFrame(cov)


    def compare_coupled_compartments(
        self,
        program_list: list,
        coh: pd.Series,
        cmap: dict,
        scores: pd.DataFrame,
        n_perm: int = 2000,
        control_theta: bool = True,
        force: bool = False,
        verbose: bool = False,
    ) -> pd.DataFrame:

        if not isinstance(program_list, list) or len(program_list) != 2:
            print("Error: program_list should be like ['TCGA', 'CPTAC']")
            return pd.DataFrame()

        program_a = program_list[0]
        program_b = program_list[1]

        fname = f"corr_coupled_{program_a}_x_{program_b}_n_perm={n_perm}_ctr_theta_{control_theta}.tsv"
        fname = title_replace(fname)
        filename = self.root_prism / fname

        if filename.exists() and not force:
            return pdreadcsv(fname, self.root_prism, verbose=verbose)

        from scipy.stats import norm
        from statsmodels.stats.multitest import multipletests

        full = {}
        for name in program_list:
            print(name)
            idx = list(scores.index[coh == name])
            sc_c, _ = self.program_scores(compartment_map=cmap, samples=idx)

            fname_corr=f"prog_corr_scores_x_malig_for_{name}_n{len(sc_c.index)}"
            full[name] = self.couple_compartments(scores=sc_c, fname_corr=fname_corr, cell_name=self.mal_cell_name, n_perm=n_perm, force=False, verbose=verbose)

        print("----------- end --------------")
        key = ["program_a", "program_b"]

        df_meta = full[program_a].merge(full[program_b], on=key, suffixes=("_t", "_c"))

        z1, z2 = np.arctanh(df_meta.r_t), np.arctanh(df_meta.r_c)
        w1, w2 = df_meta.n_t - 3, df_meta.n_c - 3

        se = np.sqrt(1/w1 + 1/w2)
        df_meta["z_het"] = (z1 - z2) / se
        df_meta["p_het"] = 2 * norm.sf(np.abs(df_meta.z_het))

        df_meta["z_meta"]   = (z1*w1 + z2*w2) / (w1 + w2)
        df_meta["r_meta"]   = np.tanh(df_meta.z_meta)                       # identical to before
        df_meta["p_meta"]   = 2 * norm.sf(np.abs(df_meta.z_meta) * np.sqrt(w1 + w2))
        df_meta["fdr_meta"] = multipletests(df_meta.p_meta, method="fdr_bh")[1]

        Q = w1*(z1 - df_meta.z_meta)**2 + w2*(z2 - df_meta.z_meta)**2       # Cochran's Q, ~chi2_1
        df_meta["I2"] = np.clip((Q - 1) / Q, 0, 1)

        df_meta["concordant"] = np.sign(df_meta.r_t) == np.sign(df_meta.r_c)

        _ = pdwritecsv(df_meta, fname, self.root_prism, verbose=verbose)    

        return df_meta
    
    
    def couple_compartments(
        self,
        scores: pd.DataFrame,
        fname_corr: str,
        cell_name: str,
        n_perm: int = 2000,
        method: str = "spearman",
        control_theta: bool = True,
        random_state: int = 0,
        force: bool = False,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """
        Test coupling between programs in DIFFERENT compartments.

        Two confounds make the naive correlation untrustworthy, and both are handled here:

        1. Compositional coupling. Z_mal and Z_fib are constrained to sum to the
           observed bulk, so they are negatively coupled by construction.
           Compartment-internal CPM removes most of the scale component; setting
           `control_theta=True` additionally partials out theta_mal, which is
           the main remaining driver.
        2. Shared reference. Both compartments come from one BayesPrism fit
           against one scRNA-seq reference, so a misspecified reference can
           induce apparent coupling. The permutation null (shuffling samples in
           one compartment only) preserves each program's marginal distribution
           while destroying the sample pairing, which is the correct null for
           "do these co-vary ACROSS tumours".

        Returns pairs with r, permutation p, and BH FDR. Cross-compartment
        pairs only -- within-compartment correlations are not interesting here
        and would dominate the FDR.
        """

        fname = fname_corr + f"_n_perm={n_perm}_ctr_theta_{control_theta}.tsv"
        fname = title_replace(fname)
        filename = self.root_prism / fname

        if filename.exists() and not force:
            return pdreadcsv(fname, self.root_prism, verbose=verbose)

        rng = np.random.default_rng(random_state)
        S = scores.dropna(axis=1, how="all").copy()

        if control_theta:
            df_th = self.df_theta[cell_name].reindex(S.index)
            D = np.column_stack([np.ones(len(S)), df_th.fillna(df_th.mean()).values])
            beta, *_ = np.linalg.lstsq(D, S.fillna(S.mean()).values, rcond=None)
            S = pd.DataFrame(S.values - D @ beta, index=S.index, columns=S.columns)

        def _corr(x, y):
            if method == "spearman":
                x, y = stats.rankdata(x), stats.rankdata(y)
            return float(np.corrcoef(x, y)[0, 1])

        cols = list(S.columns)
        rows = []
        for i, a in enumerate(cols):
            for b in cols[i + 1:]:
                ca, cb = a.split(".")[0], b.split(".")[0]
                if ca == cb:
                    continue                      # cross-compartment only
                x = S[a].values
                y = S[b].values
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 10:
                    continue
                r = _corr(x[ok], y[ok])
                null = np.empty(n_perm)
                for p_ in range(n_perm):
                    # correlation good x with permutation good y
                    null[p_] = _corr(x[ok], rng.permutation(y[ok]))
                pval = (np.sum(np.abs(null) >= abs(r)) + 1) / (n_perm + 1)
                rows.append({"program_a": a, "program_b": b, "r": round(r, 4),
                             "p_perm": pval, "n": int(ok.sum())})

        df_corr = pd.DataFrame(rows)
        if len(df_corr):
            df_corr["fdr"] = self._bh(df_corr["p_perm"])
            df_corr = df_corr.sort_values("p_perm").reset_index(drop=True)

        pdwritecsv(df_corr, fname, self.root_prism, verbose=verbose)
                    
        return df_corr


    def program_genes(self, compartments: Optional[Sequence[str]] = None) -> List[str]:
        """Every marker gene used to define program scores.

        Exclude these before running DE on program-derived labels: the groups
        were defined by these genes, so they top the table by construction and
        tell you nothing. What is informative is which OTHER genes track the
        labels.
        """

        if self.program_gene_list:
            return self.program_gene_list
    
        keys = compartments or list(self.PROGRAMS)
        out = set()
        for k in keys:
            for gs in self.PROGRAMS.get(k, {}).values():
                out.update(gs)

        self.program_gene_list = list(out)
        return list(self.program_gene_list )

    def state_signatures(
        self,
        X: pd.DataFrame,
        joint_label: pd.Series,
        exclude_program_genes: bool = True,
        min_per_group: int = 10,
        verbose: bool = False,
        **sig_kw,
    ) -> pd.DataFrame:
        """DE table for program-defined joint states, marker genes removed.

        Groups here come from thresholding marker scores, so leaving those
        markers in the feature set guarantees they rank first and proves
        nothing. Dropping them turns the question into the useful one: what
        else distinguishes these states?

        Note this does not escape the circularity that afflicts cluster DE
        generally -- the labels still came from this cohort's own expression --
        but unlike consensus clustering, the thresholds are fixed a priori by
        marker biology rather than chosen to maximise separation.
        """
        lab = joint_label.reindex(X.index).dropna()
        counts = lab.value_counts()
        keep = counts[counts >= min_per_group].index
        if len(keep) < 2:
            raise ValueError(
                f"only {len(keep)} state(s) with >= {min_per_group} samples: "
                f"{counts.to_dict()}")
        if len(keep) < len(counts):
            warnings.warn(f"dropping sparse states: "
                          f"{counts[~counts.index.isin(keep)].to_dict()}")
        lab = lab[lab.isin(keep)]

        Xs = X.loc[lab.index]
        if exclude_program_genes:
            drop = set(self.program_genes())
            n0 = Xs.shape[1]
            Xs = Xs.loc[:, [g for g in Xs.columns if g not in drop]]
            if verbose: print(f"excluded {n0 - Xs.shape[1]} program marker genes; {Xs.shape[1]} remain")

        codes = pd.Series(pd.Categorical(lab).codes, index=lab.index)
        mapping = dict(enumerate(pd.Categorical(lab).categories))
        T = self.signature_table(Xs, codes, **sig_kw)
        T["state"] = T["cluster"].map(mapping)
        return T


    def factorial_state_de(
        self,
        X: pd.DataFrame,
        disc: pd.DataFrame,
        axis_a: Optional[str] = None,
        axis_b: Optional[str] = None,
        exclude_program_genes: bool = True,
        high_label: str = "high",
        verbose: bool = False,
    ) -> pd.DataFrame:
        """Two-way factorial DE for a 2x2 state design.

        Fits, per gene:  expr ~ 1 + A + B + A:B   (A, B coded 0/1)

        - Axis scores (continuous, mean-centred):
        - A_c = basal_minus_classical   (malignant compartment)
        - B_c = myCAF_minus_iCAF        (fibroblast compartment)

        **expr ~ 1 + A_c + B_c + A_c:B_c + cohort**

        - Parameters:
        - beta_A   tumour-axis slope, evaluated at mean stroma score
        - beta_B   stroma-axis slope, evaluated at mean tumour score
        - beta_AB  interaction: change in the A slope per unit B
                    (candidate crosstalk — but compositional/θ effects
                    produce interactions too, so this is not evidence
                    of signalling on its own)

        One-vs-rest on four cells cannot separate these: a gene driven purely by
        the tumour axis appears "significant" in two of the four contrasts and
        looks like a state marker when it is not.

        Returns one row per gene with all three coefficients, their p-values and
        BH FDR computed within each term.
        """
        cols = list(disc.columns)
        axis_a = axis_a or cols[0]
        axis_b = axis_b or cols[1]
        if len(cols) < 2:
            raise ValueError(f"need >=2 discretised axes, got {cols}")

        d = disc[[axis_a, axis_b]].dropna()
        # common samples
        idx = X.index.intersection(d.index)
        d, Xs = d.loc[idx], X.loc[idx]

        if exclude_program_genes:
            drop = set(self.program_genes())
            n0 = Xs.shape[1]
            Xs = Xs.loc[:, [g for g in Xs.columns if g not in drop]]
            if verbose:
                print(f"excluded {n0 - Xs.shape[1]} program marker genes;  {Xs.shape[1]} remain")

        A = (d[axis_a] == high_label).astype(float).values
        B = (d[axis_b] == high_label).astype(float).values
        D = np.column_stack([np.ones(len(A)), A, B, A * B])
        n, p = D.shape
        if n <= p:
            raise ValueError(f"{n} samples cannot fit {p} coefficients")

        Y = Xs.values
        beta, *_ = np.linalg.lstsq(D, Y, rcond=None)
        resid = Y - D @ beta
        dof = n - p
        sigma2 = (resid ** 2).sum(axis=0) / dof
        XtX_inv = np.linalg.pinv(D.T @ D)
        se = np.sqrt(np.outer(np.diag(XtX_inv), sigma2))       # 4 x n_genes

        names = ["intercept", f"A_{axis_a}", f"B_{axis_b}", "interaction"]
        out = {"gene": Xs.columns}
        for i, nm in enumerate(names):
            if nm == "intercept":
                continue
            with np.errstate(invalid="ignore", divide="ignore"):
                t = beta[i] / se[i]
            pv = pd.Series(2 * stats.t.sf(np.abs(t), dof),
                           index=Xs.columns).fillna(1.0)
            out[f"beta_{nm}"] = beta[i]
            out[f"p_{nm}"] = pv.values
            out[f"fdr_{nm}"] = self._bh(pv).values

        dfr = pd.DataFrame(out).set_index("gene")
        dfr["n_samples"] = n
        dfr["cell_counts"] = str(d.groupby(list(d.columns)).size().to_dict())
        dfr = dfr.sort_values("fdr_interaction")
        return dfr

    # ===========================================================================
    # 3d. Joint multi-compartment states (instead of a Cartesian product)
    # ===========================================================================
    @staticmethod
    def state_feasibility(n_samples: int, n_axes: int, levels: int = 2,
                          min_per_cell: int = 15) -> pd.DataFrame:
        """How many joint states your n can actually support.

        A Cartesian product of per-compartment clusterings grows as
        levels**n_axes while n stays fixed. Run this BEFORE choosing how many
        axes to cross: it is the cheapest way to see that 3 axes x 3 levels at
        n=110 means ~4 samples per cell.
        """
        rows = []
        for a in range(1, n_axes + 1):
            for L in range(2, levels + 1):
                cells = L ** a
                rows.append({"n_axes": a, "levels": L, "cells": cells,
                             "mean_per_cell": round(n_samples / cells, 1),
                             "usable": n_samples / cells >= min_per_cell})
        return pd.DataFrame(rows)

    def discretize_axes(
        self,
        scores: pd.DataFrame,
        axes: Optional[Sequence[str]] = None,
        method: str = "median",
        min_frac: float = 0.20,
    ) -> pd.DataFrame:
        """Cut continuous program axes into high/low labels.

        Two levels only, deliberately. Three-level cuts on n~110 leave cells too
        small to interpret, and the middle level of a continuous axis is the
        least meaningful part of it.

        method="median" splits at the median (balanced by construction);
        method="gmm" fits a 2-component Gaussian mixture and warns if the
        smaller component falls below `min_frac`, which indicates the axis is
        unimodal and should stay continuous.
        """
        axes = axes or [c for c in scores.columns if ".axis_" in c]
        if not axes:
            raise ValueError(f"no axis columns found in {list(scores.columns)}")

        out = {}
        for a in axes:
            v = scores[a].dropna()
            if method == "median":
                lab = (v > v.median()).map({True: "high", False: "low"})
            elif method == "gmm":
                from sklearn.mixture import GaussianMixture
                g = GaussianMixture(2, random_state=0).fit(v.values.reshape(-1, 1))
                comp = pd.Series(g.predict(v.values.reshape(-1, 1)), index=v.index)
                hi = comp.map(dict(enumerate(g.means_.ravel()))).idxmax()
                hi_comp = comp.loc[hi]
                lab = comp.map(lambda x: "high" if x == hi_comp else "low")
                frac = min((lab == "high").mean(), (lab == "low").mean())
                if frac < min_frac:
                    warnings.warn(
                        f"{a}: smaller GMM component is {frac:.0%} of samples "
                        "-- the axis looks unimodal, so a hard cut invents a "
                        "group. Prefer the continuous score.")
            else:
                raise ValueError("method must be 'median' or 'gmm'")
            out[a.split(".")[-1].replace("axis_", "")] = lab
        return pd.DataFrame(out)

    def axis_crosstab(self, disc: pd.DataFrame, min_per_cell: int = 10) -> Dict:
        """Cross-tabulate discretised axes and TEST independence.

        The count of occupied states is an empirical result, not `levels**n_axes`.
        Correlated axes (e.g. basal tumour with myCAF stroma) concentrate samples
        on the diagonal, so the realised number of states is smaller -- often
        much smaller -- than the product. This reports what you actually have.
        """
        cols = list(disc.columns)
        joint = disc.apply(lambda r: " | ".join(f"{c}:{r[c]}" for c in cols), axis=1)
        counts = joint.value_counts()

        out = {"joint_label": joint,
               "counts": counts,
               "n_cells_possible": int(np.prod([disc[c].nunique() for c in cols])),
               "n_cells_occupied": int((counts > 0).sum()),
               "n_cells_usable": int((counts >= min_per_cell).sum()),
               "sparse_cells": counts[counts < min_per_cell].to_dict()}

        if len(cols) == 2:
            from scipy.stats import chi2_contingency
            ct = pd.crosstab(disc[cols[0]], disc[cols[1]])
            chi2, pval, dof, exp = chi2_contingency(ct)
            n = ct.values.sum()
            out["crosstab"] = ct
            out["chi2_p"] = float(pval)
            out["cramers_v"] = float(np.sqrt(chi2 / (n * (min(ct.shape) - 1))))
            out["independence_note"] = (
                "axes are NOT independent -- the product overcounts states"
                if pval < 0.05 else
                "axes look independent -- the product is a reasonable upper bound")
        else:
            out["independence_note"] = (
                "pairwise independence not tested for >2 axes; inspect `counts` "
                "for diagonal concentration")
        return out



    @staticmethod
    def harmonize_scores(scores: pd.DataFrame, batch: pd.Series,
                         scale: bool = True, report: bool = True) -> pd.DataFrame:
        """Standardise program scores within batch (location and, optionally, scale).

        Cohort effects on deconvolved scores are often heteroscedastic rather
        than a mean shift: two cohorts can have identical means while one has
        1.7x the SD. Mean-centring leaves that untouched, and every extreme --
        hence every tail cluster -- then comes from the wider cohort.

        `scale=True` z-scores within batch, equalising both. Use it when the
        extra dispersion is technical. Be careful if it might be biological: a
        cohort with genuinely more heterogeneous tumours will have its real
        spread flattened too, and with no overlapping samples the two cannot be
        distinguished. Reporting both versions is the honest option.
        """
        from scipy.stats import levene
        batch = batch.reindex(scores.index)
        out = scores.copy()
        rows = []
        for c in scores.columns:
            g = scores[c].groupby(batch)
            mu, sd = g.transform("mean"), g.transform("std")
            out[c] = (scores[c] - mu) / (sd if scale else 1.0)
            if report:
                grp = [scores.loc[batch == b, c].dropna() for b in batch.unique()]
                grp = [x for x in grp if len(x) > 2]
                rows.append({
                    "score": c,
                    "sd_ratio": round(float(max(x.std() for x in grp) /
                                            max(min(x.std() for x in grp), 1e-12)), 2),
                    "mean_range": round(float(max(x.mean() for x in grp) -
                                              min(x.mean() for x in grp)), 3),
                    "levene_p": round(float(levene(*grp).pvalue), 5) if len(grp) > 1 else np.nan,
                })
        if report and rows:
            print(pd.DataFrame(rows).to_string(index=False))
        return out

    def axis_modality(self, scores: pd.DataFrame,
                      axes: Optional[Sequence[str]] = None,
                      random_state: int = 0) -> pd.DataFrame:
        """Is each axis one cloud or two? BIC of a 1- vs 2-component Gaussian.

        Decisive for interpretation. If an axis is unimodal, states carved from
        it are tails of a continuum, not phenotypes: valid as a stratification,
        but you cannot call them subtypes, and a sample near the boundary has no
        stable assignment.

        delta_bic = BIC(1 comp) - BIC(2 comp). Positive favours two components;
        the usual reading is >10 strong, 2-10 weak, <2 negligible. `min_weight`
        guards the case where the 2-component fit wins only by devoting a tiny
        component to outliers.
        """
        from sklearn.mixture import GaussianMixture

        axes = axes or [c for c in scores.columns if ".axis_" in c] or list(scores.columns)
        rows = []
        for a in axes:
            v = scores[a].dropna().values.reshape(-1, 1)
            if len(v) < 20:
                rows.append({"axis": a, "n": len(v), "verdict": "too few samples"})
                continue
            g1 = GaussianMixture(1, random_state=random_state).fit(v)
            g2 = GaussianMixture(2, random_state=random_state, n_init=5).fit(v)
            d = g1.bic(v) - g2.bic(v)
            w = float(min(g2.weights_))
            sep = float(abs(np.diff(g2.means_.ravel()))[0] /
                        np.sqrt(g2.covariances_.ravel().mean()))
            verdict = ("bimodal" if d > 10 and w >= 0.15 and sep > 1.5 else
                       "weak/outlier-driven" if d > 2 else "unimodal")
            rows.append({"axis": a, "n": len(v), "delta_bic": round(d, 1),
                         "min_component_weight": round(w, 3),
                         "separation_sd": round(sep, 2),
                         "skew": round(float(stats.skew(v.ravel())), 2),
                         "verdict": verdict})
        return pd.DataFrame(rows)

    def joint_states(
        self,
        scores: pd.DataFrame,
        ks: Iterable[int] = range(2, 7),
        use: Optional[Sequence[str]] = None,
        **cc_kw,
    ) -> Dict:
        """Cluster samples on the multi-compartment PROGRAM SCORE matrix.

        This is the alternative to a Cartesian product. Instead of clustering
        each compartment separately and crossing the labels -- which multiplies
        groups, multiplies the ways clustering can fail, and assumes the axes are
        independent -- it finds states in the joint space directly. Occupied
        states emerge from the data at whatever number the samples support.

        The feature space is a handful of interpretable program scores rather
        than thousands of genes, so it is far less vulnerable to the quality
        gradient that derails gene-level clustering: a globally depressed sample
        scores near zero on every axis rather than defining its own cluster.
        """
        use = use or [c for c in scores.columns if ".axis_" in c] or list(scores.columns)
        M = scores[use].dropna()
        if M.shape[0] < 10:
            raise ValueError(f"only {M.shape[0]} samples with complete scores")
        Ms = (M - M.mean()) / M.std().replace(0, np.nan)
        Ms = Ms.dropna(axis=1, how="all")

        # Feature subsampling is meant for thousands of genes. On a handful of
        # program axes, dropping 20% of features means clustering on a random
        # 2-of-3 each resample, which inflates PAC and understates stability.
        cc_kw.setdefault("feature_frac", 1.0)
        if Ms.shape[1] <= 10 and cc_kw["feature_frac"] < 1.0:
            warnings.warn(
                f"feature_frac={cc_kw['feature_frac']} on only {Ms.shape[1]} "
                "axes will inflate PAC; use 1.0.")

        cc = self.consensus_cluster(Ms, ks=ks, metric="euclidean",
                                    method="ward", **cc_kw)
        return {"features": list(Ms.columns), "matrix": Ms, "consensus": cc,
                "summary": self.cluster_summary(cc)}

    def state_profile(self, js: Dict, k: Optional[int] = None,
                      min_cluster_frac: float = 0.10) -> pd.DataFrame:
        """Mean axis score per joint state -- how you name the states.

        A state is only interpretable if its axes separate: e.g. one state
        basal-high/myCAF-high, another classical-high/iCAF-high. If every state
        differs mainly in overall magnitude rather than in the pattern across
        axes, you are looking at a signal-strength gradient, not distinct
        phenotypes.
        """
        cc = js["consensus"]
        k = k or self.choose_k(cc, min_cluster_frac=min_cluster_frac)
        lab = cc[k]["labels"]
        M = js["matrix"]
        prof = M.groupby(lab.reindex(M.index)).mean().round(3)
        prof["n"] = lab.value_counts().sort_index()
        prof["mean_abs_score"] = M.abs().groupby(
            lab.reindex(M.index)).mean().mean(axis=1).round(3)
        return prof

    def load_tahoe_de(
        self,
        organs: Optional[Sequence[str]] = None,
        cell_lines: Optional[Sequence[str]] = None,
        drugs: Optional[Sequence[str]] = None,
        genes: Optional[Sequence[str]] = None,
        mode: str = "stream",
        save_to=None,
        force: bool = False,
        memory_limit: str = "8GB",
        temp_dir=None,
        threads: Optional[int] = 4,
        hf_token: Optional[str] = None,
        http_retries: int = 8,
        http_retry_wait_s: float = 5.0,
        checkpoint_dir=None,
        verbose: bool = True,
        _shard_limit: Optional[int] = None,
        _shard_subset: Optional[Sequence[str]] = None,
        dtype=np.float32,
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

        Parameters
        ----------
        mode
            "stream" (default) queries the DE parquet shards directly over
            `hf://` with DuckDB and only materialises the filtered subset --
            nothing large touches disk. "download" fetches all 21 shards into
            `self.root_tahoe` first, which is worth it only if you will re-query many
            times offline.
        save_to
            Directory to persist the derived `S` / `cond` as parquet. Defaults
            to `<root_tahoe>/derived`, so `run()` and ad-hoc
            calls share one cache. Files are keyed by a hash of
            organs/cell_lines/drugs/genes, so changing any filter re-queries
            rather than returning a stale matrix. Pass `save_to=False` to
            disable caching entirely.

        organs
            Defaults to `self.organ` (set on the constructor). Must match
            Tahoe's `Organ` vocabulary exactly; a mismatch now lists the valid
            values instead of reporting an empty result.
        genes
            Restrict the gene universe. **Pass this.** Without it the pancreas
            subset alone is ~5e8 rows (~6 lines x 379 drugs x ~4 doses x ~54k
            genes), which is >100 GB as an object-dtype pandas frame before
            pivoting -- the usual cause of an OOM here. Passing `X.columns`
            (your HVG set) cuts it by 25-30x.

            Restricting the universe does change the enrichment scores: ES is
            computed over the ranked list you supply, so a 2k-gene universe is
            a different null than a 54k one. This is the same trade CMap makes
            with its 978 landmark genes, and it is defensible provided you say
            so -- but do not mix scores computed over different universes.
        memory_limit, temp_dir
            Passed to DuckDB. With a temp dir set, DuckDB spills to disk rather
            than being OOM-killed.

        Notes
        -----
        Aggregation and pivoting happen per cell line inside DuckDB; only the
        genes x conditions block for one line is ever in memory. Column names
        have shifted between revisions of the dataset card, so the schema is
        sniffed at runtime; override with `stat_col` / `gene_col` if the
        sniffing picks the wrong one.
        """

        # DuckDB - https://duckdb.org/docs/current/guides/overview
        import duckdb

        if mode not in ("stream", "download"):
            raise ValueError(f"mode must be 'stream' or 'download', got {mode!r}")

        if organs is None and not cell_lines:
            organs = (self.organ,)

        # --- reuse a previously saved result --------------------------------
        # save_to=False disables caching; None means "use the default dir".
        if save_to is not False:
            sd = self.tahoe_derived_dir(save_to)
            key = self._tahoe_cache_key(organs, cell_lines, drugs, genes)
            f_S = sd / f"tahoe_S_{key}.parquet"
            f_c = sd / f"tahoe_cond_{key}.parquet"
            if f_S.exists() and f_c.exists() and not force:
                return pd.read_parquet(f_S), pd.read_parquet(f_c)

        # --- small tables always local (a few MB) ---------------------------
        need = [self.root_tahoe / f for f in self._TAHOE_SMALL[:2]]
        if force or not all(f.exists() for f in need):
            self.download_tahoe(de=False, force=force)

        con = duckdb.connect()
        con.execute(f"SET memory_limit='{memory_limit}'")
        con.execute("SET preserve_insertion_order=false")
        if temp_dir is not None:
            create_dir(temp_dir)
            con.execute(f"SET temp_directory='{temp_dir}'")

        if threads:
            con.execute(f"SET threads={int(threads)}")

        if mode == "stream":
            con.execute("INSTALL httpfs; LOAD httpfs;")
            # Anonymous HF requests are rate-limited far more aggressively than
            # authenticated ones; a token is the single most effective fix for
            # HTTP 429 on this table.
            tok = hf_token or os.environ.get("HF_TOKEN") or os.environ.get(
                "HUGGING_FACE_HUB_TOKEN")
            if tok:
                con.execute(
                    f"CREATE OR REPLACE SECRET hf_tok "
                    f"(TYPE huggingface, TOKEN '{tok}')")
            else:
                warnings.warn(
                    "No HF_TOKEN found. Streaming ~1026 remote parquet shards "
                    "anonymously will likely hit HTTP 429. Set HF_TOKEN, or "
                    "use mode='download'.")
            de_glob = ("'hf://datasets/" + self._TAHOE_REPO +
                       "/metadata/pseudobulk_differential_expression/*.parquet'")
        elif mode == "download":
            de_dir = self.root_tahoe / "metadata" / "pseudobulk_differential_expression"
            if _shard_subset:
                local = [self.root_tahoe / f for f in _shard_subset]
                missing = [f for f in local if not f.exists()]
                if missing:
                    raise FileNotFoundError(
                        f"{len(missing)} of {len(local)} shards not downloaded "
                        f"(e.g. {missing[0]}). Run download_shards(shards) first.")
                de_glob = "', '".join(str(f) for f in local)
                de_glob = f"['{de_glob}']"
            else:
                if force or not any(de_dir.glob("*.parquet")):
                    self.download_tahoe(de=True, force=force)
                de_glob = f"'{de_dir / '*.parquet'}'"

        cl_meta = con.execute(
            f"SELECT * FROM read_parquet('{self.root_tahoe / 'metadata' / 'cell_line_metadata.parquet'}')"
        ).df()
        drug_meta = con.execute(
            f"SELECT * FROM read_parquet('{self.root_tahoe /'metadata'/'drug_metadata.parquet'}')"
        ).df()

        # --- resolve which cell lines we want ----------------------------------
        sel = cl_meta.copy()
        if organs:
            sel = sel[sel["Organ"].isin(organs)]
        if cell_lines:
            # Accept either Cellosaurus ids or cell names -- find_shards_for and
            # organ_cell_lines both speak CVCL, so requiring names here made the
            # natural call fail with a misleading "no cell lines match".
            want = {str(x) for x in cell_lines}
            m = sel["Cell_ID_Cellosaur"].astype(str).isin(want)
            if "cell_name" in sel.columns:
                m = m | sel["cell_name"].astype(str).isin(want)
            matched = set(sel.loc[m, "Cell_ID_Cellosaur"].astype(str))
            if "cell_name" in sel.columns:
                matched |= set(sel.loc[m, "cell_name"].astype(str))
            unknown = want - matched
            if unknown:
                warnings.warn(
                    f"cell_lines not found in cell_line_metadata: "
                    f"{sorted(unknown)}. Note the metadata catalogues ~102 "
                    "lines but only ~50 were profiled, so presence here still "
                    "does not guarantee rows in the DE table.")
            sel = sel[m]
        if sel.empty:
            raise ValueError(
                f"no cell lines match organs={organs} lines={cell_lines}.\n"
                f"cell_lines accepts Cell_ID_Cellosaur ids or cell_name values.\n"
                f"Organ vocabulary (exact, case-sensitive): "
                f"{sorted(cl_meta['Organ'].dropna().unique())}")
        cvcl = sorted(sel["Cell_ID_Cellosaur"].dropna().unique().tolist())

        # --- sniff schema -------------------------------------------------------
        cols = con.execute(
            f"SELECT * FROM read_parquet({de_glob}) LIMIT 1"
        ).df().columns.tolist()

        gcol = self._pick_col(cols, self._DE_GENE_CANDS, gene_col, "gene")
        scol = self._pick_col(cols, self._DE_STAT_CANDS, stat_col, "statistic")
        ccol = self._pick_col(cols, self._DE_CELL_CANDS, None, "cell line")
        dcol = self._pick_col(cols, self._DE_DRUG_CANDS, None, "drug")

        # --- per-cell-line fetch: one genes x conditions block at a time ----
        gene_filter = ""
        if genes is not None:
            g_list = sorted({str(g) for g in genes})
            con.execute("CREATE TEMP TABLE _qgenes(gene VARCHAR)")
            con.register("_qg_df", pd.DataFrame({"gene": g_list}))
            con.execute("INSERT INTO _qgenes SELECT gene FROM _qg_df")
            gene_filter = f" AND {gcol} IN (SELECT gene FROM _qgenes)"
        else:
            warnings.warn(
                "load_tahoe_de(genes=None): pulling the full ~54k-gene universe. "
                "For the pancreas subset this is ~5e8 rows and will exhaust a "
                "64 GB machine. Pass genes=X.columns.")

        drug_filter = ""
        if drugs:
            drug_filter = f" AND {dcol} IN ({','.join(repr(str(x)) for x in drugs)})"

        # One scan of the shard set for every cell line at once. The DE table
        # is ~1026 remote parquet shards; looping per cell line re-scans all of
        # them each time and gets you rate-limited (HTTP 429). With `genes`
        # restricted the combined result is small enough to pivot in pandas.
        cl_in = ",".join(repr(str(x)) for x in cvcl)
        select_sql = (f"SELECT {ccol} AS cell_line_id, {dcol} AS drug, "
                      f"{gcol} AS gene, avg({scol}) AS stat")
        where_sql = (f"WHERE {ccol} IN ({cl_in}){gene_filter}{drug_filter} "
                     f"GROUP BY 1, 2, 3")

        if mode == "stream":
            shards = _shard_subset or self.list_de_shards()
            if _shard_limit:
                shards = shards[:_shard_limit]
            ckdir = Path(checkpoint_dir) if checkpoint_dir is not None \
                else self.root_tahoe / "_scan_ck" / \
                     self._tahoe_cache_key(organs, cell_lines, drugs, genes)
            if verbose:
                print(f"scanning {len(shards)} DE shards -> {ckdir}", flush=True)
            df = self._scan_de_sharded(con, shards, select_sql, where_sql,
                                       ckdir, http_retries, http_retry_wait_s,
                                       verbose=verbose)
            if not df.empty:
                df = (df.groupby(["cell_line_id", "drug", "gene"], observed=True,
                                 as_index=False)["stat"].mean())
        else:
            q = f"{select_sql} FROM read_parquet({de_glob}) {where_sql}"
            df = self._duck_retry(con, q, retries=http_retries,
                                  wait_s=http_retry_wait_s)

        if df.empty and (_shard_limit or _shard_subset):
            warnings.warn(
                f"0 rows from the {_shard_limit or len(_shard_subset)} sampled "
                "shard(s). Expected when "
                "shards are clustered by cell line -- probe more shards, or "
                "use find_shards_for() to locate the ones holding your lines.")
            empty = pd.DataFrame(dtype=dtype)
            return empty, pd.DataFrame(columns=["cell_line_id", "drug"])

        if df.empty:
            raise ValueError(
                "Tahoe DE query returned 0 rows. This is almost always an "
                "identifier-vocabulary mismatch, not an absent organ:\n"
                f"  cell-line filter used {ccol} IN {cvcl[:3]}... "
                f"({len(cvcl)} values from cell_line_metadata)\n"
                f"  gene filter used {gcol} IN <"
                f"{0 if genes is None else len(set(map(str, genes)))} symbols>\n"
                "Run inspect_de_schema() to see what those columns actually "
                "contain (CVCL ids vs cell names; HGNC symbols vs Ensembl ids "
                "vs integer token ids).")

        frac_neg = float((df["stat"] < 0).mean())
        if not (0.05 < frac_neg < 0.95):
            warnings.warn(
                f"stat column {scol!r} is {frac_neg:.1%} negative -- it looks "
                "unsigned (a p-value, |t|, or magnitude). WTCS needs a signed "
                "statistic; pass stat_col= with the signed column, e.g. 'stat' "
                "or 'log2FoldChange'.")

        df["stat"] = df["stat"].astype(dtype)
        df["condition"] = df["cell_line_id"].astype(str) + "|" + df["drug"].astype(str)
        S = (df.pivot(index="gene", columns="condition", values="stat")
               .astype(dtype))
        blocks = [S]

        del blocks

        cond = (df[["condition", "cell_line_id", "drug"]]
                .drop_duplicates()
                .merge(cl_meta[["Cell_ID_Cellosaur", "cell_name", "Organ"]],
                    left_on="cell_line_id", right_on="Cell_ID_Cellosaur",
                    how="left")
                .merge(drug_meta[["drug", "moa-fine", "moa-broad", "targets",
                                "human-approved"]],
                    on="drug", how="left")
                .set_index("condition"))
        cond = cond.loc[S.columns]

        if save_to is not False:
            S.to_parquet(f_S)
            cond.to_parquet(f_c)
            (sd / f"tahoe_manifest_{key}.json").write_text(json.dumps({
                "organs": list(organs) if organs else None,
                "cell_lines": list(cell_lines) if cell_lines else None,
                "drugs": list(drugs) if drugs else None,
                "n_genes": None if genes is None else len(set(map(str, genes))),
                "genes_sha1": key,
                "n_conditions": int(S.shape[1]),
                "repo": self._TAHOE_REPO,
            }, indent=2))

        return S, cond


    # --- CMap-style connectivity ------------------------------------------------
    def _weighted_es(self, metric_desc: np.ndarray, hit_mask: np.ndarray,
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


    def _wtcs(self, ref: pd.Series, up: Sequence[str], down: Sequence[str]) -> float:
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
        es_u = self._weighted_es(m, up_m)
        es_d = self._weighted_es(m, dn_m)
        return (es_u - es_d) / 2.0 if np.sign(es_u) != np.sign(es_d) else 0.0


    def score_clusters_vs_tahoe(
        self,
        sig: Dict[int, Dict],
        S: pd.DataFrame,
        cond: pd.DataFrame,
        normalize_within: str = "cell_line_id",
    ) -> pd.DataFrame:
        """
        WTCS (Weighted Connectivity Score) and NCS (Normalized Connectivity Score)
        are core metrics used in the LINCS and Connectivity Map (CMap) pipelines 
        to compare query gene signatures against reference expression profiles. 
        WTCS measures signature similarity from −1 to 1, 
        while NCS normalizes these scores within specific cell lines and perturbagen types.
        
        WTCS + within-cell-line normalised connectivity (NCS) for every (cluster, condition) pair.

        Interpretation: rank by ascending `ncs` to get compounds predicted to
        *reverse* a cluster's malignant program; descending to find compounds whose
        transcriptional footprint *is* that program (mechanistic hypothesis for
        what the cluster's cells are doing).

        Direction convention for Tahoe. 
        score_clusters_vs_tahoe returns WTCS/NCS where 
        - negative = the compound reverses the cluster's program (therapeutic hypothesis) and 
        - positive = the compound's footprint is the program (mechanistic hypothesis for what those cells are doing). 
        
        Pancreas is one of the better-represented organs among Tahoe's 47 usable lines, 
        and the DE query filters at the parquet scan so you never materialise the 4.1B-row table.
        """
        rows = []
        for c, s in sig.items():
            up, dn = s["up"], s["down"]
            for condition in S.columns:
                rows.append({
                    "cluster": c,
                    "condition": condition,
                    "wtcs": self._wtcs(S[condition], up, dn),
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
    # 6. scGPT embedding space
    # ===========================================================================
    
    def scgpt_embed_pseudobulk(
        self,
        X_counts: pd.DataFrame,
        model_dir: str,
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
        self,
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

