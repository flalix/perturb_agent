#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/07/31
# Udated  on 2026/07/31
# @author: Flavio Lichtenstein
# @local: Home sweet home

"""
Claude - Opus 5 High code

example: PAAD subtyping on the deconvolved compartments
=======================================================

Consumes a `DeconvResult` (theta, Z, states, genes) and answers, in order:

  STEP 0  composition_test()      Is the bulk Program-1/2 split just theta?
                                  If yes, everything downstream is composition.
  STEP 1  celltype_lfc()          Marker LFC per compartment (+ a check that it
                                  is not merely re-reading the reference).
  STEP 1b celltype_de()           Tumor vs normal INSIDE a compartment, via a
                                  theta-interaction model. This is the contrast
                                  most people actually want, and the ONLY one
                                  here that attributes an effect to a cell type.
  STEP 2  cluster()/run_all()     Subtypes within a compartment, k = 2..5,
                                  PCA -> Ward (UMAP for display only), with
                                  diagnostics that say whether the split is
                                  identifiable or an echo of bulk/theta.

Standalone class, not a mixin — nothing gets pasted into PRISM, so no name
collisions. Only hard deps are numpy/pandas/scipy/sklearn; umap and matplotlib
are imported lazily.

--------------------------------------------------------------------------
IDENTIFIABILITY — read before interpreting any cluster
--------------------------------------------------------------------------
The posterior is

    Z[s,k,g] = X[s,g] * w[s,k,g],   w[s,k,g] = theta[s,k] phi[k,g] / sum_k' ...

so the compartment-k profile of sample s is the BULK profile of that sample
reweighted per gene by w[s,k,:], and phi is FIXED across samples. If theta were
constant across samples, w would be constant and every compartment would be one
fixed diagonal rescaling of the bulk — clustering any of them would return the
bulk clusters exactly.

Consequences:
  * Non-malignant compartments carry no sample-specific reference. Their
    between-sample variation comes only from X and from theta. "Fibroblast
    subtypes" recovered this way can be bulk structure in a new coordinate
    system, not fibroblast biology.
  * The malignant compartment IS different, because stage 2 fits a
    sample-specific psi_mal[s,:]. That is the one compartment where the
    clustering has its own degrees of freedom.
  * Therefore every cluster reported here ships with three correlations:
    PC1 vs theta_k, PC1 vs purity, PC1 vs the leading bulk PCs. A compartment
    "subtype" that tracks any of them is not independent evidence.

VERIFIED ON SIMULATION. With a program present ONLY in the malignant
compartment, per-compartment clustering of Z returns ARI = 1.00 against that
program in EVERY compartment, and the cross-compartment ARI matrix is 1.00
everywhere. T cells "have" a malignant program. Projecting out bulk PCs does not
fix this — it deletes the true signal too, because a real malignant program is
itself a leading bulk PC (ARI drops to ~0.00 in all compartments including the
right one). So `regress_bulk` is a DIAGNOSTIC, not a remedy.

WHAT TO DO INSTEAD
  * Malignant subtypes: cluster `psi_mal` (S x G) from stage 2, not Z. That
    matrix has a genuinely sample-specific profile per sample and is the only
    compartment with its own degrees of freedom.  -> cluster_psi()
  * Composition axis: cluster theta directly.     -> composition_test()
  * Gene-level attribution: check whether a gene's compartment assignment is
    supported by the reference at all.            -> gene_reference_support()
  * Always report concordance(): off-diagonal ARI near 1.0 means the
    per-compartment partitions are the same partition and carry no
    compartment-specific information.            -> degeneracy()
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr, ttest_ind, ttest_rel
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             silhouette_score)
from sklearn.preprocessing import StandardScaler


# --------------------------------------------------------------------------- #
# statistics helpers                                                          #
# --------------------------------------------------------------------------- #
def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg. Local implementation: no statsmodels dependency."""
    p = np.asarray(p, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    order = np.argsort(p)
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def _ttest_matrix(A: np.ndarray, B: np.ndarray, paired: bool = False):
    """Vectorised over genes. Degenerate genes -> p = 1."""
    with np.errstate(invalid="ignore", divide="ignore"):
        if paired:
            t, p = ttest_rel(A, B, axis=0, nan_policy="omit")
        else:
            t, p = ttest_ind(A, B, axis=0, equal_var=False, nan_policy="omit")
    t = np.nan_to_num(np.asarray(t, dtype=float), nan=0.0)
    p = np.asarray(p, dtype=float)
    p[~np.isfinite(p)] = 1.0
    return t, p


def _clr(P: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """
    Centred log-ratio. theta lives on the simplex: 
    
    its components are forced to sum to 1, so Euclidean distance between theta rows is not a valid metric and
    PCA on raw theta manufactures spurious negative correlations. 
    
    CLR maps the simplex to a real space where it is.
    """
    P = np.clip(np.asarray(P, dtype=float), eps, None)
    P = P / P.sum(axis=1, keepdims=True)
    L = np.log(P)
    return L - L.mean(axis=1, keepdims=True)


def _regress_out(Y: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Remove the column space of C (samples x m) from Y (samples x genes)."""
    C = np.column_stack([np.ones(len(C)), C])
    beta, *_ = np.linalg.lstsq(C, Y, rcond=None)
    return Y - C @ beta


def _safe_spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 4 or np.std(a[m]) == 0 or np.std(b[m]) == 0:
        return np.nan
    return float(spearmanr(a[m], b[m]).statistic)


# --------------------------------------------------------------------------- #
class CellTypePrograms:
    """
    Parameters
    ----------
    res      : DeconvResult with .theta, .Z, .states, .genes
    df_bulk  : optional genes x samples raw counts. Required for the bulk-echo
               diagnostics and for regress_bulk.
    Z, genes : pass the output of PRISM.full_Z() to work over ALL shared genes
               instead of the marker subset used to fit theta — necessary to
               reach the HOX-antisense lncRNA loci.
    """

    def __init__(self, res, df_bulk: pd.DataFrame | None = None,
                 Z: np.ndarray | None = None, genes: Sequence[str] | None = None):

        self.res = res
        self.theta = res.theta
        self.states = list(res.states)
        self.samples = res.theta.index.astype(str)

        self.Z = np.asarray(Z if Z is not None else res.Z)

        self.genes = pd.Index([str(g) for g in (genes if genes is not None else res.genes)])

        if self.Z is None:
            raise ValueError("No Z. Re-run with return_Z=True, or pass PRISM.full_Z() output.")
        
        if self.Z.shape != (len(self.samples), len(self.states), len(self.genes)):
            raise ValueError(f"Z {self.Z.shape} != (S,K,G) "
                             f"({len(self.samples)},{len(self.states)},{len(self.genes)})")

        self.purity = getattr(res, "tumor_purity", None)
        if isinstance(self.purity, pd.DataFrame):
            self.purity = self.purity.iloc[:, 0]

        self.df_bulk = df_bulk
        self._bulk_pcs = None
        self._cache: dict = {}
        self.results: dict = {}

    # ===================================================================== #
    # matrices                                                              #
    # ===================================================================== #
    def bulk_pcs(self, n: int = 5) -> pd.DataFrame | None:
        """Leading PCs of the bulk log2-CPM matrix, on the same gene set."""
        if self.df_bulk is None:
            return None
        if self._bulk_pcs is None or self._bulk_pcs.shape[1] < n:
            g = self.genes.intersection(self.df_bulk.index)
            B = self.df_bulk.loc[g].T.to_numpy(dtype=float)
            B = np.log2(B / np.clip(B.sum(1, keepdims=True), 1, None) * 1e6 + 1)
            B = B[:, B.std(0) > 0]
            k = int(min(n, B.shape[0] - 1, B.shape[1]))
            Y = PCA(n_components=k, random_state=42).fit_transform(StandardScaler().fit_transform(B))
            self._bulk_pcs = pd.DataFrame(
                Y, index=self.df_bulk.columns.astype(str),
                columns=[f"bulkPC{i+1}" for i in range(k)]).reindex(self.samples)
        return self._bulk_pcs.iloc[:, :n]

    def compartment_matrix(self,
                           state: str,
                           min_theta: float = 0.02,
                           min_frac_expressed: float = 0.25,
                           log: bool = True) -> tuple[pd.DataFrame, dict]:
        """
        samples x genes for one compartment, CPM-normalised WITHIN the
        compartment, then log2(+1).

        The within-compartment normalisation is not cosmetic: sum_g Z[s,k,g] is
        approximately theta[s,k] * library_size[s], so without it PC1 is the
        cell-type fraction and every "subtype" is abundance.

        Samples with theta[s,k] < min_theta are dropped — below that the row is
        dominated by the prior rather than by data.
        """
        key = (state, min_theta, min_frac_expressed, log)
        if key in self._cache:
            return self._cache[key]

        k = self.states.index(state)
        raw = pd.DataFrame(self.Z[:, k, :], index=self.samples, columns=self.genes)

        f = self.theta[state].reindex(self.samples)
        keep = (f.fillna(0).to_numpy() >= min_theta) & (raw.sum(axis=1).to_numpy() > 0)
        raw = raw.loc[keep]

        lib = raw.to_numpy(dtype=float).sum(axis=1, keepdims=True)
        Y = raw.to_numpy(dtype=float) / np.clip(lib, 1e-12, None) * 1e6
        if log:
            Y = np.log2(Y + 1)
        df = pd.DataFrame(Y, index=raw.index, columns=raw.columns)
        df = df.loc[:, (df > 0).mean(axis=0) >= min_frac_expressed]

        info = dict(state=state, n_samples=int(df.shape[0]),
                    n_dropped=int((~keep).sum()), n_genes=int(df.shape[1]),
                    median_theta=float(np.nanmedian(f.loc[df.index])))
        self._cache[key] = (df, info)
        return df, info

    # ===================================================================== #
    # STEP 0 — is the bulk split just composition?                          #
    # ===================================================================== #
    def composition_test(self,
                         ref_labels: pd.Series,
                         k_list: Iterable[int] = (2, 3, 4, 5),
                         use_types: bool = False,
                         verbose: bool = True) -> pd.DataFrame:
        """
        Cluster samples on theta ALONE (CLR-transformed) and score against the
        bulk Program-1/Program-2 labels.

        This is the decisive test for the Program-1 caveat. High ARI means the
        bulk clusters are recoverable from cell-type fractions with no
        expression information at all — i.e. Program-1 is a composition axis
        (parenchyma-depleted, immune-cold) and not an epithelial program. Low
        ARI means the bulk split carries expression structure beyond
        composition, and per-compartment subtyping is worth doing.
        """
        T = self.res.theta_type if use_types else self.theta
        T = T.reindex(self.samples)
        Y = _clr(T.to_numpy(dtype=float))
        n_comp = int(min(Y.shape[1] - 1, Y.shape[0] - 1, 10))
        P = PCA(n_components=n_comp, random_state=42).fit_transform(Y)
        Zl = linkage(P, method="ward")

        ref = ref_labels.copy()
        ref.index = ref.index.astype(str)
        common = self.samples.intersection(ref.dropna().index)
        rvals = ref.loc[common].astype(str).to_numpy()

        rows = []
        for k in k_list:
            lab = pd.Series(fcluster(Zl, t=k, criterion="maxclust"), index=self.samples)
            rows.append(dict(
                k=k,
                silhouette=float(silhouette_score(P, lab.to_numpy())),
                ari_vs_bulk=float(adjusted_rand_score(rvals, lab.loc[common].to_numpy())),
                nmi_vs_bulk=float(normalized_mutual_info_score(rvals, lab.loc[common].to_numpy())),
            ))
        df = pd.DataFrame(rows)

        # how well does purity alone separate the bulk labels?
        if self.purity is not None and len(common) > 5:
            grp = pd.Series(rvals, index=common)
            top = grp.value_counts().index[:2]
            if len(top) == 2:
                a = self.purity.reindex(grp[grp == top[0]].index).dropna()
                b = self.purity.reindex(grp[grp == top[1]].index).dropna()
                if len(a) > 2 and len(b) > 2:
                    auc = (np.greater.outer(a.to_numpy(), b.to_numpy()).mean()
                           + 0.5 * np.equal.outer(a.to_numpy(), b.to_numpy()).mean())
                    df.attrs["purity_auc_bulk_split"] = float(max(auc, 1 - auc))
                    df.attrs["purity_auc_groups"] = tuple(top)

        if verbose:
            print(df.to_string(index=False))
            if "purity_auc_bulk_split" in df.attrs:
                print(f"purity alone separates {df.attrs['purity_auc_groups']} "
                      f"with AUC = {df.attrs['purity_auc_bulk_split']:.3f}")
            best = df.ari_vs_bulk.max()
            print("\n-> " + ("bulk clusters are largely COMPOSITIONAL "
                             "(theta reproduces them); treat Program-1/2 as a TME axis"
                             if best > 0.5 else
                             "bulk clusters are NOT reproduced by theta alone; "
                             "expression structure beyond composition is present"))
        return df

    # ===================================================================== #
    # STEP 1 — compartment identity LFC                                     #
    # ===================================================================== #
    def celltype_lfc(self,
                     paired: bool = True,
                     min_theta: float = 0.02,
                     lfc_cutoff: float = 1.0,
                     fdr_cutoff: float = 0.05,
                     df_gene_annot: pd.DataFrame | None = None,
                     ref: pd.DataFrame | None = None,
                     verbose: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Each compartment vs the mean of the others, across samples.

        paired=True: every sample contributes every compartment, so pairing
        removes between-sample (cohort, purity, library) variance.

        WHAT THIS DOES AND DOES NOT SHOW. Because phi is fixed, this LFC is
        largely a re-reading of the reference through the bulk. Pass `ref` and
        the returned frame gains `lfc_ref` plus a per-compartment correlation:
        r ~ 0.9 means the plumbing is correct and nothing more; a low r means
        the reference and the data disagree, which is a real problem to chase.
        Treat this as QC, not as a finding.
        """
        mats = {s: self.compartment_matrix(s, min_theta=min_theta,
                                           min_frac_expressed=0.0)[0]
                for s in self.states}
        common = mats[self.states[0]].columns
        for s in self.states[1:]:
            common = common.intersection(mats[s].columns)

        ref_log = None
        if ref is not None:
            g = common.intersection(ref.columns)
            common = g
            R = ref.loc[:, g]
            ref_log = np.log2(R.div(R.sum(axis=1), axis=0) * 1e6 + 1)

        rows = []
        for s in self.states:
            A = mats[s].loc[:, common]
            others = [c for c in self.states if c != s]

            if paired:
                idx = A.index
                for c in others:
                    idx = idx.intersection(mats[c].index)
                if len(idx) < 3:
                    continue
                Amat = A.loc[idx].to_numpy()
                Bmat = np.mean([mats[c].loc[idx, common].to_numpy() for c in others], axis=0)
                t, p = _ttest_matrix(Amat, Bmat, paired=True)
                n_in = n_out = len(idx)
            else:
                Amat = A.to_numpy()
                Bmat = pd.concat([mats[c].loc[:, common] for c in others]).to_numpy()
                t, p = _ttest_matrix(Amat, Bmat, paired=False)
                n_in, n_out = Amat.shape[0], Bmat.shape[0]

            d = pd.DataFrame({
                "gene": common, "state": s, "n_in": n_in, "n_out": n_out,
                "mean_in": Amat.mean(0), "mean_out": Bmat.mean(0),
                "lfc": Amat.mean(0) - Bmat.mean(0),
                "tstat": t, "pvalue": p, "fdr": _bh(p),
            })
            if ref_log is not None and s in ref_log.index:
                rest = ref_log.drop(index=s)
                d["lfc_ref"] = (ref_log.loc[s] - rest.mean(axis=0)).reindex(common).to_numpy()
            rows.append(d)

        dfall = pd.concat(rows, ignore_index=True)
        dfall["abs_lfc"] = dfall["lfc"].abs()
        if df_gene_annot is not None:
            key = "symbol" if "symbol" in df_gene_annot.columns else df_gene_annot.columns[0]
            dfall = dfall.merge(df_gene_annot.rename(columns={key: "symbol"}),
                                left_on="gene", right_on="symbol", how="left")

        dfsig = (dfall.query("abs_lfc >= @lfc_cutoff and fdr < @fdr_cutoff")
                      .sort_values(["state", "lfc"], ascending=[True, False])
                      .reset_index(drop=True))

        if verbose:
            print(dfsig.groupby("state").size().rename("n_sig").to_string())
            if "lfc_ref" in dfall.columns:
                r = (dfall.groupby("state")
                          .apply(lambda d: _safe_spearman(d.lfc, d.lfc_ref), include_groups=False)
                          .rename("rho(lfc, lfc_reference)"))
                print("\nrecovery of the reference (QC only, high is expected):")
                print(r.round(3).to_string())
        return dfall, dfsig

    # ===================================================================== #
    # STEP 2 — subtypes inside one compartment                              #
    # ===================================================================== #
    def cluster(self,
                state: str,
                k_list: Iterable[int] = (2, 3, 4, 5),
                method: str = "ward",
                min_theta: float = 0.02,
                n_top_var: int = 2000,
                regress_bulk: int = 0,
                ref_labels: pd.Series | None = None,
                lfc_cutoff: float = 1.0,
                fdr_cutoff: float = 0.05,
                do_umap: bool = False,
                verbose: bool = True) -> dict:
        """
        regress_bulk : project out this many leading bulk PCs. DIAGNOSTIC ONLY.
                       On simulation this removes genuine compartment signal
                       along with the artifact, so a split vanishing under
                       regression is NOT evidence that it was spurious. Read
                       concordance()/degeneracy() instead.
        """
        df_full, info = self.compartment_matrix(state, min_theta=min_theta)
        if df_full.shape[0] < 6:
            if verbose:
                print(f"[{state}] skipped: {df_full.shape[0]} usable samples")
            return dict(info=info, skipped=True)

        Y = df_full.to_numpy(dtype=float)
        if regress_bulk:
            bp = self.bulk_pcs(regress_bulk)
            if bp is None:
                raise ValueError("regress_bulk requires df_bulk.")
            Y = _regress_out(Y, bp.reindex(df_full.index).to_numpy())
        df_res = pd.DataFrame(Y, index=df_full.index, columns=df_full.columns)

        hvg = df_res.loc[:, df_res.var(axis=0).sort_values(ascending=False)
                                 .index[:min(n_top_var, df_res.shape[1])]]
        Xs = StandardScaler().fit_transform(hvg.to_numpy())
        n_comp = int(min(10, Xs.shape[0] - 1, Xs.shape[1]))
        pca = PCA(n_components=n_comp, random_state=42)
        P = pca.fit_transform(Xs)
        df_pca = pd.DataFrame(P[:, :3], index=hvg.index,
                              columns=[f"PC{i+1}" for i in range(min(3, n_comp))])

        # ---- identifiability diagnostics --------------------------------
        pc1 = df_pca["PC1"]
        diag = dict(
            rho_pc1_theta=_safe_spearman(pc1, self.theta[state].reindex(pc1.index)),
            rho_pc1_purity=(_safe_spearman(pc1, self.purity.reindex(pc1.index))
                            if self.purity is not None else np.nan),
            var_pc1=float(pca.explained_variance_ratio_[0]),
        )
        bp = self.bulk_pcs(3)
        if bp is not None:
            rhos = [abs(_safe_spearman(pc1, bp[c].reindex(pc1.index))) for c in bp.columns]
            diag["rho_pc1_bulkPC_max"] = float(np.nanmax(rhos)) if len(rhos) else np.nan

        Zl = linkage(df_pca.to_numpy(), method=method)
        clusters = pd.DataFrame(index=df_pca.index)
        rows, markers = [], {}

        for k in k_list:
            if k >= df_pca.shape[0]:
                continue
            lab = fcluster(Zl, t=k, criterion="maxclust")
            lab_km = KMeans(n_clusters=k, random_state=42,
                            n_init="auto").fit_predict(df_pca.to_numpy()) + 1
            clusters[f"k{k}"] = lab
            clusters[f"k{k}_km"] = lab_km

            row = dict(state=state, k=k, n=int(df_pca.shape[0]),
                       silhouette=float(silhouette_score(df_pca.to_numpy(), lab)),
                       ari_hca_km=float(adjusted_rand_score(lab, lab_km)),
                       sizes=np.bincount(lab)[1:].tolist(), **diag)

            # does the split just separate high- from low-theta samples?
            th = self.theta[state].reindex(df_pca.index).to_numpy()
            row["theta_kruskal_rho"] = _safe_spearman(lab, th)

            if ref_labels is not None:
                r = ref_labels.copy()
                r.index = r.index.astype(str)
                com = df_pca.index.intersection(r.dropna().index)
                if len(com) >= 5:
                    s_lab = pd.Series(lab, index=df_pca.index).loc[com]
                    row["ari_vs_bulk"] = float(adjusted_rand_score(
                        r.loc[com].astype(str).to_numpy(), s_lab.to_numpy()))
            rows.append(row)

            markers[k] = self.cluster_markers(
                df_full, pd.Series(lab, index=df_pca.index),
                lfc_cutoff=lfc_cutoff, fdr_cutoff=fdr_cutoff)

        df_eval = pd.DataFrame(rows)
        df_umap = self._umap(df_pca) if do_umap else None

        if verbose and len(df_eval):
            cols = [c for c in ["k", "n", "silhouette", "ari_hca_km", "ari_vs_bulk",
                                "rho_pc1_theta", "rho_pc1_purity", "rho_pc1_bulkPC_max"]
                    if c in df_eval.columns]
            print(f"\n[{state}] n={info['n_samples']} "
                  f"(dropped {info['n_dropped']} theta<{min_theta}) "
                  f"genes={info['n_genes']} regress_bulk={regress_bulk}")
            print(df_eval[cols].round(3).to_string(index=False))
            self._flag(diag, state)

        out = dict(info=info, diag=diag, df_pca=df_pca, df_umap=df_umap,
                   linkage=Zl, clusters=clusters, eval=df_eval,
                   markers=markers, regress_bulk=regress_bulk, skipped=False)
        self.results[(state, regress_bulk)] = out
        return out

    @staticmethod
    def _flag(diag: dict, state: str):
        msgs = []
        if abs(diag.get("rho_pc1_theta", 0) or 0) > 0.6:
            msgs.append(f"PC1 tracks theta[{state}] (rho={diag['rho_pc1_theta']:.2f}) "
                        "-> abundance, not state")
        if abs(diag.get("rho_pc1_purity", 0) or 0) > 0.6:
            msgs.append(f"PC1 tracks tumor purity (rho={diag['rho_pc1_purity']:.2f})")
        if abs(diag.get("rho_pc1_bulkPC_max", 0) or 0) > 0.7:
            msgs.append(f"PC1 echoes a leading bulk PC "
                        f"(|rho|={diag['rho_pc1_bulkPC_max']:.2f}) "
                        "-> rerun with regress_bulk=2")
        for m in msgs:
            print("   ! " + m)

    def run_all(self, states: Sequence[str] | None = None, **kw):
        states = states or self.states
        evals = []
        for s in states:
            r = self.cluster(s, **kw)
            if not r.get("skipped", False):
                evals.append(r["eval"])
        return (pd.concat(evals, ignore_index=True) if evals else pd.DataFrame())

    # ===================================================================== #
    # markers within a compartment                                          #
    # ===================================================================== #
    def cluster_markers(self, df_mat: pd.DataFrame, labels: pd.Series,
                        lfc_cutoff: float = 1.0, fdr_cutoff: float = 0.05):
        """df_mat: samples x genes (log2 CPM within compartment)."""
        genes = df_mat.columns.to_numpy()
        out = []
        for cid in np.unique(labels):
            s_in = labels.index[labels == cid]
            s_out = labels.index[labels != cid]
            if len(s_in) < 2 or len(s_out) < 2:
                continue
            A = df_mat.loc[s_in].to_numpy()
            B = df_mat.loc[s_out].to_numpy()
            t, p = _ttest_matrix(A, B)
            out.append(pd.DataFrame({
                "gene": genes, "cluster": cid,
                "n_in": len(s_in), "n_out": len(s_out),
                "mean_in": A.mean(0), "mean_out": B.mean(0),
                "lfc": A.mean(0) - B.mean(0),
                "tstat": t, "pvalue": p, "fdr": _bh(p),
            }))
        if not out:
            return pd.DataFrame(), pd.DataFrame()
        dfall = pd.concat(out, ignore_index=True)
        dfall["abs_lfc"] = dfall["lfc"].abs()
        dfsig = (dfall.query("abs_lfc >= @lfc_cutoff and fdr < @fdr_cutoff")
                      .sort_values(["cluster", "lfc", "fdr"], ascending=[True, False, True])
                      .reset_index(drop=True))
        return dfall, dfsig

    # ===================================================================== #
    # is a gene's compartment assignment real, or prior-driven?             #
    # ===================================================================== #
    def gene_reference_support(self, ref: pd.DataFrame,
                               genes: Sequence[str] | None = None,
                               verbose: bool = True) -> pd.DataFrame:
        """
        For each gene: is its compartment allocation driven by the reference, or
        just by theta?

        If phi[:,g] is flat (gene absent or non-specific in the scRNA reference,
        which is the normal situation for HOX-antisense lncRNAs on 10x 3' data),
        then Z[s,k,g] ∝ theta[s,k] and the "compartment share" of that gene is
        the mean theta, carrying no information. Columns:

          ref_counts      total counts for the gene in the reference
          phi_max_share   max_k phi[k,g] / sum_k phi[k,g]   (1/K = flat)
          js_vs_theta     Jensen-Shannon distance between the observed
                          compartment share and mean theta.
                          ~0 -> allocation is the prior. DO NOT INTERPRET.
          top_state       argmax of the observed share
        """
        genes = list(genes) if genes is not None else list(self.genes)
        gi = {g: i for i, g in enumerate(self.genes)}
        genes = [g for g in genes if g in gi]

        R = ref.reindex(columns=self.genes).fillna(0.0)
        R = R.reindex(index=self.states)
        phi = (R + 1e-8).div((R + 1e-8).sum(axis=1), axis=0)
        theta_mean = self.theta.reindex(columns=self.states).mean(axis=0).to_numpy()
        theta_mean = theta_mean / theta_mean.sum()

        rows = []
        for g in genes:
            j = gi[g]
            share = self.Z[:, :, j].sum(axis=0)
            tot = share.sum()
            share = share / tot if tot > 0 else np.full(len(self.states), np.nan)
            pg = phi.iloc[:, j].to_numpy()
            pg = pg / pg.sum() if pg.sum() > 0 else pg
            rows.append(dict(
                gene=g,
                ref_counts=float(R.iloc[:, j].sum()),
                phi_max_share=float(np.nanmax(pg)),
                js_vs_theta=float(jensenshannon(share, theta_mean, base=2))
                if np.isfinite(share).all() else np.nan,
                top_state=self.states[int(np.nanargmax(share))]
                if np.isfinite(share).any() else None,
                top_share=float(np.nanmax(share)),
            ))
        df = pd.DataFrame(rows)
        df["flat_reference"] = df.phi_max_share <= (1.5 / len(self.states))
        df["prior_driven"] = (df.js_vs_theta < 0.15) | df.flat_reference | (df.ref_counts < 10)

        if verbose and len(df):
            n = int(df.prior_driven.sum())
            print(f"{n}/{len(df)} genes have prior-driven compartment allocation "
                  "(no reference support) — their compartment assignment is not evidence.")
        return df.sort_values("js_vs_theta", ascending=False).reset_index(drop=True)

    # ===================================================================== #
    # tumor vs normal, INSIDE a compartment (interaction model)              #
    # ===================================================================== #
    def celltype_de(self,
                    groups: pd.Series,
                    df_bulk: pd.DataFrame | None = None,
                    case: str | None = None,
                    control: str | None = None,
                    covariates: pd.DataFrame | None = None,
                    fdr_cutoff: float = 0.05,
                    min_theta_sd: float = 0.02,
                    df_gene_annot: pd.DataFrame | None = None,
                    verbose: bool = True) -> pd.DataFrame:
        """
        Cell-type-specific differential expression: tumor vs normal WITHIN each
        compartment. This is the TOAST / csSAM / CellDMC framework.

        Bulk CPM is LINEAR in the compartment profiles, so

            y[s,g] = sum_k theta[s,k]*mu[k,g]
                   + sum_k theta[s,k]*case[s]*delta[k,g] + e[s,g]

        and delta[k,g] is the case-vs-control effect inside compartment k. It is
        identified by the VARIATION of theta across samples — real information —
        rather than by the fixed phi, which is why this attributes and the
        per-compartment LFC on Z does not.

        Verified on simulation (DE in fibroblast only): mean |t| 35.7 in
        fibroblast vs 0.5-0.8 elsewhere. The same contrast computed as a
        per-compartment LFC on Z gave 0.850 vs 0.843 — no attribution at all.

        IDENTIFIABILITY. delta[k] needs theta[:,k] to vary across samples AND to
        be non-zero in both groups. Compartments failing that get flagged
        `identifiable=False` — read nothing into their effects. In PAAD the
        malignant compartment ALWAYS fails a tumor-vs-normal contrast, because
        theta_malignant ~ 0 in normal pancreas: there is no malignant
        compartment in normal tissue for the tumor one to be compared against.

        Returns long-format: gene, state, delta, se, tstat, pvalue, fdr.
        """
        df_bulk = df_bulk if df_bulk is not None else self.df_bulk
        if df_bulk is None:
            raise ValueError("celltype_de needs df_bulk (genes x samples raw counts).")

        g = groups.copy()
        g.index = g.index.astype(str)
        samples = self.samples.intersection(g.dropna().index).intersection(
            df_bulk.columns.astype(str))
        g = g.loc[samples].astype(str)

        levels = list(pd.unique(g))
        if case is None or control is None:
            if len(levels) != 2:
                raise ValueError(f"Need exactly 2 groups or explicit case/control; got {levels}")
            control = control or (["normal", "control", "adjacent"] and
                                  next((l for l in levels if l.lower() in
                                        ("normal", "control", "adjacent", "nat")), levels[0]))
            case = case or next(l for l in levels if l != control)
        mask = g.isin([case, control])
        samples, g = samples[mask.to_numpy()], g[mask]

        genes = self.genes.intersection(df_bulk.index)
        B = df_bulk.loc[genes, samples].T.to_numpy(dtype=float)
        Y = B / np.clip(B.sum(1, keepdims=True), 1, None) * 1e6      # linear CPM

        Th = self.theta.reindex(index=samples, columns=self.states).to_numpy(dtype=float)
        T = (g.to_numpy() == case).astype(float)

        # identifiability per compartment
        ident = []
        for k, s in enumerate(self.states):
            sd = float(np.std(Th[:, k]))
            in_case = float(np.mean(Th[T == 1, k]))
            in_ctrl = float(np.mean(Th[T == 0, k]))
            ident.append(dict(state=s, theta_sd=sd, mean_case=in_case, mean_control=in_ctrl,
                              identifiable=bool(sd >= min_theta_sd
                                                and in_case >= 0.01 and in_ctrl >= 0.01)))
        df_ident = pd.DataFrame(ident)

        D = np.hstack([Th, Th * T[:, None]])
        if covariates is not None:
            D = np.hstack([D, covariates.reindex(samples).to_numpy(dtype=float)])
        K = len(self.states)

        cond = float(np.linalg.cond(D))
        beta, *_ = np.linalg.lstsq(D, Y, rcond=None)
        resid = Y - D @ beta
        dof = D.shape[0] - np.linalg.matrix_rank(D)
        if dof < 5:
            raise ValueError(f"Only {dof} residual dof — too few samples for {D.shape[1]} terms.")
        sig2 = (resid ** 2).sum(0) / dof
        XtXi = np.linalg.pinv(D.T @ D)
        delta = beta[K:2 * K]
        se = np.sqrt(np.clip(np.outer(np.diag(XtXi)[K:2 * K], sig2), 1e-30, None))
        tstat = delta / se

        from scipy.stats import t as _tdist
        pval = 2 * _tdist.sf(np.abs(tstat), dof)

        rows = []
        for k, s in enumerate(self.states):
            rows.append(pd.DataFrame({
                "gene": genes, "state": s,
                "delta_cpm": delta[k], "se": se[k],
                "tstat": tstat[k], "pvalue": pval[k], "fdr": _bh(pval[k]),
            }))
        out = pd.concat(rows, ignore_index=True).merge(df_ident, on="state", how="left")
        if df_gene_annot is not None:
            key = "symbol" if "symbol" in df_gene_annot.columns else df_gene_annot.columns[0]
            out = out.merge(df_gene_annot.rename(columns={key: "symbol"}),
                            left_on="gene", right_on="symbol", how="left")
        out.attrs.update(case=case, control=control, n_case=int(T.sum()),
                         n_control=int((1 - T).sum()), dof=int(dof), cond=cond)

        if verbose:
            print(f"contrast: {case} (n={int(T.sum())}) vs {control} (n={int((1-T).sum())}) "
                  f"| {len(genes)} genes | dof={dof} | cond(D)={cond:.1f}")
            if cond > 1e3:
                print("   ! design is ill-conditioned; theta columns are near-collinear")
            print(df_ident.round(3).to_string(index=False))
            sig = (out.query("fdr < @fdr_cutoff and identifiable")
                      .groupby("state").size().rename("n_sig"))
            print("\nsignificant genes per compartment (identifiable only):")
            print(sig.to_string() if len(sig) else "  none")
        return out

    # ===================================================================== #
    # malignant subtypes from the stage-2 sample-specific profile            #
    # ===================================================================== #
    def cluster_psi(self,
                    psi_mal: np.ndarray | pd.DataFrame,
                    genes: Sequence[str] | None = None,
                    k_list: Iterable[int] = (2, 3, 4, 5),
                    method: str = "ward",
                    n_top_var: int = 2000,
                    ref_labels: pd.Series | None = None,
                    lfc_cutoff: float = 1.0,
                    fdr_cutoff: float = 0.05,
                    verbose: bool = True) -> dict:
        """
        Cluster the stage-2 malignant profile psi_mal (S x G, rows sum to 1).

        THIS is the legitimate object for PDAC subtyping. Unlike Z, psi_mal
        carries one profile PER SAMPLE that was estimated rather than imposed,
        so its between-sample structure is not a reweighting of a fixed phi.
        It is what `update_reference()` returns and what `subtype_malignant()`
        should score Moffitt signatures on.
        """
        if isinstance(psi_mal, pd.DataFrame):
            df = psi_mal.copy()
        else:
            g = list(genes) if genes is not None else list(self.genes)
            df = pd.DataFrame(np.asarray(psi_mal), index=self.samples, columns=g)
        df = np.log2(df.div(df.sum(axis=1), axis=0) * 1e6 + 1)

        hvg = df.loc[:, df.var(axis=0).sort_values(ascending=False)
                            .index[:min(n_top_var, df.shape[1])]]
        Xs = StandardScaler().fit_transform(hvg.to_numpy())
        n_comp = int(min(10, Xs.shape[0] - 1, Xs.shape[1]))
        pca = PCA(n_components=n_comp, random_state=42)
        P_ = pca.fit_transform(Xs)
        df_pca = pd.DataFrame(P_[:, :3], index=hvg.index,
                              columns=[f"PC{i+1}" for i in range(min(3, n_comp))])

        diag = dict(var_pc1=float(pca.explained_variance_ratio_[0]),
                    rho_pc1_purity=(_safe_spearman(df_pca["PC1"],
                                                   self.purity.reindex(df_pca.index))
                                    if self.purity is not None else np.nan))
        Zl = linkage(df_pca.to_numpy(), method=method)
        clusters, rows, markers = pd.DataFrame(index=df_pca.index), [], {}

        for k in k_list:
            if k >= df_pca.shape[0]:
                continue
            lab = fcluster(Zl, t=k, criterion="maxclust")
            clusters[f"k{k}"] = lab
            row = dict(source="psi_mal", k=k, n=int(df_pca.shape[0]),
                       silhouette=float(silhouette_score(df_pca.to_numpy(), lab)),
                       sizes=np.bincount(lab)[1:].tolist(), **diag)
            if ref_labels is not None:
                r = ref_labels.copy(); r.index = r.index.astype(str)
                com = df_pca.index.intersection(r.dropna().index)
                if len(com) >= 5:
                    row["ari_vs_bulk"] = float(adjusted_rand_score(
                        r.loc[com].astype(str).to_numpy(),
                        pd.Series(lab, index=df_pca.index).loc[com].to_numpy()))
            rows.append(row)
            markers[k] = self.cluster_markers(df, pd.Series(lab, index=df_pca.index),
                                              lfc_cutoff=lfc_cutoff, fdr_cutoff=fdr_cutoff)

        df_eval = pd.DataFrame(rows)
        if verbose and len(df_eval):
            print("[psi_mal] stage-2 malignant profile")
            print(df_eval.drop(columns=["source"]).round(3).to_string(index=False))
            if abs(diag["rho_pc1_purity"] or 0) > 0.6:
                print(f"   ! PC1 tracks purity (rho={diag['rho_pc1_purity']:.2f}) — "
                      "psi_mal is less well estimated in low-purity samples")
        out = dict(info=dict(state="psi_mal", n_samples=int(df_pca.shape[0])),
                   diag=diag, df_pca=df_pca, df_umap=None, linkage=Zl,
                   clusters=clusters, eval=df_eval, markers=markers, skipped=False)
        self.results[("psi_mal", 0)] = out
        return out

    def degeneracy(self, k: int = 2, regress_bulk: int = 0) -> float:
        """
        Mean off-diagonal ARI of concordance(). Near 1.0 means every compartment
        returned the SAME partition, so per-compartment clustering of Z is
        carrying no compartment-specific information and must not be reported as
        compartment biology.
        """
        M = self.concordance(k=k, regress_bulk=regress_bulk).to_numpy(dtype=float)
        off = M[~np.eye(len(M), dtype=bool)]
        v = float(np.nanmean(off)) if off.size else np.nan
        print(f"mean off-diagonal ARI at k={k}: {v:.3f}  -> "
              + ("DEGENERATE: one partition wearing K hats" if v > 0.8 else
                 "partially compartment-specific" if v > 0.5 else
                 "compartments give distinct partitions"))
        return v

    # ===================================================================== #
    # cross-compartment concordance                                         #
    # ===================================================================== #
    def concordance(self, k: int = 2, regress_bulk: int = 0) -> pd.DataFrame:
        """
        ARI between compartments at fixed k. Under a fixed phi these are
        expected to be HIGH by construction — the informative reading is which
        pairs are unusually low (genuinely private structure) and whether the
        malignant compartment stands apart from the rest.
        """
        keys = [s for s in self.states if (s, regress_bulk) in self.results]
        M = pd.DataFrame(np.nan, index=keys, columns=keys)
        for i, a in enumerate(keys):
            for b in keys[i:]:
                ca = self.results[(a, regress_bulk)]["clusters"].get(f"k{k}")
                cb = self.results[(b, regress_bulk)]["clusters"].get(f"k{k}")
                if ca is None or cb is None:
                    continue
                com = ca.index.intersection(cb.index)
                if len(com) >= 5:
                    v = float(adjusted_rand_score(ca.loc[com], cb.loc[com]))
                    M.loc[a, b] = M.loc[b, a] = v
        return M

    # ===================================================================== #
    # display                                                               #
    # ===================================================================== #
    def _umap(self, df_pca: pd.DataFrame, n_neighbors: int = 15, min_dist: float = 0.2):
        try:
            import umap
        except ImportError:
            return None
        nn = int(max(2, min(n_neighbors, df_pca.shape[0] - 1)))
        Y = umap.UMAP(n_neighbors=nn, min_dist=min_dist, random_state=42
                      ).fit_transform(df_pca.to_numpy())
        return pd.DataFrame(Y, index=df_pca.index, columns=["UMAP1", "UMAP2"])

    def plot(self, result: dict, k: int = 2, figsize=(11, 4)):
        import matplotlib.pyplot as plt
        from scipy.cluster.hierarchy import dendrogram
        lab = result["clusters"][f"k{k}"]
        emb = result["df_umap"] if result.get("df_umap") is not None else result["df_pca"]
        names = list(emb.columns[:2])

        fig, ax = plt.subplots(1, 3, figsize=figsize)
        for c in np.unique(lab):
            m = (lab == c).to_numpy()
            ax[0].scatter(emb.iloc[m, 0], emb.iloc[m, 1], s=45, label=f"c{c}")
        ax[0].set(xlabel=names[0], ylabel=names[1],
                  title=f"{result['info']['state']} k={k}")
        ax[0].legend(fontsize=7)

        th = self.theta[result["info"]["state"]].reindex(emb.index)
        sc = ax[1].scatter(emb.iloc[:, 0], emb.iloc[:, 1], c=th, s=45, cmap="viridis")
        ax[1].set(xlabel=names[0], title="coloured by theta")
        fig.colorbar(sc, ax=ax[1])

        dendrogram(result["linkage"], labels=emb.index.tolist(), leaf_rotation=90, ax=ax[2])
        ax[2].set_title("Ward (PCA space)")
        ax[2].tick_params(labelsize=5)
        plt.tight_layout()
        plt.show()


# --------------------------------------------------------------------------- #
"""
Z_full, genes_full = prism.full_Z(res, df_bulk, ref)
cp = CellTypePrograms(res, df_bulk=df_bulk, Z=Z_full, genes=genes_full)

ref_labels = df_hca.set_index('sample')['cluster']          # bulk Program-1/2

# STEP 0 — the decisive test
cp.composition_test(ref_labels)

# STEP 1 — QC
dfall, dfsig = cp.celltype_lfc(ref=ref, df_gene_annot=prism.df_gene_annot)

# STEP 1b — tumor vs normal inside each compartment (needs BOTH groups in df_bulk)
groups = df_meta['sample_type'].map(lambda v: 'tumor' if 'Tumor' in str(v) else 'normal')
de = cp.celltype_de(groups, case='tumor', control='normal',
                    df_gene_annot=prism.df_gene_annot)
de.query("fdr < 0.05 and identifiable").sort_values('tstat').head(30)

# STEP 2a — check the per-compartment clustering is not degenerate FIRST
ev = cp.run_all(k_list=(2,3,4,5), ref_labels=ref_labels)
cp.degeneracy(k=2)          # if > 0.8, do not report per-compartment subtypes

# STEP 2b — the legitimate malignant subtyping, on the stage-2 profile
psi_mal, psi_env = prism.update_reference(res.Z, mal_idx)
r = cp.cluster_psi(psi_mal, genes=res.genes, ref_labels=ref_labels)

# the lncRNA program: is its compartment assignment even identifiable?
cp.gene_reference_support(ref, genes=['FAM83A-AS1','HOXA10-AS','HOXB-AS3',
                                      'HOXB-AS4','MIR7-3HG','TFAP2A-AS2'])
"""
