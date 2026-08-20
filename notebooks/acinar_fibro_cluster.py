"""
acinar_fibro_cluster.py
=======================

Two compartments (acinar, fibroblast), DEG-filtered against normal, clustered over
tumour samples.

Matrix
------
rows    : (compartment, gene)   -- same gene may appear under both compartments,
                                   because psi_acinar,g and psi_fibro,g are two
                                   different measurements
columns : tumour samples
values  : log1p(psi), z-scored within compartment

Two things are built in deliberately.

**A null for the clustering.** Discrete structure was already rejected on this data
by GMM, consensus clustering and PCA. So a partition here only means something
relative to what unimodal data of the same shape produces. :func:`null_pac` draws
from a Gaussian with the observed covariance spectrum and reruns the identical
consensus procedure; compare against that, not against an absolute PAC cutoff.

**Acinar in and out.** :func:`run` clusters twice. If structure appears only with
the acinar block, the acinar artifacts (RNA integrity, ADM, zero-inflation at the
psi >= 0 boundary) are the first hypothesis, not new biology. The diagnostics that
adjudicate this are computed up front by :func:`cluster_diagnostics`, not after
someone has seen an interesting partition.

Normal reference
----------------
``psi_normal`` must be *deconvolved* normal samples -- compartment-matched. Passing
bulk normals compares a pure compartment against ~85% acinar tissue and makes every
fibroblast gene "up" by composition alone. The bulk path is allowed and warned about
so you can see what it does, but it is not an analysis.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

__all__ = [
    "deg_vs_normal",
    "select_degs",
    "build_matrix",
    "consensus_cluster",
    "null_pac",
    "cluster_diagnostics",
    "run",
]

_RNG = np.random.default_rng


# --------------------------------------------------------------------------- #
# 1. differential expression, per compartment
# --------------------------------------------------------------------------- #

def deg_vs_normal(
    psi_tumor: pd.DataFrame,
    psi_normal: pd.DataFrame,
    *,
    compartment_matched: bool = True,
    min_detect: float = 0.20,
) -> pd.DataFrame:
    """Per-gene tumour-vs-normal contrast within one compartment.

    psi is a posterior expression estimate, not counts, so DESeq2/edgeR negative
    binomial machinery does not apply. Welch's t on log1p(psi) with BH correction,
    plus a rank-based check that does not assume normality.

    ``compartment_matched=False`` records that the denominator is bulk normal and
    warns -- the resulting LFC then confounds abundance with regulation.

    Returns
    -------
    DataFrame indexed by gene: ``lfc, t, p, fdr, p_wilcox, detect_t, detect_n,
    mean_t, mean_n, zero_frac_t``.
    """
    if not compartment_matched:
        warnings.warn(
            "psi_normal is not compartment-matched: LFC will confound cell-type "
            "abundance with disease-associated regulation. Every marker of a "
            "tumour-enriched compartment will appear up.",
            stacklevel=2,
        )

    genes = psi_tumor.index.intersection(psi_normal.index)
    T = np.log1p(psi_tumor.loc[genes].astype(float))
    N = np.log1p(psi_normal.loc[genes].astype(float))

    det_t = (psi_tumor.loc[genes] > 0).mean(axis=1)
    det_n = (psi_normal.loc[genes] > 0).mean(axis=1)
    keep = (det_t >= min_detect) | (det_n >= min_detect)
    T, N = T.loc[keep], N.loc[keep]

    tv, nv = T.to_numpy(), N.to_numpy()
    t_stat, p = stats.ttest_ind(tv, nv, axis=1, equal_var=False)
    p_w = np.array([stats.mannwhitneyu(a, b, alternative="two-sided")[1]
                    if np.ptp(np.r_[a, b]) > 0 else 1.0
                    for a, b in zip(tv, nv)])

    out = pd.DataFrame({
        "lfc": T.mean(axis=1) - N.mean(axis=1),
        "t": t_stat,
        "p": p,
        "p_wilcox": p_w,
        "detect_t": det_t.loc[keep],
        "detect_n": det_n.loc[keep],
        "mean_t": T.mean(axis=1),
        "mean_n": N.mean(axis=1),
        "zero_frac_t": (psi_tumor.loc[keep.index[keep]] == 0).mean(axis=1),
    })
    out["fdr"] = _bh(out["p"].to_numpy())
    out.attrs["compartment_matched"] = compartment_matched
    return out.sort_values("fdr")


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.clip(ranked, 0, 1)
    q[ok] = tmp
    return q


def select_degs(
    deg: pd.DataFrame,
    *,
    fdr: float = 0.05,
    min_abs_lfc: float = 1.0,
    max_zero_frac: float = 0.30,
    top_n: int | None = None,
) -> pd.Index:
    """Filter to usable DEGs.

    ``max_zero_frac`` is the load-bearing one. Genes at the detection floor in a
    large fraction of tumour samples hit the psi >= 0 boundary, and truncation there
    manufactures apparent bimodality -- in simulation, dBIC > 1800 from a
    distribution with no biological structure at all. Those genes will split your
    samples cleanly and the split will mean nothing.
    """
    m = (deg["fdr"] <= fdr) & (deg["lfc"].abs() >= min_abs_lfc) & \
        (deg["zero_frac_t"] <= max_zero_frac)
    sel = deg.loc[m]
    if top_n is not None:
        sel = sel.reindex(sel["lfc"].abs().sort_values(ascending=False).index[:top_n])
    return pd.Index(sel.index, name="gene")


# --------------------------------------------------------------------------- #
# 2. matrix
# --------------------------------------------------------------------------- #

def build_matrix(
    psi: Mapping[str, pd.DataFrame],
    genes: Mapping[str, Sequence[str]],
    *,
    samples: Sequence[str] | None = None,
) -> pd.DataFrame:
    """(compartment, gene) x sample, log1p then z-scored within compartment.

    Per-compartment standardisation is not cosmetic: acinar and fibroblast psi have
    very different variances, and without it the first PC tracks whichever
    compartment varies more rather than any shared structure.
    """
    blocks = []
    for comp, g in genes.items():
        idx = psi[comp].index.intersection(pd.Index(g))
        V = np.log1p(psi[comp].loc[idx].astype(float))
        if samples is not None:
            V = V.reindex(columns=list(samples))
        sd = V.std(axis=1).replace(0, np.nan)
        Z = V.sub(V.mean(axis=1), axis=0).div(sd, axis=0)
        Z.index = pd.MultiIndex.from_product([[comp], idx], names=["compartment", "gene"])
        blocks.append(Z)
    M = pd.concat(blocks).dropna(how="any")
    if M.empty:
        raise ValueError("empty matrix after filtering")
    return M


# --------------------------------------------------------------------------- #
# 3. consensus clustering
# --------------------------------------------------------------------------- #

def _cluster_once(X: np.ndarray, k: int, method: str = "average") -> np.ndarray:
    """X is samples x features. Correlation distance, hierarchical."""
    C = np.corrcoef(X)
    C = np.clip(np.nan_to_num(C, nan=0.0), -1, 1)
    D = squareform(1.0 - C, checks=False)
    return fcluster(linkage(D, method=method), k, criterion="maxclust")


def consensus_cluster(
    M: pd.DataFrame,
    k_range: Sequence[int] = (2, 3, 4, 5, 6),
    *,
    n_resample: int = 200,
    frac_features: float = 0.80,
    frac_samples: float = 0.80,
    seed: int = 0,
) -> dict:
    """Monti-style consensus clustering over tumour samples.

    Returns ``{k: {'consensus': DataFrame, 'pac': float, 'labels': Series}}``.

    PAC = proportion of consensus entries in (0.1, 0.9): the fraction of sample
    pairs the procedure cannot decide about. Lower is crisper. **PAC has no absolute
    threshold** -- compare it against :func:`null_pac`.
    """
    X = M.to_numpy(float).T                     # samples x features
    n, p = X.shape
    rng = _RNG(seed)
    out = {}

    for k in k_range:
        co = np.zeros((n, n))
        cnt = np.zeros((n, n))
        for _ in range(n_resample):
            si = rng.choice(n, max(3, int(frac_samples * n)), replace=False)
            fi = rng.choice(p, max(2, int(frac_features * p)), replace=False)
            lab = _cluster_once(X[np.ix_(si, fi)], k)
            same = (lab[:, None] == lab[None, :]).astype(float)
            co[np.ix_(si, si)] += same
            cnt[np.ix_(si, si)] += 1.0
        with np.errstate(invalid="ignore", divide="ignore"):
            Cm = np.where(cnt > 0, co / cnt, np.nan)
        np.fill_diagonal(Cm, 1.0)

        off = Cm[~np.eye(n, dtype=bool)]
        off = off[np.isfinite(off)]
        pac = float(np.mean((off > 0.1) & (off < 0.9)))

        out[k] = {
            "consensus": pd.DataFrame(Cm, index=M.columns, columns=M.columns),
            "pac": pac,
            "labels": pd.Series(_cluster_once(X, k), index=M.columns, name=f"k{k}"),
        }
    return out


def null_pac(
    M: pd.DataFrame,
    k_range: Sequence[int] = (2, 3, 4, 5, 6),
    *,
    n_null: int = 20,
    seed: int = 0,
    **cc_kw,
) -> pd.DataFrame:
    """PAC under a single-Gaussian null with the observed covariance spectrum.

    SigClust logic: the relevant question is not "does clustering return clusters"
    (it always does) but "is this partition crisper than one drawn from unimodal
    data of the same shape". Null samples are generated in PC space with matched
    eigenvalues, then run through the identical consensus procedure.

    Returns a table with observed PAC, null mean/sd, z, and an empirical p-value.
    """
    X = M.to_numpy(float).T
    n, p = X.shape
    Xc = X - X.mean(0)
    # eigen-spectrum of the sample covariance in sample space
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = (s ** 2) / max(n - 1, 1)

    obs = consensus_cluster(M, k_range, seed=seed, **cc_kw)
    rng = _RNG(seed + 1)

    rows = []
    null = {k: [] for k in k_range}
    for b in range(n_null):
        Zn = rng.normal(size=(n, len(var))) * np.sqrt(var)
        Xn = Zn @ Vt
        Mn = pd.DataFrame(Xn.T, columns=M.columns)
        cn = consensus_cluster(Mn, k_range, seed=seed + 100 + b, **cc_kw)
        for k in k_range:
            null[k].append(cn[k]["pac"])

    for k in k_range:
        nv = np.asarray(null[k], float)
        o = obs[k]["pac"]
        z = (o - nv.mean()) / (nv.std(ddof=1) if nv.std(ddof=1) > 0 else np.nan)
        rows.append({
            "k": k, "pac_obs": o,
            "pac_null_mean": float(nv.mean()), "pac_null_sd": float(nv.std(ddof=1)),
            "z": float(z), "p_emp": float((np.sum(nv <= o) + 1) / (len(nv) + 1)),
        })
    res = pd.DataFrame(rows)
    res.attrs["observed"] = obs
    return res


# --------------------------------------------------------------------------- #
# 4. artifact diagnostics
# --------------------------------------------------------------------------- #

def cluster_diagnostics(
    labels: pd.Series,
    covariates: pd.DataFrame,
) -> pd.DataFrame:
    """Test cluster labels against technical covariates.

    Supply everything you have: acinar zero fraction, RIN / transcript integrity,
    tissue source site, cohort, theta_acinar, theta_fibroblast, library size,
    ADM marker score. Categorical covariates get chi-square, numeric get
    Kruskal-Wallis.

    Compute this **before** looking at the partition. Running it afterwards, on a
    partition that already looks interesting, is how a technical split becomes a
    subtype.
    """
    rows = []
    lab = labels.reindex(covariates.index)
    ok = lab.notna()
    for c in covariates.columns:
        v = covariates.loc[ok, c]
        g = lab[ok]
        if v.dtype == object or str(v.dtype).startswith("category") or v.nunique() <= 6:
            ct = pd.crosstab(g, v)
            stat, p = (stats.chi2_contingency(ct)[:2] if ct.shape[0] > 1 and ct.shape[1] > 1
                       else (np.nan, np.nan))
            test = "chi2"
        else:
            groups = [v[g == u].dropna().to_numpy() for u in g.unique()]
            groups = [x for x in groups if len(x) > 1]
            stat, p = (stats.kruskal(*groups) if len(groups) > 1 else (np.nan, np.nan))
            test = "kruskal"
        rows.append({"covariate": c, "test": test, "stat": stat, "p": p})
    d = pd.DataFrame(rows)
    d["fdr"] = _bh(d["p"].to_numpy())
    return d.sort_values("p")


# --------------------------------------------------------------------------- #
# 5. orchestration
# --------------------------------------------------------------------------- #

def run(
    psi_tumor: Mapping[str, pd.DataFrame],
    psi_normal: Mapping[str, pd.DataFrame],
    *,
    covariates: pd.DataFrame | None = None,
    compartment_matched: bool = True,
    k_range: Sequence[int] = (2, 3, 4, 5, 6),
    deg_kw: dict | None = None,
    sel_kw: dict | None = None,
    cc_kw: dict | None = None,
) -> dict:
    """Full pass: DEGs -> matrix -> clustering with and without acinar.

    The two-way comparison is the point. If structure appears only in the run that
    includes acinar, treat the acinar artifacts as the leading explanation and go to
    ``diagnostics`` before considering biology.
    """
    deg_kw = deg_kw or {}
    sel_kw = sel_kw or {}
    cc_kw = cc_kw or {}

    degs, sel = {}, {}
    for comp in psi_tumor:
        degs[comp] = deg_vs_normal(psi_tumor[comp], psi_normal[comp],
                                   compartment_matched=compartment_matched, **deg_kw)
        sel[comp] = select_degs(degs[comp], **sel_kw)
        print(f"{comp:12s} {len(sel[comp]):5d} DEGs retained "
              f"of {len(degs[comp])} tested")

    samples = psi_tumor["fibroblast"].columns
    results = {}
    for name, comps in [("fibroblast_only", ["fibroblast"]),
                        ("with_acinar", list(psi_tumor))]:
        if not all(c in psi_tumor for c in comps):
            continue
        M = build_matrix({c: psi_tumor[c] for c in comps},
                         {c: sel[c] for c in comps}, samples=samples)
        tab = null_pac(M, k_range, **cc_kw)
        results[name] = {"matrix": M, "null_table": tab,
                         "consensus": tab.attrs["observed"]}
        best = tab.loc[tab["z"].idxmin()] if tab["z"].notna().any() else None
        print(f"\n[{name}] rows={M.shape[0]} samples={M.shape[1]}")
        print(tab.to_string(index=False))
        if best is not None:
            print(f"  most-structured k={int(best['k'])}: "
                  f"PAC {best['pac_obs']:.3f} vs null {best['pac_null_mean']:.3f} "
                  f"(z={best['z']:.2f}, p={best['p_emp']:.3f})")

    if covariates is not None:
        for name, r in results.items():
            k = int(r["null_table"].loc[r["null_table"]["z"].idxmin(), "k"])
            r["diagnostics"] = cluster_diagnostics(
                r["consensus"][k]["labels"], covariates)
            print(f"\n[{name}] diagnostics at k={k}")
            print(r["diagnostics"].to_string(index=False))

    return {"degs": degs, "selected": sel, "results": results}
