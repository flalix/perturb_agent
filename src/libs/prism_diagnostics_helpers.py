"""
Diagnostic helpers for prism_malig_lib.py.

Place at MODULE level (outside the class). Everything else in the file is a
method taking `self`; pasting these into the class body would make `A` bind to
`self` on every call.

Requires at module top:
    import numpy as np
    import pandas as pd
    from scipy.stats import fisher_exact, hypergeom
"""


def gene_corr(A, B, min_samples=3):
    """Per-gene Pearson r between two sample x gene matrices, aligned on the
    intersection of both axes. Returns (r, n_samples, n_genes).

    Genes with zero variance in either matrix get NaN, not 0. The naive version
    returns 0.0 for them: 0/0 gives NaN in the standardisation, and pandas
    .sum() skips NaN, so a constant gene reports "uncorrelated" as if it had
    been measured. That is indistinguishable from a real null in the output.

    Genes with any missing values also get NaN rather than a skipna sum over a
    full-n denominator, which would shrink r toward 0 in proportion to how much
    data is missing.
    """
    idx = A.index.intersection(B.index)
    g = A.columns.intersection(B.columns)
    if len(idx) < min_samples:
        raise ValueError(f"only {len(idx)} shared samples (need >= {min_samples})")
    if len(g) == 0:
        raise ValueError("no shared genes -- check gene ID space (symbol vs Ensembl)")

    a, b = A.loc[idx, g], B.loc[idx, g]
    sa, sb = a.std(), b.std()                       # ddof=1, matches / (n-1)
    ok = (sa > 0) & (sb > 0) & a.notna().all() & b.notna().all()

    za = (a - a.mean()).div(sa.where(sa > 0))
    zb = (b - b.mean()).div(sb.where(sb > 0))
    r = (za * zb).sum(axis=0) / (len(idx) - 1)
    return r.where(ok), len(idx), len(g)


def enrich(hits, members, background):
    """One-sided (enrichment) Fisher / hypergeometric for `members` among
    `hits`. All three are intersected with `background` first, so the counts
    are internally consistent whatever was passed in -- the usual failure is a
    `members` set built against a different gene universe than `hits`.

    Returns NaN fold/odds when `hits` or `members` is empty rather than 1.0:
    an undefined enrichment is not a null result.
    """
    background = set(background)
    hits, members = set(hits) & background, set(members) & background
    M, n, N = len(background), len(members), len(hits)
    k = len(hits & members)

    if M == 0 or N == 0 or n == 0:
        return dict(k=k, N=N, n=n, M=M, expected=np.nan, fold=np.nan,
                    odds=np.nan, p=np.nan, p_hyper=np.nan)

    odds, p = fisher_exact([[k, N - k], [n - k, M - n - N + k]],
                           alternative="greater")
    return dict(k=k, N=N, n=n, M=M, expected=N * n / M,
                fold=(k / N) / (n / M), odds=odds, p=p,
                p_hyper=hypergeom.sf(k - 1, M, n, N))


def hits_from(R, contrast, alpha=0.05):
    """Select DE hits by NAMED contrast. Positional selection
    (`R.filter(like="fdr_").iloc[:, 0]`) silently returns a different
    hypothesis when the design gains a level or the columns reorder.
    """
    col = f"fdr_{contrast}"
    if col not in R.columns:
        raise KeyError(f"{col!r} not found; available: "
                       f"{[c for c in R.columns if c.startswith('fdr_')]}")
    return R.index[R[col] < alpha]


def contrasts_of(R):
    """Contrast names from the fdr_ prefix, not from column position."""
    return [c[4:] for c in R.columns if c.startswith("fdr_")]


def boot_median_gap(r, group, n_boot=2000, seed=0):
    """Bootstrap CI on median(in group) - median(out of group).

    For comparing a small, selected subset against the rest -- e.g. genes with
    their own BayesPrism posterior versus genes reconstructed by projection.
    Filters like min_share can leave the subset at n=30-50 while the complement
    is in the thousands, and at that size the point estimate carries a CI of
    roughly +/- 0.03 on a correlation scale. Reporting it bare invites reading
    a difference that is inside sampling noise.

    Returns (gap, lo, hi, n_in, n_out).
    """
    rng = np.random.default_rng(seed)
    m = r.index.isin(set(group))
    a, b = r[m].dropna().values, r[~m].dropna().values
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan, np.nan, len(a), len(b)

    # vectorised: (n_boot, n) index matrices beat a Python loop by ~50x
    ia = rng.integers(0, len(a), size=(n_boot, len(a)))
    ib = rng.integers(0, len(b), size=(n_boot, len(b)))
    d = np.median(a[ia], axis=1) - np.median(b[ib], axis=1)
    lo, hi = np.percentile(d, [2.5, 97.5])
    return float(np.median(a) - np.median(b)), float(lo), float(hi), len(a), len(b)


def pc1(M, return_loadings=False):
    """First principal component of a sample x gene matrix, gene-standardised.

    Returns (scores, variance_explained), or (scores, variance_explained,
    loadings) when return_loadings=True. Loadings are needed to say what the
    component IS; without them a high PC1 concordance between two matrices
    cannot be told from a shared technical axis.

    The sign is arbitrary -- SVD fixes it only up to a joint flip of scores and
    loadings -- so compare components with abs() and read loadings only
    relative to each other.

    Interpretation caution: read the participation ratio, not just the top
    genes. At an effective width of a few hundred genes the top 30 is under 10%
    of the component, and gene-length bias concentrates there, so the visible
    tail can look like a coherent program while the component is something
    else. participation ratio = 1 / sum(loadings**4).
    """
    Z = ((M - M.mean()) / M.std()).fillna(0.0)
    A = Z.values - Z.values.mean(0)
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    scores = pd.Series(U[:, 0] * S[0], index=M.index)
    var = float((S ** 2 / (S ** 2).sum())[0])
    if return_loadings:
        return scores, var, pd.Series(Vt[0], index=M.columns)
    return scores, var


def batch_axis_check(M, batch, gene_lengths=None, top=30):
    """Is this matrix's dominant axis a cohort effect?

    Run before any cross-compartment correlation on pooled data. BayesPrism
    corrects reference-vs-bulk mismatch, not batch WITHIN the bulk, so a
    two-cohort deconvolution puts the same batch axis into every compartment;
    compartments then correlate with each other more strongly than either
    correlates with the bulk, and the excess is cohort membership.

    Reports point-biserial r between PC1 and batch, variance explained,
    participation ratio, and (given gene_lengths) whether the loading tail is
    length-biased -- long transcripts co-vary as a block under differing
    library prep or RNA degradation.
    """
    idx = M.index.intersection(batch.dropna().index)
    scores, var, load = pc1(M.loc[idx], return_loadings=True)

    lv = pd.unique(batch.loc[idx])
    out = {"var_explained": round(var, 3),
           "participation_ratio": int(round(1.0 / (load ** 4).sum())),
           "n_levels": len(lv)}
    if len(lv) == 2:
        y = (batch.loc[idx] == lv[0]).astype(float)
        out["pc1_vs_batch_r"] = round(abs(scores.corr(y)), 3)
    else:
        out["pc1_vs_batch_eta2"] = round(
            float(scores.groupby(batch.loc[idx]).mean().var()
                  / scores.var()), 3)

    if gene_lengths is not None:
        L = np.log(gene_lengths.reindex(load.index).dropna())
        if len(L) > 10:
            out["abs_loading_vs_log_len_rho"] = round(
                float(load.reindex(L.index).abs().corr(L, method="spearman")), 3)
            t = load.abs().nlargest(top).index
            out["top_len_ratio"] = round(
                float(np.exp(L.reindex(t).dropna().median() - L.median())), 2)
    return out
