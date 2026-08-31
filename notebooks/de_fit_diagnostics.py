"""
Diagnostics for the OLS fit inside factorial_state_de.

Reconstructs the fit outside the method so the intermediates are inspectable.
D has only FOUR distinct rows (the 2x2 cells), so fitted values take exactly
four values per gene -- a per-gene residual-vs-fitted plot is four points and
tells you nothing. The informative views are pooled across genes, and the
p-value histogram most of all.

Usage:
    diag = refit(X_fib2, disc, mc)
    plot_fit_diagnostics(diag, title="fibroblast")
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def refit(X, disc, mc=None, axis_a=None, axis_b=None,
          high_label="high", exclude_program_genes=True):
    """Same fit as factorial_state_de, returning the intermediates."""
    cols = list(disc.columns)
    if len(cols) < 2:
        raise ValueError(f"need >=2 discretised axes, got {cols}")
    axis_a = axis_a or cols[0]
    axis_b = axis_b or cols[1]

    d = disc[[axis_a, axis_b]].dropna()
    idx = X.index.intersection(d.index)
    d, Xs = d.loc[idx], X.loc[idx]

    if exclude_program_genes and mc is not None:
        drop = set(mc.program_genes())
        n0 = Xs.shape[1]
        Xs = Xs.loc[:, [g for g in Xs.columns if g not in drop]]
        print(f"excluded {n0 - Xs.shape[1]} program marker genes; "
              f"{Xs.shape[1]} remain")

    A = (d[axis_a] == high_label).astype(float).values
    B = (d[axis_b] == high_label).astype(float).values
    D = np.column_stack([np.ones(len(A)), A, B, A * B])
    n, p = D.shape

    Y = Xs.values
    beta, *_ = np.linalg.lstsq(D, Y, rcond=None)
    fitted = D @ beta
    resid = Y - fitted
    dof = n - p
    sigma2 = (resid ** 2).sum(axis=0) / dof
    se = np.sqrt(np.outer(np.diag(np.linalg.pinv(D.T @ D)), sigma2))

    names = ["intercept", f"A_{axis_a}", f"B_{axis_b}", "interaction"]
    with np.errstate(invalid="ignore", divide="ignore"):
        tstat = beta / se
    pvals = pd.DataFrame(2 * stats.t.sf(np.abs(tstat), dof).T,
                         index=Xs.columns, columns=names).fillna(1.0)

    return dict(Xs=Xs, D=D, beta=beta, fitted=fitted, resid=resid,
                sigma2=sigma2, se=se, tstat=tstat, pvals=pvals,
                names=names, dof=dof, n=n, p=p, cells=d)


def plot_fit_diagnostics(diag, title="", max_points=40000, seed=0):
    rng = np.random.default_rng(seed)
    resid, fitted = diag["resid"], diag["fitted"]
    flat_r, flat_f = resid.ravel(), fitted.ravel()
    if flat_r.size > max_points:
        sel = rng.choice(flat_r.size, max_points, replace=False)
        flat_r, flat_f = flat_r[sel], flat_f[sel]

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))

    # 1. residual vs fitted. Four vertical bands per gene, smeared into a
    #    cloud once pooled. A fan opening to the right = mean-variance trend.
    ax[0, 0].scatter(flat_f, flat_r, s=1, alpha=.05, rasterized=True)
    ax[0, 0].axhline(0, color="k", lw=.8)
    ax[0, 0].set(xlabel="fitted (log2-CPM)", ylabel="residual",
                 title="residual vs fitted")

    # 2. mean-variance. OLS assumes sigma2 is unrelated to expression level;
    #    RNA-seq usually violates this, which is the argument for a
    #    limma-style moderated variance instead of per-gene sigma2.
    m = diag["Xs"].mean(axis=0).values
    ax[0, 1].scatter(m, np.sqrt(diag["sigma2"]), s=2, alpha=.15, rasterized=True)
    ax[0, 1].set(xlabel="mean log2-CPM", ylabel="residual SD",
                 title="mean-variance trend")

    # 3. QQ of standardised residuals against the normal the t-test assumes.
    sd = resid.std(axis=0, ddof=diag["p"])
    z = (resid / np.where(sd > 0, sd, np.nan)).ravel()
    z = z[np.isfinite(z)]
    if z.size > max_points:
        z = rng.choice(z, max_points, replace=False)
    z.sort()
    q = stats.norm.ppf((np.arange(1, z.size + 1) - .5) / z.size)
    ax[0, 2].plot(q, z, ",", alpha=.4)
    lim = [min(q[0], z[0]), max(q[-1], z[-1])]
    ax[0, 2].plot(lim, lim, "r-", lw=.8)
    ax[0, 2].set(xlabel="normal quantile", ylabel="std residual",
                 title="QQ, residuals")

    # 4-6. p-value histograms per term. THE diagnostic for a DE fit: uniform
    #      under the null, with a spike near 0 when there is signal. A hump in
    #      the middle or a rise at 1 means the model is misspecified -- the
    #      hit count is then not interpretable whatever the FDR says.
    for j, nm in enumerate([n for n in diag["names"] if n != "intercept"]):
        a = ax[1, j]
        pv = diag["pvals"][nm].values
        a.hist(pv, bins=50, range=(0, 1), color="steelblue")
        a.axhline(len(pv) / 50, color="r", ls="--", lw=.9,
                  label="uniform null")
        a.set(xlabel="p", ylabel="genes", title=nm)
        a.legend(fontsize=7)

    fig.suptitle(f"{title}   n={diag['n']}  genes={diag['Xs'].shape[1]}  "
                 f"dof={diag['dof']}", y=1.0)
    plt.tight_layout()
    plt.show()


def cell_means_plot(diag, genes, gene_map=None):
    """Per-gene 2x2 cell means -- the four values the model can produce.

    Use on the handful of genes that reached FDR < 0.05. If a 'hit' is driven
    by one cell with a couple of outlying samples, it shows up here and
    nowhere in the summary table.
    """
    d = diag["cells"]
    key = d.iloc[:, 0].astype(str) + "/" + d.iloc[:, 1].astype(str)
    fig, axes = plt.subplots(1, len(genes), figsize=(3 * len(genes), 3.2),
                             squeeze=False)
    for a, g in zip(axes[0], genes):
        v = diag["Xs"][g]
        lab = (gene_map["symbol"].get(g, g) if gene_map is not None else g)
        for i, k in enumerate(sorted(key.unique())):
            y = v[key == k].values
            a.scatter(np.full(len(y), i) + np.random.uniform(-.1, .1, len(y)),
                      y, s=12, alpha=.6)
            a.plot([i - .2, i + .2], [y.mean()] * 2, "k-", lw=2)
        a.set_xticks(range(len(key.unique())))
        a.set_xticklabels(sorted(key.unique()), rotation=45, ha="right",
                          fontsize=7)
        a.set_title(lab, fontsize=9)
    plt.tight_layout()
    plt.show()


def interaction_plot(diag, genes, gene_map=None, ncol=4, R=None):
    """THE plot for `expr ~ 1 + A + B + A:B`.

    x = factor A (low/high); one line per level of factor B; y = cell mean.
    The four points are everything the model can fit -- with two binary
    factors, OLS is exactly a reparameterisation of the 2x2 cell means:

        cell(0,0) = b0
        cell(1,0) = b0 + bA
        cell(0,1) = b0 + bB
        cell(1,1) = b0 + bA + bB + bAB

    So bA is the horizontal shift of the LOWER line, bB the vertical gap at
    A=0, and bAB the extent to which the two lines are NOT parallel. Parallel
    lines mean the axes act additively; crossing or diverging lines are the
    tumour-stroma crosstalk the interaction term is meant to catch.

    Reading the raw points matters as much as the lines: with ~29 samples per
    cell a "hit" can be two outliers in one corner, and the summary table
    cannot show that.
    """
    d = diag["cells"]
    a_col, b_col = d.columns[0], d.columns[1]
    A = (d[a_col] == "high").values
    B = (d[b_col] == "high").values

    genes = list(genes)
    nrow = int(np.ceil(len(genes) / ncol))
    fig, axes = plt.subplots(nrow, min(ncol, len(genes)),
                             figsize=(3.4 * min(ncol, len(genes)), 3.2 * nrow),
                             squeeze=False)

    short_a = a_col.split(".")[-1].replace("axis_", "")
    short_b = b_col.split(".")[-1].replace("axis_", "")
    lo_lab, hi_lab = short_b.split("_minus_")[::-1]   # iCAF, myCAF

    for ax, g in zip(axes.ravel(), genes):
        y = diag["Xs"][g].values
        lab = gene_map["symbol"].get(g, g) if gene_map is not None else g

        for bval, colour, mark in [(False, "tab:cyan", "o"),
                                   (True, "tab:red", "s")]:

            means = []
            for i, aval in enumerate([False, True]):
                m = (A == aval) & (B == bval)
                pts = y[m]
                ax.scatter(np.full(pts.size, i) + np.random.uniform(-.07, .07, pts.size),
                           pts, s=10, alpha=.45, color=colour)
                means.append(pts.mean() if pts.size else np.nan)

            ax.plot([0, 1], means, marker=mark, color=colour, lw=2,
                   label=f"{hi_lab if bval else lo_lab}-high")

        j = diag["names"]
        bA = diag["beta"][j.index([n for n in j if n.startswith("A_")][0])]
        bB = diag["beta"][j.index([n for n in j if n.startswith("B_")][0])]
        bAB = diag["beta"][j.index("interaction")]
        k = list(diag["Xs"].columns).index(g)

        k = list(diag["Xs"].columns).index(g)
        pv = diag["pvals"].iloc[k]
        nA = [n for n in j if n.startswith("A_")][0]
        nB = [n for n in j if n.startswith("B_")][0]
        if R is not None and g in R.index:
            row = R.loc[g]
            fdrs = {c[4:]: row[c] for c in R.columns if c.startswith("fdr_")}
            sub = "  ".join(f"{k} FDR={v:.3f}" for k, v in fdrs.items() if v < 0.10)  # .split('_')[0]
        else:
            sub = "----"

        m00, m10 = y[(~A) & (~B)].mean(), y[A & (~B)].mean()
        m01, m11 = y[(~A) & B].mean(), y[A & B].mean()
        s_diff = (f"A|Blo {m10-m00:+.2f}  A|Bhi {m11-m01:+.2f}  "
                f"B|Alo {m01-m00:+.2f} (A|Bhi-A|Blo) {(m11-m01)-(m10-m00):+.2f}")

        ax.set_title(f"{lab} ({g})\n{sub}\n{s_diff}", fontsize=7)
        # ax.set_title(f"{lab}\nbA={bA[k]:+.2f} bB={bB[k]:+.2f} bAB={bAB[k]:+.2f}", fontsize=8)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["low", "high"], fontsize=8)
        ax.set_xlabel(short_a, fontsize=8)
        ax.set_ylabel("log2-CPM", fontsize=8)

    for ax in axes.ravel()[len(genes):]:
        ax.axis("off")

    h, l = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False,
            bbox_to_anchor=(0.5, 1.08), fontsize=9,
            title=short_b, title_fontsize=8)

    plt.tight_layout()
    plt.show()
