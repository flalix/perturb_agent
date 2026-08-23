"""
plot_compartment_resid.py
=========================

Visualise the shared tumour-vs-normal trend between two deconvolved
compartments and the residuals against it.

Why the residual is the quantity of interest
--------------------------------------------
The two compartments' LFCs correlate at r ~ 0.955 on genes expressed in both.
A raw per-compartment DEG list is therefore ~91% the shared contrast, and a gene
being significant in the fibroblast arm says almost nothing about fibroblast
biology. What is compartment-specific is the deviation from the fitted line.

Two things the plots are designed to expose
-------------------------------------------
1. The residual is heteroscedastic in AveExpr. Genes at AveExpr 2-4 produce
   large |resid| from small absolute changes (the psi -> 0 boundary), which is
   how PNLIP, AMY2A and HRNR reached the top of a delta ranking despite being
   2000x lower in fibroblast than acinar. Panel B plots resid against expression
   so that class is visible rather than inferred.
2. Panel-level shifts are the evidence, not single genes. Panel C shows the
   distribution of residuals for a marker set against the transcriptome-wide
   null, which is the form the quiescent-fibroblast result takes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def fit_shared_trend(df, x="lfc_aci", y="lfc_fib"):
    """OLS of one compartment's LFC on the other; returns (slope, intercept)."""
    d = df[[x, y]].dropna()
    slope, icpt = np.polyfit(d[x], d[y], 1)
    return float(slope), float(icpt)


def plot_compartment_resid(
    df,
    panels: dict[str, list[str]] | None = None,
    x: str = "lfc_aci",
    y: str = "lfc_fib",
    expr: str = "AveExpr_fib",
    symbol: str = "symbol",
    x_label: str = "acinar",
    y_label: str = "fibroblast",
    expr_floor: float = 2.0,
    n_label: int = 12,
    figsize=(15, 4.6),
):
    """Three-panel diagnostic.

    A  lfc_y vs lfc_x with the fitted shared trend
    B  residual vs expression, with the low-expression zone shaded
    C  residual distribution: all genes vs each named panel
    """
    d = df.dropna(subset=[x, y]).copy()
    slope, icpt = fit_shared_trend(d, x, y)
    d["resid"] = d[y] - (slope * d[x] + icpt)
    r = d[x].corr(d[y])

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    colors = plt.cm.tab10.colors

    # ---- A: the shared trend -------------------------------------------
    ax = axes[0]
    ax.scatter(d[x], d[y], s=4, alpha=0.15, color="0.5", linewidths=0, rasterized=True)
    xx = np.linspace(d[x].min(), d[x].max(), 100)
    ax.plot(xx, slope * xx + icpt, color="crimson", lw=1.5,
            label=f"fit: y = {slope:.2f}x + {icpt:.2f}")
    ax.plot(xx, xx, color="0.3", lw=1, ls=":", label="y = x")
    ax.axhline(0, color="0.8", lw=0.6); ax.axvline(0, color="0.8", lw=0.6)

    if panels:
        for i, (name, genes) in enumerate(panels.items()):
            sub = d[d[symbol].isin(genes)]
            ax.scatter(sub[x], sub[y], s=26, color=colors[i % 10],
                       edgecolors="white", linewidths=0.5, label=f"{name} (n={len(sub)})")
    ax.set_xlabel(f"LFC {x_label}"); ax.set_ylabel(f"LFC {y_label}")
    ax.set_title(f"A. shared trend   r = {r:.3f}   n = {len(d):,}")
    ax.legend(fontsize=7, frameon=False, loc="upper left")

    # ---- B: residual vs expression -------------------------------------
    ax = axes[1]
    ax.scatter(d[expr], d["resid"], s=4, alpha=0.15, color="0.5",
               linewidths=0, rasterized=True)
    ax.axhline(0, color="crimson", lw=1)
    ax.axvspan(d[expr].min(), expr_floor, color="orange", alpha=0.12)
    ax.text(expr_floor, ax.get_ylim()[1] * 0.92,
            f"  AveExpr < {expr_floor:g}\n  boundary artefacts",
            fontsize=7, va="top", color="darkorange")

    # label the largest residuals -- these are the ones a delta ranking returns
    ext = d.reindex(d["resid"].abs().sort_values(ascending=False).index).head(n_label)
    # alternate the offset so dense clusters of extreme genes stay readable
    for k, (_, row) in enumerate(ext.iterrows()):
        dy = 6 if k % 2 == 0 else -9
        ax.annotate(row[symbol], (row[expr], row["resid"]), fontsize=6,
                    xytext=(4, dy), textcoords="offset points",
                    ha="left", va="center")
    ax.scatter(ext[expr], ext["resid"], s=18, facecolors="none",
               edgecolors="crimson", linewidths=0.8)
    if panels:
        for i, (name, genes) in enumerate(panels.items()):
            sub = d[d[symbol].isin(genes)]
            ax.scatter(sub[expr], sub["resid"], s=26, color=colors[i % 10],
                       edgecolors="white", linewidths=0.5)
    ax.set_xlabel(f"AveExpr ({y_label})")
    ax.set_ylabel(f"residual  (LFC {y_label} - fitted)")
    ax.set_title("B. residual vs expression")

    # ---- C: residual distributions -------------------------------------
    ax = axes[2]
    ax.axvline(0, color="0.7", lw=0.8)
    ax.hist(d["resid"], bins=80, color="0.7", density=True,
            label=f"all genes (median {d['resid'].median():+.3f})")
    if panels:
        for i, (name, genes) in enumerate(panels.items()):
            sub = d[d[symbol].isin(genes)]
            if sub.empty:
                continue
            med = sub["resid"].median()
            n_neg = int((sub["resid"] < 0).sum())
            # report whichever direction dominates; "13/13 neg" and "0/4 neg"
            # describe opposite findings and should not share a label format
            dirn = (f"{n_neg}/{len(sub)} neg" if n_neg * 2 >= len(sub)
                    else f"{len(sub) - n_neg}/{len(sub)} pos")
            ax.axvline(med, color=colors[i % 10], lw=2,
                       label=f"{name}: median {med:+.3f}, {dirn}")
            ax.plot(sub["resid"], np.full(len(sub), -0.03 - 0.05 * i), "|",
                    color=colors[i % 10], ms=8, mew=1.2)
    ax.set_xlabel("residual"); ax.set_ylabel("density")
    ax.set_title("C. panel shift vs transcriptome")
    ax.legend(fontsize=7, frameon=False)

    fig.tight_layout()
    return fig, d


if __name__ == "__main__":
    # Panels whose expected direction is known in advance, so they test the
    # pipeline rather than only describing it.
    QUIESCENT = ["C7", "MFAP4", "COL14A1", "THBS4", "ABCA6", "ABCA8", "ABCA9",
                 "ABCA10", "CHRDL1", "SFRP1", "TNXB", "ADAMTSL3", "GPC3"]
    ACINAR_CORE = ["PRSS1", "CPA1", "CPB1", "CELA3A", "CTRB1", "PNLIP", "CLPS",
                   "AMY2A", "CEL", "CTRC", "GP2", "PLA2G1B"]
    BMP_TARGET = ["ID1", "ID2", "ID3", "BAMBI"]

    # df_merge2 must carry: symbol, lfc_fib, lfc_aci, AveExpr_fib
    fig, d = plot_compartment_resid(
        df_merge2,
        panels={"quiescent fibroblast": QUIESCENT,
                "acinar core": ACINAR_CORE,
                "BMP targets": BMP_TARGET},
    )
    fig.savefig("compartment_resid.png", dpi=200, bbox_inches="tight")
    print(d.reindex(d["resid"].abs().sort_values(ascending=False).index)
            .head(15)[["symbol", "lfc_fib", "lfc_aci", "resid", "AveExpr_fib"]]
            .round(3).to_string(index=False))
