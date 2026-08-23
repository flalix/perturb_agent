#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
filter_and_lfc.py
=================

Filter deconvolved samples on df_theta plausibility, then run limma-trend.

Why filtering is needed at all
------------------------------
In the PAAD run the normal arm shows **vertex collapse**: in each normal sample
almost all mass lands on exactly one of {Acinar, Fibroblast, Ductal type 2}, with
the other two at underflow (1e-150 and below, sometimes exact 0). The histograms
have an empty middle -- no normal sample sits between 0.3 and 0.7 acinar -- which
is an optimiser signature, not a tissue composition. Seven normals have no df_theta
at all (all-NaN rows).

So filtering here is **selection on model fit**, not on biology. State that in the
writeup. It is a salvage step; the real fix is upstream (marker coverage for the
collapsing states, and flooring df_theta each iteration so components cannot be
absorbed).

What the filter does and does not buy
-------------------------------------
Keeping samples where the compartment of interest has usable df_theta removes the
worst artefacts. It does not certify the survivors: a sample at acinar 0.986 may
be genuinely acinar-rich or may have been pushed to that vertex. Nothing in df_theta
distinguishes those.

Empirically, in this run, profile quality degrades *gradually* with df_theta rather
than falling off a cliff -- correlation to the reference profile ran 0.760 at
df_theta=0.005 up to 0.851 at df_theta=0.581, monotone. Low-df_theta samples are noisier
(psi = Z/df_theta amplifies assignment noise by 1/df_theta), not prior-collapsed. Noise
inflates within-group variance, which makes the test conservative rather than
anticonservative -- but it is heteroscedastic across samples, so pass
``array_weights=True`` and let limma downweight them.
"""

from __future__ import annotations

import warnings
from pathlib import Path
import sys
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

ROOT_SRC = Path("../src")
sys.path.insert(0, str(ROOT_SRC))

from libs.Basic import pdwritecsv, pdreadcsv, title_replace

__all__ = [
    "theta_report",
    "filter_samples_by_theta",
    "drop_dead_columns",
    "restrict_to_cohort",
    "run_lfc",
]


__version__ = "0.1.0"

# --------------------------------------------------------------------------- #
# diagnostics
# --------------------------------------------------------------------------- #

def theta_report(df_theta: pd.DataFrame, samples: Sequence[str],
                 states: Sequence[str] | None = None) -> pd.DataFrame:
    """Per-sample view of the states that can absorb large mass.

    ``max_component`` near 1 with ``n_underflow`` high is the vertex-collapse
    signature. A well-behaved sample has moderate values across several states.
    """
    states = list(states) if states is not None else list(df_theta.columns)
    sub = df_theta.reindex(list(samples))
    out = pd.DataFrame(index=sub.index)
    out["n_states_gt_001"] = (sub[states] > 0.01).sum(axis=1)
    out["n_underflow"] = (sub[states] < 1e-6).sum(axis=1)
    out["max_component"] = sub[states].max(axis=1)
    out["argmax"] = sub[states].idxmax(axis=1)
    out["all_nan"] = sub[states].isna().all(axis=1)
    return out.join(sub[states])


def drop_dead_columns(df: pd.DataFrame, label: str = "") -> pd.DataFrame:
    """Remove all-NaN and zero-sum columns.

    ``state_expression`` divides Z by df_theta and emits NaN where df_theta collapsed,
    so these appear as all-NaN columns. They must go before any matrix assembly --
    ``build_psi_matrix_and_metadata`` raises on NaN, and a zero column would
    silently distort a fit that tolerated it.
    """
    dead = df.isna().all(axis=0) | (df.fillna(0).sum(axis=0) == 0)
    if dead.any():
        print(f"  {label}dropping {int(dead.sum())}/{df.shape[1]} dead columns")
    return df.loc[:, ~dead]


# --------------------------------------------------------------------------- #
# filtering
# --------------------------------------------------------------------------- #

def filter_samples_by_theta(
    df_theta: pd.DataFrame,
    samples: Sequence[str],
    compartment: str,
    min_theta: float = 0.001,
    require_states: Mapping[str, float] | None = None,
    max_state: Mapping[str, float] | None = None,
    verbose: bool = True,
) -> pd.Index:
    """Keep samples whose df_theta is plausible for ``compartment``.

    Parameters
    ----------
    min_theta
        Floor on the compartment of interest. Below this, psi = Z/df_theta amplifies
        assignment noise so much the profile carries little information.
    require_states
        ``{state: floor}`` -- other states that must also be present. For a normal
        pancreas sample, ``{'Acinar cell': 0.001}`` excludes samples where the
        acinar mass was absorbed elsewhere, which is the collapse signature.
    max_state
        ``{state: ceiling}`` -- e.g. ``{'Ductal cell type 2': 0.30}`` to exclude
        normals modelled as majority malignant.
    """
    sub = df_theta.reindex(list(samples))
    keep = pd.Series(True, index=sub.index)
    reasons: dict[str, int] = {}

    def _apply(mask, why):
        nonlocal keep
        lost = int((keep & ~mask).sum())
        if lost:
            reasons[why] = lost
        keep = keep & mask

    _apply(~sub.isna().all(axis=1), "no df_theta (all-NaN row)")
    if compartment not in sub.columns:
        raise KeyError(f"{compartment!r} not in df_theta columns: {list(sub.columns)}")
    _apply(sub[compartment].fillna(0) > min_theta, f"{compartment} <= {min_theta}")

    for st, floor in (require_states or {}).items():
        _apply(sub[st].fillna(0) > floor, f"{st} <= {floor}")
    for st, ceil in (max_state or {}).items():
        _apply(sub[st].fillna(1) < ceil, f"{st} >= {ceil}")

    kept = pd.Index(sub.index[keep.to_numpy()])
    if verbose:
        print(f"  kept {len(kept)}/{len(sub)}")
        for why, n in reasons.items():
            print(f"    -{n:3d}  {why}")
        if len(kept):
            print(f"    df_theta[{compartment}] range: "
                  f"{sub.loc[kept, compartment].min():.4g} .. "
                  f"{sub.loc[kept, compartment].max():.4g}")
    return kept


def restrict_to_cohort(samples: Sequence[str], cohort: str) -> pd.Index:
    """Keep samples whose ID contains ``cohort`` (e.g. 'TCGA', 'C3L', 'C3N').

    Worth doing. After filtering, the normal arm here is 7 CPTAC + 1 TCGA while
    the tumour arm is mixed. Cohort is then near-confounded with condition: the
    design is not rank-deficient so the R rank check will not stop it, but the
    contrast will absorb cohort differences as disease effects. A within-cohort
    contrast (CPTAC tumours vs CPTAC normals) has no such term.
    """
    return pd.Index([s for s in samples if cohort in s])


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

def run_lfc(
    calc,                                   # CALC_DEGS_PSI instance
    df_psi_tumor: pd.DataFrame,
    df_psi_normal: pd.DataFrame,
    df_theta: pd.DataFrame,
    compartment: str,
    df_gene_map: pd.DataFrame,
    min_theta_tumor: float = 0.01,
    min_theta_normal: float = 0.001,
    require_states: Mapping[str, float] | None = None,
    max_state: Mapping[str, float] | None = None,
    cohort: str | None = None,
    sample_meta: pd.DataFrame | None = None,
    array_weights: bool = True,
    min_prop_detected: float = 0.50,
    root_prism: str | Path = '.',
    force: bool = False,
    verbose: bool = False,
    **limma_kw,
) -> tuple[pd.DataFrame, dict]:
    """Filter both arms on df_theta, then run limma-trend. df_theta
    ``min_prop_detected`` defaults to 0.50 here rather than 0.20. With a normal
    arm in single digits, 0.20 means "detected in 2 samples", which lets genes
    through on almost no evidence and puts them at the psi >= 0 boundary where
    truncation manufactures apparent structure.
    """

    root_prism = Path(root_prism)

    print(f"\n=== {compartment} ===")
    print(" tumour arm:")
    keep_t = filter_samples_by_theta(
        df_theta, df_psi_tumor.columns, compartment=compartment, min_theta=min_theta_tumor)
    print(" normal arm:")
    keep_n = filter_samples_by_theta(
        df_theta, df_psi_normal.columns, compartment=compartment, min_theta=min_theta_normal,
        require_states=require_states, max_state=max_state)

    if cohort:
        keep_t = keep_t.intersection(restrict_to_cohort(keep_t, cohort))
        keep_n = keep_n.intersection(restrict_to_cohort(keep_n, cohort))
        print(f" restricted to {cohort}: {len(keep_t)} tumour / {len(keep_n)} normal")

    T = drop_dead_columns(df_psi_tumor[keep_t], "tumour ")
    N = drop_dead_columns(df_psi_normal[keep_n], "normal ")

    th_n = df_theta.loc[N.columns, compartment]


    compt = compartment.lower()
    fname = f"prism_lfc_compartment_{compt}_min_theta_tumor_{min_theta_tumor}_x_min_theta_normal_{min_theta_normal}_req_stat_{require_states}_min_prop_detected_{min_prop_detected}.tsv"
    fname = title_replace(fname)
    filename =  root_prism / fname

    if filename.exists() and not force:
        df_lfc = pdreadcsv(fname, root_prism, verbose=verbose)

        info = {
            "compartment": compartment,
            "n_tumor": int(T.shape[1]),
            "n_normal": int(N.shape[1]),
            "theta_normal": th_n.describe().to_dict(),
            "kept_normal": list(N.columns),
            "n_sig": int((df_lfc["fdr"] < 0.05).sum()),
        }
            
        return df_lfc, info


    if N.shape[1] < 3:
        raise ValueError(
            f"only {N.shape[1]} normals survive; limma cannot estimate a group "
            "variance. Report median profiles descriptively instead of fitting.")
    if N.shape[1] < 6:
        warnings.warn(
            f"{N.shape[1]} normals: the group mean rests on very few samples and "
            "eBayes moderates variance, not means. Treat p-values as indicative.",
            stacklevel=2)

    spread = float(th_n.max() / max(th_n.min(), 1e-12))
    if spread > 20:
        print(f"  df_theta spans {spread:.0f}x across normals -- "
              f"array_weights strongly recommended" +
              ("" if array_weights else "  [NOT ENABLED]"))

    df_lfc = calc.run_limma(
        T, N,
        mode="limma_trend",
        sample_meta=sample_meta,
        min_prop_detected=min_prop_detected,
        array_weights=array_weights,
        **limma_kw,
    )

    if df_lfc.empty:
        return df_lfc, {"compartment": compartment, 
                        "n_tumor": int(T.shape[1]), 
                        "n_normal": int(N.shape[1]), 
                        "n_sig": 0, "kept_normal": list(N.columns)}

    info = {
        "compartment": compartment,
        "n_tumor": int(T.shape[1]),
        "n_normal": int(N.shape[1]),
        "theta_normal": th_n.describe().to_dict(),
        "kept_normal": list(N.columns),
        "n_sig": int((df_lfc["fdr"] < 0.05).sum()),
    }
    
    print(f"n={info['n_tumor']}T / {info['n_normal']}N fdr<0.05: {info['n_sig']}")

    id2sym = df_gene_map["symbol"].astype(str).to_dict()
    df_lfc.insert(1, "symbol", df_lfc["geneid"].map(id2sym))

    pdwritecsv(df_lfc, fname, root_prism, verbose=verbose)
    
    return df_lfc, info
