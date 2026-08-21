#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
calc_degs_psi_lib.py
====================

Subclass of :class:`CALC_DEGS` for BayesPrism psi input.

Two changes from the parent, both blocking for this use case:

1. **No integer rounding.** ``build_counts_and_metadata`` ends with
   ``dfn[count_cols].round().astype(int)``. psi values are small posterior
   expression estimates; rounding collapses most of the range to 0 and 1. Correct
   for counts, fatal for psi.

2. **Sample IDs preserved.** The parent renames columns to ``tumor_1..N`` /
   ``normal_1..N``. Those IDs are needed downstream to join RIN, tissue source
   site, cohort and theta for cluster diagnostics, and to identify matched
   tumour/normal pairs. Here the original ID is carried in ``orig_sample`` and the
   df_matrix column is prefixed only to avoid collisions.

Also adds ``--pair-col`` plumbing, so CPTAC3 matched adjacent normals can be run
as a within-patient contrast.

Method
------
``limma_trend`` is the default and the only one that should be used on psi.
DESeq2 and edgeR model counts with a negative binomial; psi is a posterior mean,
not an observed count, and feeding a rounded estimate to a NB treats an estimate
as data. limma-trend works on log-expression with an intensity-dependent variance
trend, which matches the shape of psi.

Caveat for every method here: none propagates deconvolution uncertainty. psi
enters as if measured exactly, so p-values are anticonservative. Weight the
effect-size and detection filters accordingly.
"""

from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Literal, Mapping, Sequence

import numpy as np
import pandas as pd

from libs.calc_degs_lib import CALC_DEGS


class CALC_DEGS_PSI(CALC_DEGS):

    def __init__(self, root_src: Path, run_conda: bool = True,
                 rscript_name: str = "calc_degs_v2.R"):
        super().__init__(root_src, run_conda=run_conda)
        self.rscript_calc_degs = self.libs_dir / rscript_name
        if not self.rscript_calc_degs.exists():
            raise FileNotFoundError(f"R script not found: {self.rscript_calc_degs}")

    # ------------------------------------------------------------------ #
    # df_matrix assembly
    # ------------------------------------------------------------------ #

    def build_psi_matrix_and_metadata(
        self,
        psi_tumor: pd.DataFrame,
        psi_normal: pd.DataFrame,
        pairs: Mapping[str, str] | None = None,
        sample_meta: pd.DataFrame | None = None,
        how: Literal["inner", "outer"] = "inner",
        gene_index_name: str = "geneid",
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Assemble the df_matrix and metadata for one compartment.

        Parameters
        ----------
        psi_tumor, psi_normal
            genes x samples, linear scale, **not** rounded. Both must be the same
            compartment: psi_fibroblast tumour against psi_fibroblast normal.
            Passing bulk normals here compares a pure compartment against a tissue
            that is ~85% acinar, and every marker comes back "up" by composition.
        pairs
            ``{tumor_sample_id: normal_sample_id}`` for matched patients. Produces
            a ``pair`` column; pass ``pair_col='pair'`` to :meth:`run_limma`.
            Unmatched samples get their own singleton level.
        sample_meta
            Optional per-sample covariates indexed by original sample ID; merged
            into the metadata so they can be used with ``covariates=``.

        Returns
        -------
        (df_matrix, df_meta)
            ``df_matrix`` has ``geneid`` plus one column per sample.
            ``df_meta`` has ``sample, orig_sample, condition`` and any extras.
        """
        if psi_tumor.empty or psi_normal.empty:
            print("psi_tumor and psi_normal must both be non-empty")
            return pd.DataFrame(), pd.DataFrame()

        t = psi_tumor.copy()
        n = psi_normal.copy()
        t.index = t.index.astype(str)
        n.index = n.index.astype(str)

        genes = (t.index.intersection(n.index) if how == "inner"
                 else t.index.union(n.index))
        genes = pd.Index(sorted(set(genes)), name=gene_index_name)
        if len(genes) == 0:
            print("no genes shared between tumour and normal psi")
            return pd.DataFrame(), pd.DataFrame()

        t = t.reindex(index=genes).astype(float)
        n = n.reindex(index=genes).astype(float)
        if how == "outer":
            t = t.fillna(0.0)
            n = n.fillna(0.0)

        # prefix only to guarantee uniqueness; the real ID is kept in orig_sample
        t_cols = [f"T__{s}" for s in t.columns]
        n_cols = [f"N__{s}" for s in n.columns]
        t.columns, n.columns = t_cols, n_cols

        mat = pd.concat([n, t], axis=1)
        if mat.isna().any().any():
            print("NaN in assembled df_matrix; check gene indices")
            return pd.DataFrame(), pd.DataFrame()

        if (mat < 0).any().any():
            print("negative values in psi; expected non-negative expression")
            return pd.DataFrame(), pd.DataFrame()

        df_matrix = mat.reset_index().rename(columns={"index": gene_index_name})
        if gene_index_name not in df_matrix.columns:
            df_matrix.insert(0, gene_index_name, genes)

        df_meta = pd.DataFrame({
            "sample": n_cols + t_cols,
            "orig_sample": list(psi_normal.columns) + list(psi_tumor.columns),
            "condition": ["normal"] * n.shape[1] + ["tumor"] * t.shape[1],
        })

        if pairs:
            inv = {v: k for k, v in pairs.items()}
            lab = []
            for s, cond in zip(df_meta["orig_sample"], df_meta["condition"]):
                if cond == "tumor" and s in pairs:
                    lab.append(f"P_{s}")
                elif cond == "normal" and s in inv:
                    lab.append(f"P_{inv[s]}")
                else:
                    lab.append(f"U_{cond}_{s}")
            df_meta["pair"] = lab
            n_pair = sum(1 for x in lab if x.startswith("P_")) // 2
            print(f"  paired patients: {n_pair} "
                  f"({df_meta.shape[0] - 2*n_pair} unmatched samples)")

        if sample_meta is not None:
            extra = sample_meta.reindex(df_meta["orig_sample"]).reset_index(drop=True)
            for c in extra.columns:
                if c not in df_meta.columns:
                    df_meta[c] = extra[c].to_numpy()

        self._sanity_report(mat, df_meta)
        return df_matrix, df_meta

    @staticmethod
    def _sanity_report(mat: pd.DataFrame, df_meta: pd.DataFrame) -> None:
        """Print the numbers that decide whether the contrast is estimable."""
        for cond in ("normal", "tumor"):
            cols = df_meta.loc[df_meta.condition == cond, "sample"]
            sub = mat[cols]
            zf = float((sub == 0).mean().mean())
            print(f"  {cond:6s} n={len(cols):3d}  "
                  f"median={np.median(sub.to_numpy()):.4g}  "
                  f"zero_frac={zf:.3f}")
        if df_meta.groupby("condition").size().min() < 3:
            warnings.warn("fewer than 3 samples in one condition; "
                          "moderated t will be unstable", stacklevel=3)

    # ------------------------------------------------------------------ #
    # R invocation
    # ------------------------------------------------------------------ #

    def run_limma(
        self,
        psi_tumor: pd.DataFrame,
        psi_normal: pd.DataFrame,
        mode: Literal["limma_trend", "limma_voom"] = "limma_trend",
        pairs: Mapping[str, str] | None = None,
        sample_meta: pd.DataFrame | None = None,
        covariates: Sequence[str] | None = None,
        pair_col: str | None = None,
        already_log: bool = False,
        normalize: Literal["none", "quantile", "scale"] = "none",
        array_weights: bool = False,
        min_prop_detected: float = 0.20,
        prior_count: float = 0.5,
        merge_how: Literal["inner", "outer"] = "inner",
        conda_env: str = "renv",
        keep_temp: bool = False,
        verbose: bool = True,
    ) -> pd.DataFrame:
        """Run limma for one compartment and return the results table.

        ``covariates`` are df_meta column names added to the design. A covariate
        perfectly confounded with condition -- tumours from one cohort and normals
        from another -- makes the design rank-deficient and the R script will
        refuse it rather than fit something meaningless.
        """
        df_matrix, df_meta = self.build_psi_matrix_and_metadata(
            psi_tumor, psi_normal, pairs=pairs,
            sample_meta=sample_meta, how=merge_how)

        if df_matrix.empty or df_meta.empty:
            print("Matrix assembly failed; See message above.")
            return pd.DataFrame()

        if covariates:
            miss = [c for c in covariates if c not in df_meta.columns]
            if miss:
                print(f"covariates not in metadata: {miss};  pass them via sample_meta")
                return pd.DataFrame()
            
        if pair_col and pair_col not in df_meta.columns:
            print(f"pair_col '{pair_col}' not in metadata;  pass pairs= to create it")
            return pd.DataFrame()

        tmpdir_obj = tempfile.TemporaryDirectory()
        tmpdir = Path(tmpdir_obj.name)
        try:
            fname_mat = tmpdir / "df_matrix.tsv"
            fname_meta = tmpdir / "df_meta.tsv"
            fname_out = tmpdir / "limma_results.tsv"
            
            df_matrix.to_csv(fname_mat, sep="\t", index=False)
            df_meta.to_csv(fname_meta, sep="\t", index=False)

            cmd = [
                "Rscript", str(self.rscript_calc_degs),
                "--counts", str(fname_mat),
                "--meta", str(fname_meta),
                "--out", str(fname_out),
                "--method", mode,
                "--min-prop-detected", str(min_prop_detected),
                "--prior-count", str(prior_count),
                "--normalize", normalize,
            ]
            if already_log:
                cmd.append("--already-log")
            if array_weights:
                cmd.append("--array-weights")
            if covariates:
                cmd += ["--covariates", ",".join(covariates)]
            if pair_col:
                cmd += ["--pair-col", pair_col]
            if self.run_conda and self.has_conda():
                cmd = ["conda", "run", "-n", conda_env] + cmd

            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if proc.returncode != 0:
                print(f"limma run failed.\nCommand: {' '.join(cmd)}\n\n STDOUT:\n{proc.stdout.strip()}\n\nSTDERR:\n{proc.stderr.strip()}")
                return pd.DataFrame()

            if verbose and proc.stdout.strip():
                print(proc.stdout.strip())
            if not fname_out.exists():
                print(f"R finished without error but wrote no output: {fname_out}")
                return pd.DataFrame()

            df = pd.read_csv(fname_out, sep="\t")

            try:
                df = df.rename(columns={"padj": "fdr", "log2FoldChange": "lfc"})
                df["abs_lfc"] = df["lfc"].abs()
            except KeyError as e:
                print(f"unexpected R output schema: {e}; columns={list(df.columns)}")
                return pd.DataFrame()
            
            if keep_temp:
                print(f"temp kept at: {tmpdir}")
                tmpdir_obj = None
            return df
        finally:
            if tmpdir_obj is not None and not keep_temp:
                tmpdir_obj.cleanup()

    # ------------------------------------------------------------------ #
    # multi-compartment driver
    # ------------------------------------------------------------------ #

    def run_all_compartments(
        self,
        psi_tumor: Mapping[str, pd.DataFrame],
        psi_normal: Mapping[str, pd.DataFrame],
        **kw,
    ) -> dict[str, pd.DataFrame]:
        """Run limma once per compartment. Returns ``{compartment: results}``.

        Each compartment is a separate model; there is no pooling across them, and
        FDR is controlled within compartment. That is the right unit -- a gene in
        acinar and the same gene in fibroblast are two different measurements.
        """
        out = {}
        for comp in psi_tumor:
            if comp not in psi_normal:
                warnings.warn(f"no normal psi for compartment {comp!r}; skipped", stacklevel=2)
                continue

            print(f"\n=== {comp} ===")
            df = self.run_limma(psi_tumor[comp], psi_normal[comp], **kw)
            if df.empty:
                out[comp] = df
                continue
        
            df["compartment"] = comp

            n_degs = df[ (df["fdr"] < 0.05) & (df["abs_lfc"] >= 1)].shape[0]
            n_up   = df[ (df["fdr"] < 0.05) & (df["lfc"] >= 1) ].shape[0]
            n_dw   = df[ (df["fdr"] < 0.05) & (df["lfc"] <= -1) ].shape[0]
            print(f" cutoffs: fdr=0.05 lfc=1: {n_degs} DEGs (up {n_up} / down {n_dw})")
            out[comp] = df
        return out


# --------------------------------------------------------------------------- #
# annotation against the marker panels
# --------------------------------------------------------------------------- #

def annotate_panels(
    deg: pd.DataFrame,
    panels: Mapping[str, Sequence[str]],
    *,
    id_col: str = "geneid",
) -> pd.DataFrame:
    """Label DEG rows by marker-panel membership.

    The DEG scan stays transcriptome-wide; this only flags which results fall in
    a curated panel so those can be inspected first. Panel genes are worth looking
    at because their direction is *predictable*, which makes them a test of the
    pipeline rather than a discovery:

    * acinar core in the acinar compartment -> should be near zero after proper
      deconvolution. Strongly down means residual abundance confounding.
    * quiescent_stellate in fibroblast -> should be down. If not, psi_fibroblast
      is not capturing activation.
    * acinar core in the *fibroblast* compartment -> should show nothing at all.
      Any signal there is leakage, with unambiguous ground truth.
    """
    d = deg.copy()
    lookup: dict[str, list[str]] = {}
    for name, genes in panels.items():
        for g in genes:
            lookup.setdefault(str(g), []).append(name)
    d["panel"] = d[id_col].astype(str).map(lambda g: ",".join(lookup.get(g, [])) or None)
    d["in_panel"] = d["panel"].notna()
    return d


def panel_summary(deg: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Per-panel direction summary: n, n significant, median LFC, fraction up."""
    d = deg.loc[deg["in_panel"]].copy()
    if d.empty:
        return pd.DataFrame()
    d = d.assign(panel=d["panel"].str.split(",")).explode("panel")
    g = d.groupby("panel")
    return pd.DataFrame({
        "n": g.size(),
        "n_degs": g.apply(lambda x: int((x["fdr"] < alpha).sum()), include_groups=False),
        "median_lfc": g["lfc"].median(),
        "frac_up": g.apply(lambda x: float((x["lfc"] > 0).mean()), include_groups=False),
    }).sort_values("median_lfc")
