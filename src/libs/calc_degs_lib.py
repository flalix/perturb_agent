#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Literal

from numpy.testing import verbose
import pandas as pd

from libs.Basic import pdwritecsv


class CALC_DEGS(object):
    def __init__(self, root_src: Path, run_conda: bool = False):

        self.COMMONS = ["geneid", "symbol", "biotype"]
        self.COMMONS_NORMAL = ["geneid", "symbol"]

        self.root_src = Path(root_src)
        self.libs_dir = root_src / "libs"

        self.run_conda = run_conda

        # R script path
        self.rscript_calc_degs = self.libs_dir / "calc_degs.R"

        if not self.rscript_calc_degs.exists():
            raise FileNotFoundError(f"R script not found: {self.rscript_calc_degs}")

    """
    def running_on_render(self) -> bool:
        return "RENDER" in os.environ or "PORT" in os.environ
    """

    def has_conda(self):
        return shutil.which("conda") is not None

    def deduplicate_by_max_reads(self, df: pd.DataFrame) -> pd.DataFrame:

        if df is None or df.empty:
            return pd.DataFrame()
        
        count_cols = [c for c in df.columns if c not in self.COMMONS]
        if not count_cols:
            raise ValueError("No count columns found.")

        df2 = df.copy()

        for c in count_cols:
            df2[c] = pd.to_numeric(df2[c], errors="coerce").fillna(0)

        df2["_total_reads"] = df2[count_cols].sum(axis=1)

        df2 = (
            df2.sort_values(["geneid", "_total_reads"], ascending=[True, False])
            .drop_duplicates(subset="geneid", keep="first")
            .drop(columns="_total_reads")
            .reset_index(drop=True)
        )

        return df2

    def _find_count_columns(self, df: pd.DataFrame) -> list[str]:
        return [c for c in df.columns if c not in self.COMMONS]

    def _validate_expression_df(self, df: pd.DataFrame, name: str, required_columns: list[str]) -> None:
        missing = [c for c in required_columns if c not in df.columns]
        if missing:
            raise ValueError(f"{name} is missing required columns: {missing}")

        count_cols = self._find_count_columns(df)
        if not count_cols:
            raise ValueError(f"{name} has no count columns besides {self.COMMONS}")

        if df["geneid"].duplicated().any():
            dup_n = int(df["geneid"].duplicated().sum())
            raise ValueError(f"{name} has duplicated geneid values ({dup_n} duplicates).")

    def _rename_count_columns(self, df: pd.DataFrame, prefix: str) -> pd.DataFrame:
        count_cols = self._find_count_columns(df)
        rename_map = {old: f"{prefix}_{i + 1}" for i, old in enumerate(count_cols)}
        return df.rename(columns=rename_map).copy()

    def build_counts_and_metadata(
        self,
        df_tumor: pd.DataFrame,
        df_normal: pd.DataFrame,
        how: Literal["inner", "outer"] = "inner",
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Build:
        1) counts matrix: gene_id, symbol, gene_type, normal_1..., tumor_1...
        2) metadata: sample, condition

        Parameters
        ----------
        how : "inner" or "outer"
            Merge strategy on genes. Usually "inner" is safer.
        """
        self._validate_expression_df(df_tumor, "df_tumor", self.COMMONS)
        self._validate_expression_df(df_normal, "df_normal", self.COMMONS_NORMAL)

        df_turmor_rev = self._rename_count_columns(df_tumor, "tumor")
        df_normal_rev = self._rename_count_columns(df_normal, "normal")

        if 'biotype' in df_turmor_rev.columns:
            df_turmor_rev = df_turmor_rev.drop(columns=['biotype'])
 
        if 'symbol' in df_turmor_rev.columns:
            df_turmor_rev = df_turmor_rev.drop(columns=['symbol'])

        if 'biotype' in df_normal_rev.columns:
            df_normal_rev = df_normal_rev.drop(columns=['biotype'])

        if 'symbol' in df_normal_rev.columns:
            df_normal_rev = df_normal_rev.drop(columns=['symbol'])

        dfn = pd.merge(
            df_turmor_rev, df_normal_rev, on='geneid', how=how, validate="one_to_one"
        )

        count_cols = [c for c in dfn.columns if c not in self.COMMONS]
        if not count_cols:
            raise ValueError("No count columns after merging tumor and normal dataframes.")

        for c in count_cols:
            dfn[c] = pd.to_numeric(dfn[c], errors="coerce").fillna(0)

        # round and convert to int for RNA-seq counts
        dfn[count_cols] = dfn[count_cols].round().astype(int)

        sample_cols = [c for c in dfn.columns if c.startswith("normal_")] + [
            c for c in dfn.columns if c.startswith("tumor_")
        ]

        meta = pd.DataFrame(
            {
                "sample": sample_cols,
                "condition": ["normal"] * sum(c.startswith("normal_") for c in sample_cols)
                + ["tumor"] * sum(c.startswith("tumor_") for c in sample_cols),
            }
        )

        counts = dfn[ ['geneid'] + sample_cols].copy()
        return counts, meta

    def run_deg_rscript(
        self,
        df_tumor: pd.DataFrame,
        df_normal: pd.DataFrame,
        method: str = "auto",  # Literal["auto", "deseq2", "edger"]
        manual_dispersion: float = 0.1,
        min_total_count: int = 10,
        merge_how: Literal["inner", "outer"] = "inner",
        keep_temp: bool = False,
        conda_env: str = "renv",
    ) -> pd.DataFrame:
        
        """
        Run DEG analysis in R using DESeq2 or edgeR through a conda environment.

        Parameters
        ----------
        df_tumor, df_normal
            Expression tables with columns:
            gene_id, symbol, gene_type, counts1, counts2, ...
        method
            "auto", "deseq2", or "edger"
        manual_dispersion
            Used by edgeR fallback when replication is limited
        min_total_count
            Minimum total counts across all samples to keep a gene
        merge_how
            How to merge tumor and normal tables: "inner" or "outer"
        keep_temp
            Whether to keep temporary files
        conda_env
            Conda environment name containing Rscript and required R packages

        Returns
        -------
        pd.DataFrame
            DEG results table
        """

        df_counts, df_meta = self.build_counts_and_metadata(
            df_tumor=df_tumor,
            df_normal=df_normal,
            how=merge_how
        )

        # _ = pdwritecsv(df_counts, "counts.tsv")
        # _ = pdwritecsv(df_meta, "meta.tsv")

        tmpdir_obj = tempfile.TemporaryDirectory()
        tmpdir = Path(tmpdir_obj.name)

        try:
            counts_file = tmpdir / "counts.tsv"
            meta_file = tmpdir / "meta.tsv"
            out_file = tmpdir / "lfc_results.tsv"

            df_counts.to_csv(counts_file, sep="\t", index=False)
            df_meta.to_csv(meta_file, sep="\t", index=False)

            cmd = [
                "Rscript",
                str(self.rscript_calc_degs),
                "--counts",
                str(counts_file),
                "--meta",
                str(meta_file),
                "--out",
                str(out_file),
                "--method",
                method,
                "--manual-dispersion",
                str(manual_dispersion),
                "--min-total-count",
                str(min_total_count),
            ]

            if self.run_conda and self.has_conda():
                print(f">>> LFC calc with {method} - running inside conda environment.")
                cmd = [
                    "conda",
                    "run",
                    "-n",
                    conda_env,
                ] + cmd

            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )

            if proc.returncode != 0:
                raise RuntimeError(
                    f"R DEG script failed -> {str(self.rscript_calc_degs)}.\n"
                    f"Command: {' '.join(cmd)}\n\n"
                    f"STDOUT:\n{proc.stdout.strip()}\n\n"
                    f"STDERR:\n{proc.stderr.strip()}"
                )
            else:
                if verbose:
                    print(f"R DEG script OK -> {str(self.rscript_calc_degs)}.\n")
                    print(f"Command: {' '.join(cmd)}\n\n")

            if not out_file.exists():
                raise RuntimeError(
                    "R script finished without error, but the output file was not created: "
                    f"{out_file}"
                )
            else:
                if verbose:
                    print(f">>> File saved at: {out_file}")

            df_lfc = pd.read_csv(out_file, sep="\t")

            if keep_temp:
                print(f"Temporary DEG files kept at: {tmpdir}")
                tmpdir_obj = None  # prevent cleanup below

            return df_lfc

        finally:
            if (tmpdir_obj is not None) and (not keep_temp):
                tmpdir_obj.cleanup()
