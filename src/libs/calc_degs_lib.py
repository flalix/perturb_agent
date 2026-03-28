#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/19
# Udated  on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Literal

import pandas as pd


class CALC_DEGS(object):
    def __init__(self, root_data:str='../data/'):

        self.GENE_COLS = ["gene_id", "symbol", "gene_type"]
            
        self.root_data = Path(root_data).resolve()

        # resolve project structure robustly
        self.src_dir = Path(__file__).resolve().parent.parent
        self.libs_dir = self.src_dir / "libs"

        # R script path
        self.rscript_path = self.libs_dir / "calc_degs.R"

        if not self.rscript_path.exists():
            raise FileNotFoundError(f"R script not found: {self.rscript_path}")



    def _find_count_columns(self, df: pd.DataFrame) -> list[str]:
        return [c for c in df.columns if c not in self.GENE_COLS]


    def _validate_expression_df(self, df: pd.DataFrame, name: str) -> None:
        missing = [c for c in self.GENE_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"{name} is missing required columns: {missing}")

        count_cols = self._find_count_columns(df)
        if not count_cols:
            raise ValueError(f"{name} has no count columns besides {self.GENE_COLS}")

        if df["gene_id"].duplicated().any():
            dup_n = int(df["gene_id"].duplicated().sum())
            raise ValueError(f"{name} has duplicated gene_id values ({dup_n} duplicates).")


    def _rename_count_columns(self, df: pd.DataFrame, prefix: str) -> pd.DataFrame:
        count_cols = self._find_count_columns(df)
        rename_map = {old: f"{prefix}_{i+1}" for i, old in enumerate(count_cols)}
        return df.rename(columns=rename_map).copy()


    def build_counts_and_metadata(self, df_tumor: pd.DataFrame, df_normal: pd.DataFrame, 
                                 how: Literal["inner", "outer"] = "inner") -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Build:
        1) counts matrix: gene_id, symbol, gene_type, normal_1..., tumor_1...
        2) metadata: sample, condition

        Parameters
        ----------
        how : "inner" or "outer"
            Merge strategy on genes. Usually "inner" is safer.
        """
        self._validate_expression_df(df_tumor, "df_tumor")
        self._validate_expression_df(df_normal, "df_normal")

        df_turmor_rev = self._rename_count_columns(df_tumor, "tumor")
        df_normal_rev = self._rename_count_columns(df_normal, "normal")

        dfn = pd.merge(df_normal_rev, df_turmor_rev, on=self.GENE_COLS, how=how, validate="one_to_one")

        count_cols = [c for c in dfn.columns if c not in self.GENE_COLS]
        if not count_cols:
            raise ValueError("No count columns after merging tumor and normal dataframes.")

        for c in count_cols:
            dfn[c] = pd.to_numeric(dfn[c], errors="coerce").fillna(0)

        # round and convert to int for RNA-seq counts
        dfn[count_cols] = dfn[count_cols].round().astype(int)

        sample_cols = [c for c in dfn.columns if c.startswith("normal_")] + [
            c for c in dfn.columns if c.startswith("tumor_")
        ]

        meta = pd.DataFrame({
            "sample": sample_cols,
            "condition": ["normal"] * sum(c.startswith("normal_") for c in sample_cols)
                        + ["tumor"] * sum(c.startswith("tumor_") for c in sample_cols),
        })

        counts = dfn[self.GENE_COLS + sample_cols].copy()
        return counts, meta


    def run_deg_rscript(
        self, df_tumor: pd.DataFrame, df_normal: pd.DataFrame,
        method: Literal["auto", "deseq2", "edger"] = "auto",
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


        counts_df, meta_df = self.build_counts_and_metadata(
            df_tumor=df_tumor,
            df_normal=df_normal,
            how=merge_how,
        )

        tmpdir_obj = tempfile.TemporaryDirectory()
        tmpdir = Path(tmpdir_obj.name)

        try:
            counts_file = tmpdir / "counts.tsv"
            meta_file = tmpdir / "meta.tsv"
            out_file = tmpdir / "deg_results.tsv"

            counts_df.to_csv(counts_file, sep="\t", index=False)
            meta_df.to_csv(meta_file, sep="\t", index=False)

            cmd = [
                "conda", "run", "-n", conda_env,
                "Rscript", str(self.rscript_path),
                "--counts", str(counts_file),
                "--meta", str(meta_file),
                "--out", str(out_file),
                "--method", method,
            "--manual-dispersion", str(manual_dispersion),
            "--min-total-count", str(min_total_count),
            ]

            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )

            if proc.returncode != 0:
                raise RuntimeError(
                    "R DEG script failed.\n"
                    f"Command: {' '.join(cmd)}\n\n"
                    f"STDOUT:\n{proc.stdout.strip()}\n\n"
                    f"STDERR:\n{proc.stderr.strip()}"
                )

            if not out_file.exists():
                raise RuntimeError(
                    "R script finished without error, but the output file was not created: "
                    f"{out_file}"
                )

            deg = pd.read_csv(out_file, sep="\t")

            if keep_temp:
                print(f"Temporary DEG files kept at: {tmpdir}")
                tmpdir_obj = None  # prevent cleanup below

            return deg

        finally:
            if (tmpdir_obj is not None) and (not keep_temp):
                tmpdir_obj.cleanup()