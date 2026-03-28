from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Literal

import pandas as pd


GENE_COLS = ["gene_id", "symbol", "gene_type"]


def _find_count_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in GENE_COLS]


def _validate_expression_df(df: pd.DataFrame, name: str) -> None:
    missing = [c for c in GENE_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")

    count_cols = _find_count_columns(df)
    if not count_cols:
        raise ValueError(f"{name} has no count columns besides {GENE_COLS}")

    if df["gene_id"].duplicated().any():
        dup_n = int(df["gene_id"].duplicated().sum())
        raise ValueError(f"{name} has duplicated gene_id values ({dup_n} duplicates).")


def _rename_count_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    count_cols = _find_count_columns(df)
    rename_map = {old: f"{prefix}_{i+1}" for i, old in enumerate(count_cols)}
    return df.rename(columns=rename_map).copy()


def build_counts_and_metadata(
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
    _validate_expression_df(df_tumor, "df_tumor")
    _validate_expression_df(df_normal, "df_normal")

    tumor = _rename_count_columns(df_tumor, "tumor")
    normal = _rename_count_columns(df_normal, "normal")

    merged = pd.merge(
        normal,
        tumor,
        on=GENE_COLS,
        how=how,
        validate="one_to_one",
    )

    count_cols = [c for c in merged.columns if c not in GENE_COLS]
    if not count_cols:
        raise ValueError("No count columns after merging tumor and normal dataframes.")

    for c in count_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0)

    # round and convert to int for RNA-seq counts
    merged[count_cols] = merged[count_cols].round().astype(int)

    sample_cols = [c for c in merged.columns if c.startswith("normal_")] + [
        c for c in merged.columns if c.startswith("tumor_")
    ]

    meta = pd.DataFrame({
        "sample": sample_cols,
        "condition": ["normal"] * sum(c.startswith("normal_") for c in sample_cols)
                    + ["tumor"] * sum(c.startswith("tumor_") for c in sample_cols),
    })

    counts = merged[GENE_COLS + sample_cols].copy()
    return counts, meta


def run_deg_rscript(
    df_tumor: pd.DataFrame,
    df_normal: pd.DataFrame,
    method: Literal["auto", "deseq2", "edger"] = "auto",
    manual_dispersion: float = 0.1,
    min_total_count: int = 10,
    merge_how: Literal["inner", "outer"] = "inner",
    keep_temp: bool = False,
) -> pd.DataFrame:
    """
    Run DEG analysis in R using DESeq2 or edgeR.

    Returns
    -------
    pd.DataFrame
        DEG results table
    """
    if shutil.which("Rscript") is None:
        raise EnvironmentError("Rscript was not found in PATH.")

    rscript_path = Path("./calc_degs.R")
    if not rscript_path.exists():
        raise FileNotFoundError(f"R script not found: {rscript_path}")

    counts_df, meta_df = build_counts_and_metadata(df_tumor, df_normal, how=merge_how)

    tmpdir_obj = tempfile.TemporaryDirectory()
    tmpdir = Path(tmpdir_obj.name)

    counts_file = tmpdir / "counts.tsv"
    meta_file = tmpdir / "meta.tsv"
    out_file = tmpdir / "deg_results.tsv"

    counts_df.to_csv(counts_file, sep="\t", index=False)
    meta_df.to_csv(meta_file, sep="\t", index=False)

    cmd = [
        "Rscript",
        str(rscript_path),
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
        stderr = proc.stderr.strip()
        stdout = proc.stdout.strip()
        msg = (
            f"R DEG script failed with exit code {proc.returncode}\n"
            f"STDOUT:\n{stdout}\n\nSTDERR:\n{stderr}"
        )
        if not keep_temp:
            tmpdir_obj.cleanup()
        raise RuntimeError(msg)

    if not out_file.exists():
        if not keep_temp:
            tmpdir_obj.cleanup()
        raise RuntimeError("R script finished but output file was not created.")

    deg = pd.read_csv(out_file, sep="\t")

    if keep_temp:
        print(f"Temporary DEG files kept at: {tmpdir}")
    else:
        tmpdir_obj.cleanup()

    return deg