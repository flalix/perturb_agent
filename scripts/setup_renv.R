#!/usr/bin/env Rscript

# scripts/setup_renv.R

install.packages("renv", repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "limma",
  "apeglm",
  "tximport",
  "biomaRt",
  "optparse",
  "data.table"
))

renv::snapshot()


