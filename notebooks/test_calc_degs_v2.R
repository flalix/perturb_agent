#!/usr/bin/env Rscript
#
# test_calc_degs_v2.R
# -------------------
# Verify the environment and exercise every method in calc_degs_v2.R on synthetic
# data before pointing it at real psi. Run inside the renv conda environment:
#
#   conda activate renv
#   Rscript test_calc_degs_v2.R /path/to/calc_degs_v2.R
#
# Checks, in order:
#   1. R version and package availability
#   2. limma_trend on continuous (psi-like) input, with 200 planted DEGs
#   3. limma_trend with a paired design
#   4. limma_voom / edgeR / DESeq2 on integer counts
#   5. that a confounded covariate is REFUSED rather than silently fitted

args <- commandArgs(trailingOnly = TRUE)
script <- if (length(args) >= 1) args[1] else "calc_degs_v2.R"
if (!file.exists(script)) stop(sprintf("script not found: %s", script))

cat("=== 1. environment ===\n")
cat(sprintf("R: %s\n", R.version.string))
pkgs <- c("optparse", "data.table", "limma", "edgeR", "DESeq2")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-12s %s%s\n", p, if (ok) "OK  " else "MISSING",
              if (ok) as.character(packageVersion(p)) else ""))
}
if (requireNamespace("BiocManager", quietly = TRUE))
  cat(sprintf("  Bioconductor %s\n", as.character(BiocManager::version())))
if (!all(sapply(pkgs, requireNamespace, quietly = TRUE)))
  stop("install missing packages before continuing")

tmp <- tempfile("degtest"); dir.create(tmp)
set.seed(1)
G <- 3000; nT <- 40; nN <- 20

# ---- synthetic psi-like matrix with 200 planted DEGs ----------------------
base  <- 2^rnorm(G, mean = 4, sd = 1.5)
mkmat <- function(n, mult) {
  m <- matrix(rgamma(G * n, shape = 8, scale = base * mult / 8), nrow = G)
  m[sample.int(length(m), length(m) * 0.05)] <- 0      # dropout
  m
}
mult <- rep(1, G); truth <- sample.int(G, 200)
mult[truth] <- rep(c(4, 0.25), each = 100)

psi <- cbind(mkmat(nN, 1)[, seq_len(nN)], mkmat(nT, 1)[, seq_len(nT)])
psi[, (nN + 1):(nN + nT)] <- psi[, (nN + 1):(nN + nT)] * mult
genes <- sprintf("ENSG%08d", seq_len(G))
cols  <- c(sprintf("N_%02d", seq_len(nN)), sprintf("T_%02d", seq_len(nT)))

write_inputs <- function(mat, meta, tag, integer_mode = FALSE) {
  m <- if (integer_mode) round(mat * 20) else mat
  df <- data.frame(geneid = genes, m, check.names = FALSE)
  colnames(df) <- c("geneid", cols)
  cf <- file.path(tmp, sprintf("counts_%s.tsv", tag))
  mf <- file.path(tmp, sprintf("meta_%s.tsv", tag))
  data.table::fwrite(df, cf, sep = "\t")
  data.table::fwrite(meta, mf, sep = "\t")
  list(counts = cf, meta = mf, out = file.path(tmp, sprintf("res_%s.tsv", tag)))
}

meta_basic <- data.frame(
  sample    = cols,
  condition = c(rep("normal", nN), rep("tumor", nT)),
  cohort    = sample(c("A", "B"), nN + nT, replace = TRUE),
  stringsAsFactors = FALSE
)
# perfectly confounded: every normal in one cohort, every tumour in the other
meta_conf <- meta_basic
meta_conf$cohort <- c(rep("A", nN), rep("B", nT))
# paired: 20 tumours matched to the 20 normals, 20 unmatched
meta_pair <- meta_basic
meta_pair$pair <- c(sprintf("P%02d", seq_len(nN)),
                    sprintf("P%02d", seq_len(nN)),
                    sprintf("U%02d", seq_len(nT - nN)))

run <- function(io, ..., expect_fail = FALSE) {
  a <- c(script, "--counts", io$counts, "--meta", io$meta, "--out", io$out, ...)
  p <- system2("Rscript", a, stdout = TRUE, stderr = TRUE)
  status <- attr(p, "status")
  failed <- !is.null(status) && status != 0
  if (expect_fail) {
    cat(if (failed) "  REFUSED (correct)\n" else "  ACCEPTED -- PROBLEM\n")
    if (failed) cat("   ", tail(p, 2)[1], "\n")
    return(invisible(NULL))
  }
  if (failed) { cat(paste(p, collapse = "\n"), "\n"); stop("run failed") }
  cat(paste0("  ", grep("^\\[calc_degs\\]", p, value = TRUE), collapse = "\n"), "\n")
  res <- data.table::fread(io$out, data.table = FALSE)
  hits <- res$geneid[which(res$padj < 0.05)]
  tp <- sum(hits %in% genes[truth]); fp <- length(hits) - tp
  cat(sprintf("  sig=%d  recovered=%d/200  false=%d\n", length(hits), tp, fp))
  invisible(res)
}

cat("\n=== 2. limma_trend, continuous psi ===\n")
io <- write_inputs(psi, meta_basic, "trend")
run(io, "--method", "limma_trend", "--min-prop-detected", "0.2")

cat("\n=== 3. limma_trend, paired design ===\n")
io <- write_inputs(psi, meta_pair, "pair")
run(io, "--method", "limma_trend", "--pair-col", "pair")

cat("\n=== 4. count-based methods ===\n")
for (m in c("limma_voom", "edger", "deseq2")) {
  cat(sprintf(" -- %s\n", m))
  io <- write_inputs(psi, meta_basic, m, integer_mode = TRUE)
  run(io, "--method", m, "--min-total-count", "10")
}

cat("\n=== 5. confounded covariate must be refused ===\n")
io <- write_inputs(psi, meta_conf, "conf")
run(io, "--method", "limma_trend", "--covariates", "cohort", expect_fail = TRUE)

cat(sprintf("\ntemp dir: %s\n", tmp))
