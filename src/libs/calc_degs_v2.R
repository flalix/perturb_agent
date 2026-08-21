#!/usr/bin/env Rscript
#
# calc_degs_v2.R
# --------------
# Adds limma to the existing DESeq2 / edgeR paths.
#
#   --method limma_trend   log-expression input (BayesPrism psi). eBayes(trend=TRUE).
#   --method limma_voom    integer counts. voom mean-variance weights.
#   --method deseq2        integer counts (unchanged)
#   --method edger         integer counts (unchanged)
#
# Why limma-trend for psi
# -----------------------
# DESeq2 and edgeR model counts with a negative binomial. BayesPrism psi is a
# posterior expression estimate, not an observed count -- rounding it to feed a
# NB treats an estimate as data and understates uncertainty. limma-trend works on
# log-expression with an intensity-dependent variance trend, which is the right
# shape for psi.
#
# Caveat that applies to every method here: none propagates deconvolution
# uncertainty. psi enters as if measured exactly, so p-values are anticonservative.
# Weight effect size accordingly.
#
# New options
# -----------
#   --already-log        input is already log-scale (skip the log2 transform)
#   --normalize          none | quantile | scale   (between-sample, trend mode)
#   --covariates         comma-separated meta columns to add to the design
#   --pair-col           meta column identifying matched pairs -> ~ pair + condition
#   --min-prop-detected  detection filter for continuous input (replaces
#                        --min-total-count, which is meaningless on a psi scale)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--counts", type = "character", help = "Expression/counts matrix TSV"),
  make_option("--meta", type = "character", help = "Sample metadata TSV"),
  make_option("--out", type = "character", help = "Output TSV"),
  make_option("--method", type = "character", default = "auto",
              help = "auto | deseq2 | edger | limma_trend | limma_voom"),
  make_option("--manual-dispersion", type = "double", default = 0.1,
              help = "Manual dispersion for edgeR when replication is insufficient"),
  make_option("--min-total-count", type = "integer", default = 10,
              help = "Count-mode filter: drop genes with total counts below this"),
  make_option("--min-prop-detected", type = "double", default = 0.20,
              help = "Continuous-mode filter: min proportion of detected (>0) samples in at least one group"),
  make_option("--already-log", action = "store_true", default = FALSE,
              help = "Input matrix is already on a log scale"),
  make_option("--normalize", type = "character", default = "none",
              help = "none | quantile | scale  (limma_trend only)"),
  make_option("--covariates", type = "character", default = NULL,
              help = "Comma-separated meta columns to include in the design"),
  make_option("--pair-col", type = "character", default = NULL,
              help = "Meta column of matched-pair IDs -> ~ pair + condition"),
  make_option("--prior-count", type = "double", default = 0.5,
              help = "Offset added before log2 in continuous mode"),
  make_option("--array-weights", action = "store_true", default = FALSE,
              help = "Estimate per-sample weights (limma::arrayWeights) and use them in lmFit.
                      Use when sample precision varies -- e.g. deconvolved psi where noise
                      scales with 1/theta, so low-abundance samples are noisier.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$counts) || is.null(opt$meta) || is.null(opt$out)) {
  stop("Required arguments: --counts, --meta, --out")
}

counts_df <- fread(opt$counts, sep = "\t", header = TRUE, data.table = FALSE)
meta      <- fread(opt$meta,   sep = "\t", header = TRUE, data.table = FALSE)

if (!("geneid" %in% colnames(counts_df))) stop("Counts table missing column: geneid")
if (!all(c("sample", "condition") %in% colnames(meta)))
  stop("Metadata must contain columns: sample, condition")

sample_cols <- meta$sample
missing_samples <- setdiff(sample_cols, colnames(counts_df))
if (length(missing_samples) > 0) {
  stop(sprintf("Samples in metadata missing from matrix: %s",
               paste(missing_samples, collapse = ", ")))
}

gene_annot <- counts_df[, "geneid", drop = FALSE]
mat <- counts_df[, sample_cols, drop = FALSE]
for (j in seq_along(mat)) mat[[j]] <- as.numeric(mat[[j]])
mat[is.na(mat)] <- 0
rownames(mat) <- gene_annot$geneid
mat <- as.matrix(mat)

meta$condition <- factor(as.character(meta$condition))
if (length(levels(meta$condition)) != 2) stop("Exactly two conditions are required.")
if (!all(c("normal", "tumor") %in% levels(meta$condition)))
  stop("Conditions must include both 'normal' and 'tumor'.")
meta$condition <- relevel(meta$condition, ref = "normal")

group_sizes <- table(meta$condition)
n_normal <- as.integer(group_sizes[["normal"]])
n_tumor  <- as.integer(group_sizes[["tumor"]])

chosen_method <- opt$method
if (chosen_method == "auto") {
  chosen_method <- if (n_normal >= 2 && n_tumor >= 2) "deseq2" else "edger"
}
continuous_mode <- chosen_method == "limma_trend"

# ---------------------------------------------------------------- filtering
if (continuous_mode) {
  det <- function(idx) rowMeans(mat[, idx, drop = FALSE] > 0)
  keep <- (det(meta$condition == "normal") >= opt$`min-prop-detected`) |
          (det(meta$condition == "tumor")  >= opt$`min-prop-detected`)
} else {
  keep <- rowSums(mat) >= opt$`min-total-count`
}
mat        <- mat[keep, , drop = FALSE]
gene_annot <- gene_annot[keep, , drop = FALSE]
if (nrow(mat) == 0) stop("No genes remain after filtering.")

cat(sprintf("[calc_degs] method=%s  genes=%d  normal=%d  tumor=%d\n",
            chosen_method, nrow(mat), n_normal, n_tumor))

# ---------------------------------------------------------------- design
build_design <- function(meta, covariates = NULL, pair_col = NULL) {
  terms <- c()
  if (!is.null(pair_col) && nzchar(pair_col)) {
    if (!(pair_col %in% colnames(meta)))
      stop(sprintf("pair-col '%s' not found in metadata", pair_col))
    meta[[pair_col]] <- factor(as.character(meta[[pair_col]]))
    # a pair level present in only one condition contributes nothing and makes
    # the design rank-deficient; drop those samples upstream if this triggers
    tab <- table(meta[[pair_col]], meta$condition)
    incomplete <- rownames(tab)[rowSums(tab > 0) < 2]
    if (length(incomplete) > 0)
      warning(sprintf("%d pair level(s) present in only one condition: %s",
                      length(incomplete), paste(head(incomplete, 5), collapse = ", ")))
    terms <- c(terms, pair_col)
  }
  if (!is.null(covariates) && nzchar(covariates)) {
    cv <- trimws(strsplit(covariates, ",")[[1]])
    miss <- setdiff(cv, colnames(meta))
    if (length(miss) > 0)
      stop(sprintf("covariates not in metadata: %s", paste(miss, collapse = ", ")))
    for (c0 in cv) if (is.character(meta[[c0]])) meta[[c0]] <- factor(meta[[c0]])
    terms <- c(terms, cv)
  }
  terms <- c(terms, "condition")
  form <- as.formula(paste("~", paste(terms, collapse = " + ")))
  design <- model.matrix(form, data = meta)
  if (qr(design)$rank < ncol(design))
    stop("Design matrix is rank-deficient. A covariate is likely confounded with condition -- ",
         "e.g. tumours from one cohort and normals from another cannot be adjusted for cohort.")
  list(design = design, formula = form)
}

coef_name <- "conditiontumor"

# ---------------------------------------------------------------- limma
run_limma <- function(mat, meta, gene_annot, out_file, mode = "trend") {
  suppressPackageStartupMessages(library(limma))

  bd <- build_design(meta, opt$covariates, opt$`pair-col`)
  design <- bd$design
  cat(sprintf("[calc_degs] design: %s\n", deparse1(bd$formula)))
  if (!(coef_name %in% colnames(design)))
    stop(sprintf("coefficient '%s' not in design columns: %s",
                 coef_name, paste(colnames(design), collapse = ", ")))

  if (mode == "voom") {
    suppressPackageStartupMessages(library(edgeR))
    y <- DGEList(counts = round(mat))
    y <- calcNormFactors(y)
    v <- voom(y, design, plot = FALSE)
    wts <- NULL
    if (opt$`array-weights`) {
      wts <- arrayWeights(v, design)
      cat(sprintf("[calc_degs] array weights: min=%.3f median=%.3f max=%.3f\n",
                  min(wts), median(wts), max(wts)))
    }
    fit <- eBayes(lmFit(v, design, weights = wts), robust = TRUE)
    method_lab <- if (opt$`array-weights`) "limma_voom_aw" else "limma_voom"
    ave <- fit$Amean
  } else {
    E <- if (opt$`already-log`) mat else log2(mat + opt$`prior-count`)
    if (opt$normalize %in% c("quantile", "scale")) {
      E <- normalizeBetweenArrays(E, method = opt$normalize)
      method_lab <- sprintf("limma_trend_%s", opt$normalize)
    } else {
      method_lab <- "limma_trend"
    }
    wts <- NULL
    if (opt$`array-weights`) {
      wts <- arrayWeights(E, design)
      method_lab <- paste0(method_lab, "_aw")
      cat(sprintf("[calc_degs] array weights: min=%.3f median=%.3f max=%.3f\n",
                  min(wts), median(wts), max(wts)))
      wr <- data.frame(sample = colnames(E), weight = as.numeric(wts))
      wr <- wr[order(wr$weight), ]
      cat("[calc_degs] five lowest-weight samples:\n")
      print(utils::head(wr, 5), row.names = FALSE)
    }
    fit <- eBayes(lmFit(E, design, weights = wts), trend = TRUE, robust = TRUE)
    ave <- fit$Amean
  }

  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")

  res <- data.frame(
    geneid         = rownames(tt),
    log2FoldChange = tt$logFC,
    AveExpr        = tt$AveExpr,
    statistic      = tt$t,
    pvalue         = tt$P.Value,
    padj           = tt$adj.P.Val,
    B              = tt$B,
    method         = method_lab,
    stringsAsFactors = FALSE
  )
  gene_sub <- gene_annot[match(res$geneid, gene_annot$geneid), , drop = FALSE]
  res <- cbind(gene_sub[, setdiff(colnames(gene_sub), "geneid"), drop = FALSE], res)
  res <- res[order(res$padj, res$pvalue), ]
  fwrite(res, file = out_file, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("[calc_degs] wrote %d rows (padj<0.05: %d)\n",
              nrow(res), sum(res$padj < 0.05, na.rm = TRUE)))
}

# ---------------------------------------------------------------- edgeR
run_edger <- function(count_mat, meta, gene_annot, out_file, manual_dispersion = 0.1) {
  suppressPackageStartupMessages(library(edgeR))
  y <- DGEList(counts = round(count_mat), group = meta$condition)
  k <- filterByExpr(y, group = meta$condition)
  y <- y[k, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  gene_sub <- gene_annot[match(rownames(y$counts), gene_annot$geneid), , drop = FALSE]
  design <- build_design(meta, opt$covariates, opt$`pair-col`)$design

  if (all(table(meta$condition) >= 2)) {
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = coef_name)
    tt <- topTags(qlf, n = Inf, sort.by = "PValue")$table
    res <- data.frame(geneid = rownames(tt), log2FoldChange = tt$logFC,
                      logCPM = tt$logCPM, statistic = tt$F, pvalue = tt$PValue,
                      padj = tt$FDR, method = "edgeR_QLF", stringsAsFactors = FALSE)
  } else {
    et <- exactTest(y, dispersion = manual_dispersion)
    tt <- topTags(et, n = Inf, sort.by = "PValue")$table
    res <- data.frame(geneid = rownames(tt), log2FoldChange = tt$logFC,
                      logCPM = tt$logCPM, statistic = NA_real_, pvalue = tt$PValue,
                      padj = p.adjust(tt$PValue, method = "BH"),
                      method = sprintf("edgeR_exact_manualDisp_%s", manual_dispersion),
                      stringsAsFactors = FALSE)
  }
  res <- merge(gene_sub, res, by = "geneid", all.y = TRUE, sort = FALSE)
  res <- res[order(res$padj, res$pvalue), ]
  fwrite(res, file = out_file, sep = "\t", quote = FALSE, na = "NA")
}

# ---------------------------------------------------------------- DESeq2
run_deseq2 <- function(count_mat, meta, gene_annot, out_file) {
  suppressPackageStartupMessages(library(DESeq2))
  bd <- build_design(meta, opt$covariates, opt$`pair-col`)
  # build_design factors these in its own local copy only; DESeq2 reads colData
  # directly, so convert here as well or a character pair column silently
  # becomes an unordered factor with a different reference level.
  for (v in all.vars(bd$formula)) {
    if (v %in% colnames(meta) && !is.factor(meta[[v]]) && !is.numeric(meta[[v]]))
      meta[[v]] <- factor(as.character(meta[[v]]))
  }
  dds <- DESeqDataSetFromMatrix(countData = round(count_mat), colData = meta,
                                design = bd$formula)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  rr <- results(dds, contrast = c("condition", "tumor", "normal"))
  res <- data.frame(geneid = rownames(rr), log2FoldChange = rr$log2FoldChange,
                    lfcSE = rr$lfcSE, statistic = rr$stat, pvalue = rr$pvalue,
                    padj = rr$padj, baseMean = rr$baseMean, method = "DESeq2",
                    stringsAsFactors = FALSE)
  gene_sub <- gene_annot[match(res$geneid, gene_annot$geneid), , drop = FALSE]
  res <- cbind(gene_sub[, setdiff(colnames(gene_sub), "geneid"), drop = FALSE], res)
  res <- res[order(res$padj, res$pvalue), ]
  fwrite(res, file = out_file, sep = "\t", quote = FALSE, na = "NA")
}

# ---------------------------------------------------------------- dispatch
if (chosen_method == "limma_trend") {
  run_limma(mat, meta, gene_annot, opt$out, mode = "trend")
} else if (chosen_method == "limma_voom") {
  run_limma(mat, meta, gene_annot, opt$out, mode = "voom")
} else if (chosen_method == "deseq2") {
  if (n_normal < 2 || n_tumor < 2)
    stop("DESeq2 requires replication. Use edger or limma_trend.")
  run_deseq2(mat, meta, gene_annot, opt$out)
} else if (chosen_method == "edger") {
  run_edger(mat, meta, gene_annot, opt$out, manual_dispersion = opt$`manual-dispersion`)
} else {
  stop("Invalid method. Use: auto, deseq2, edger, limma_trend, limma_voom")
}
