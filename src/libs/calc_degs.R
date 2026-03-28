#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--counts", type = "character", help = "Counts matrix TSV"),
  make_option("--meta", type = "character", help = "Sample metadata TSV"),
  make_option("--out", type = "character", help = "Output TSV"),
  make_option("--method", type = "character", default = "auto",
              help = "auto | deseq2 | edger"),
  make_option("--manual-dispersion", type = "double", default = 0.1,
              help = "Manual dispersion for edgeR when replication is insufficient [default %default]"),
  make_option("--min-total-count", type = "integer", default = 10,
              help = "Filter genes with total counts < this threshold [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$counts) || is.null(opt$meta) || is.null(opt$out)) {
  stop("Required arguments: --counts, --meta, --out")
}

# -----------------------------
# Read input
# -----------------------------
filename_count = opt$counts
filename_metadata = opt$meta

counts_df <- fread(filename_count, sep = "\t", header = TRUE, data.table = FALSE)
meta <- fread(filename_metadata, sep = "\t", header = TRUE, data.table = FALSE)

required_gene_cols <- c("gene_id", "symbol", "gene_type")

missing_gene_cols <- setdiff(required_gene_cols, colnames(counts_df))

if (length(missing_gene_cols) > 0) {
  stop(sprintf("Counts table missing required columns: %s",
               paste(missing_gene_cols, collapse = ", ")))
}

if (!all(c("sample", "condition") %in% colnames(meta))) {
  stop("Metadata must contain columns: sample, condition")
}

sample_cols <- meta$sample
missing_samples <- setdiff(sample_cols, colnames(counts_df))
if (length(missing_samples) > 0) {
  stop(sprintf("These samples from metadata are missing in counts matrix: %s",
               paste(missing_samples, collapse = ", ")))
}

# preserve gene annotations
gene_annot <- counts_df[, required_gene_cols, drop = FALSE]
count_mat <- counts_df[, sample_cols, drop = FALSE]

# coerce to numeric integer-like matrix
for (j in seq_along(count_mat)) {
  count_mat[[j]] <- as.numeric(count_mat[[j]])
}
count_mat[is.na(count_mat)] <- 0

rownames(count_mat) <- gene_annot$gene_id
count_mat <- as.matrix(count_mat)

# -----------------------------
# Basic checks
# -----------------------------
meta$condition <- factor(meta$condition)
if (!all(levels(meta$condition) %in% c("normal", "tumor"))) {
  # relevel later if both present
  meta$condition <- factor(as.character(meta$condition))
}

if (length(unique(meta$condition)) != 2) {
  stop("Exactly two conditions are required.")
}

if (!("normal" %in% meta$condition) || !("tumor" %in% meta$condition)) {
  stop("Conditions must include both 'normal' and 'tumor'.")
}

meta$condition <- relevel(meta$condition, ref = "normal")

group_sizes <- table(meta$condition)
n_normal <- ifelse("normal" %in% names(group_sizes), group_sizes[["normal"]], 0)
n_tumor  <- ifelse("tumor"  %in% names(group_sizes), group_sizes[["tumor"]], 0)

# filter low counts
keep <- rowSums(count_mat) >= opt$`min-total-count`
count_mat  <- count_mat[keep, , drop = FALSE]
gene_annot <- gene_annot[keep, , drop = FALSE]

if (nrow(count_mat) == 0) {
  stop("No genes remain after filtering.")
}

# -----------------------------
# Helper: edgeR
# -----------------------------
run_edger <- function(count_mat, meta, gene_annot, out_file, manual_dispersion = 0.1) {
  suppressPackageStartupMessages(library(edgeR))

  y <- DGEList(counts = count_mat, group = meta$condition)
  keep <- filterByExpr(y, group = meta$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  gene_sub <- gene_annot[match(rownames(y$counts), gene_annot$gene_id), , drop = FALSE]

  design <- model.matrix(~ condition, data = meta)

  # Replicates available in both groups -> estimate dispersion normally
  group_sizes <- table(meta$condition)
  enough_rep <- all(group_sizes >= 2)

  if (enough_rep) {
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    tt <- topTags(qlf, n = Inf, sort.by = "PValue")$table

    res <- data.frame(
      gene_id = rownames(tt),
      log2FoldChange = tt$logFC,
      logCPM = tt$logCPM,
      statistic = tt$F,
      pvalue = tt$PValue,
      padj = tt$FDR,
      method = "edgeR_QLF",
      stringsAsFactors = FALSE
    )
  } else {
    # No proper replication -> manual dispersion approximation
    et <- exactTest(y, dispersion = manual_dispersion)
    tt <- topTags(et, n = Inf, sort.by = "PValue")$table

    res <- data.frame(
      gene_id = rownames(tt),
      log2FoldChange = tt$logFC,
      logCPM = tt$logCPM,
      statistic = NA_real_,
      pvalue = tt$PValue,
      padj = p.adjust(tt$PValue, method = "BH"),
      method = sprintf("edgeR_exact_manualDisp_%s", manual_dispersion),
      stringsAsFactors = FALSE
    )
  }

  res <- merge(gene_sub, res, by = "gene_id", all.y = TRUE, sort = FALSE)
  res <- res[order(res$padj, res$pvalue), ]
  fwrite(res, file = out_file, sep = "\t", quote = FALSE, na = "NA")
}

# -----------------------------
# Helper: DESeq2
# -----------------------------
run_deseq2 <- function(count_mat, meta, gene_annot, out_file) {
  suppressPackageStartupMessages(library(DESeq2))

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = meta,
    design = ~ condition
  )

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  dds <- DESeq(dds)
  rr <- results(dds, contrast = c("condition", "tumor", "normal"))

  res <- data.frame(
    gene_id = rownames(rr),
    log2FoldChange = rr$log2FoldChange,
    lfcSE = rr$lfcSE,
    statistic = rr$stat,
    pvalue = rr$pvalue,
    padj = rr$padj,
    baseMean = rr$baseMean,
    method = "DESeq2",
    stringsAsFactors = FALSE
  )

  gene_sub <- gene_annot[match(res$gene_id, gene_annot$gene_id), , drop = FALSE]
  res <- cbind(gene_sub, res[, setdiff(colnames(res), "gene_id"), drop = FALSE])
  res <- res[order(res$padj, res$pvalue), ]
  fwrite(res, file = out_file, sep = "\t", quote = FALSE, na = "NA")
}

# -----------------------------
# Method selection
# -----------------------------
chosen_method <- opt$method
if (chosen_method == "auto") {
  # Prefer DESeq2 only when both groups have at least 2 samples
  if (n_normal >= 2 && n_tumor >= 2) {
    chosen_method <- "deseq2"
  } else {
    chosen_method <- "edger"
  }
}

if (chosen_method == "deseq2") {
  if (n_normal < 2 || n_tumor < 2) {
    stop("DESeq2 requires replication support. Use method='edger' or method='auto' for sparse designs.")
  }
  run_deseq2(count_mat, meta, gene_annot, opt$out)

} else if (chosen_method == "edger") {
  run_edger(count_mat, meta, gene_annot, opt$out, manual_dispersion = opt$`manual-dispersion`)

} else {
  stop("Invalid method. Use: auto, deseq2, edger")
}