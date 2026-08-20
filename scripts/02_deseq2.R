#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})
source("scripts/00_utils.R")

metadata <- read_sample_metadata()
metadata$sample <- as.character(metadata$sample)
rownames(metadata) <- metadata$sample
metadata$condition <- factor(metadata$condition)
if (!("WT" %in% levels(metadata$condition))) stop("WT is required as the reference condition.")
metadata$condition <- relevel(metadata$condition, ref = "WT")
metadata$timepoint <- factor(
  sub("^T", "", as.character(metadata$timepoint)),
  levels = c("0", "8", "24", "72")
)
if (anyNA(metadata$timepoint)) stop("Expected timepoints T0, T8, T24, and T72.")

use_batch <- "batch" %in% names(metadata) && length(unique(metadata$batch)) > 1L
if (use_batch) metadata$batch <- factor(metadata$batch)
design_formula <- if (use_batch) {
  ~ batch + timepoint + condition + timepoint:condition
} else {
  ~ timepoint + condition + timepoint:condition
}

txi_candidates <- c(
  repo_path("results", "tximport", "tximport_gene_level.rds"),
  repo_path("data", "counts", "tximport_gene_level.rds")
)
txi_path <- txi_candidates[file.exists(txi_candidates)][1]
if (!is.na(txi_path)) {
  txi <- readRDS(txi_path)
  missing_samples <- setdiff(metadata$sample, colnames(txi$counts))
  if (length(missing_samples) > 0L) {
    stop("Samples missing from tximport object: ", paste(missing_samples, collapse = ", "))
  }
  txi$counts <- txi$counts[, metadata$sample, drop = FALSE]
  txi$abundance <- txi$abundance[, metadata$sample, drop = FALSE]
  txi$length <- txi$length[, metadata$sample, drop = FALSE]
  dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = design_formula)
} else {
  counts_table <- read_gene_counts()
  missing_samples <- setdiff(metadata$sample, names(counts_table))
  if (length(missing_samples) > 0L) {
    stop("Samples missing from count matrix: ", paste(missing_samples, collapse = ", "))
  }
  sample_columns <- metadata$sample
  counts <- as.matrix(counts_table[, ..sample_columns])
  rownames(counts) <- counts_table$gene_id
  storage.mode(counts) <- "numeric"
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = design_formula
  )
}

dds <- DESeq(dds)

output_dir <- ensure_directory(repo_path("results", "deseq2"))
saveRDS(dds, file.path(output_dir, "dds.rds"))
fwrite(
  data.table(gene_id = rownames(dds), as.data.table(counts(dds, normalized = TRUE))),
  file.path(output_dir, "normalized_counts.tsv"),
  sep = "\t"
)
writeLines(resultsNames(dds), file.path(output_dir, "model_coefficients.txt"))
message("DESeq2 model fitted for ", ncol(dds), " samples and ", nrow(dds), " genes.")
