#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})
source("scripts/00_utils.R")

dds_path <- repo_path("results", "deseq2", "dds.rds")
require_files(dds_path)
dds <- readRDS(dds_path)
coefficient_names <- resultsNames(dds)
alpha <- 0.05
timepoints <- c("0", "8", "24", "72")
conditions <- c("WT", "BBS1", "Dark", "Light")
output_dir <- ensure_directory(repo_path("results", "deseq2", "contrasts"))

time_coefficient <- function(tp) {
  if (tp == "0") return(character())
  name <- paste0("timepoint_", tp, "_vs_0")
  if (!(name %in% coefficient_names)) stop("Missing coefficient: ", name)
  name
}

condition_coefficient <- function(condition) {
  if (condition == "WT") return(character())
  name <- paste0("condition_", condition, "_vs_WT")
  if (!(name %in% coefficient_names)) stop("Missing coefficient: ", name)
  name
}

interaction_coefficient <- function(tp, condition) {
  if (tp == "0" || condition == "WT") return(character())
  name <- paste0("timepoint", tp, ".condition", condition)
  if (!(name %in% coefficient_names)) stop("Missing coefficient: ", name)
  name
}

condition_at_time <- function(condition, tp) {
  c(condition_coefficient(condition), interaction_coefficient(tp, condition))
}

time_from_t0 <- function(condition, tp) {
  c(time_coefficient(tp), interaction_coefficient(tp, condition))
}

contrast_vector <- function(positive, negative = character()) {
  vector <- setNames(rep(0, length(coefficient_names)), coefficient_names)
  vector[positive] <- vector[positive] + 1
  vector[negative] <- vector[negative] - 1
  vector
}

write_result <- function(id, group, vector, fields) {
  result <- as.data.table(
    as.data.frame(results(dds, contrast = vector, alpha = alpha)),
    keep.rownames = "gene_id"
  )
  for (name in names(fields)) result[, (name) := fields[[name]]]
  result[, `:=`(contrast_id = id, contrast_group = group)]
  setcolorder(
    result,
    c("contrast_id", "contrast_group", names(fields), "gene_id", "baseMean",
      "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  )
  fwrite(result, file.path(output_dir, paste0(id, ".tsv")), sep = "\t")
  result
}

summaries <- list()
index <- 0L

for (tp in timepoints) {
  for (pair in combn(conditions, 2, simplify = FALSE)) {
    comparison <- paste(pair, collapse = "_vs_")
    id <- paste0("condition_", comparison, "_at_T", tp)
    result <- write_result(
      id,
      "condition_within_timepoint",
      contrast_vector(condition_at_time(pair[[1]], tp), condition_at_time(pair[[2]], tp)),
      list(comparison = comparison, timepoint = paste0("T", tp))
    )
    index <- index + 1L
    summaries[[index]] <- data.table(
      contrast_id = id,
      comparison = comparison,
      timepoint = paste0("T", tp),
      n_tested = sum(!is.na(result$pvalue)),
      n_significant = sum(!is.na(result$padj) & result$padj < alpha)
    )
  }
}

for (condition in conditions) {
  for (tp in timepoints[timepoints != "0"]) {
    id <- paste0("time_T", tp, "_vs_T0_in_", condition)
    result <- write_result(
      id,
      "time_within_condition",
      contrast_vector(time_from_t0(condition, tp)),
      list(comparison = paste0("T", tp, "_vs_T0"), condition = condition)
    )
    index <- index + 1L
    summaries[[index]] <- data.table(
      contrast_id = id,
      comparison = paste0("T", tp, "_vs_T0"),
      condition = condition,
      n_tested = sum(!is.na(result$pvalue)),
      n_significant = sum(!is.na(result$padj) & result$padj < alpha)
    )
  }
}

fwrite(rbindlist(summaries, fill = TRUE), file.path(output_dir, "contrast_summary.tsv"), sep = "\t")
message("Exported ", length(summaries), " contrasts.")

