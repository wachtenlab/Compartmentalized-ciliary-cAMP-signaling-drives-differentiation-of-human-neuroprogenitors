suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(
  Sys.getenv("RNA_SEQ_REPOSITORY", unset = "."),
  mustWork = TRUE
)

repo_path <- function(...) file.path(project_root, ...)

ensure_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required files:\n - ", paste(missing, collapse = "\n - "))
  }
  invisible(paths)
}

read_sample_metadata <- function() {
  path <- repo_path("data", "metadata", "sample_metadata.tsv")
  require_files(path)
  metadata <- fread(path, sep = "\t", data.table = FALSE)
  required <- c("sample", "condition", "timepoint")
  missing <- setdiff(required, names(metadata))
  if (length(missing) > 0L) {
    stop("Sample metadata is missing columns: ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(metadata$sample)) stop("Sample identifiers must be unique.")
  metadata
}

read_gene_counts <- function() {
  path <- repo_path("data", "counts", "gene_counts.csv")
  require_files(path)
  counts <- fread(path)
  if (ncol(counts) < 2L) stop("The count matrix must contain a gene column and samples.")
  setnames(counts, 1L, "gene_id")
  if (anyNA(counts$gene_id) || any(counts$gene_id == "")) stop("Gene identifiers must not be missing.")
  counts[, gene_id := make.unique(as.character(gene_id))]
  counts
}

contrast_path <- function(comparison, timepoint) {
  repo_path(
    "results", "deseq2", "contrasts",
    paste0("condition_", comparison, "_at_", timepoint, ".tsv")
  )
}

read_light_vs_dark <- function(timepoints = c("T0", "T8", "T24", "T72")) {
  paths <- vapply(
    timepoints,
    function(tp) contrast_path("Dark_vs_Light", tp),
    character(1)
  )
  require_files(paths)
  rbindlist(lapply(seq_along(paths), function(i) {
    x <- fread(paths[[i]])
    x[, `:=`(
      timepoint = timepoints[[i]],
      gene = toupper(gene_id),
      log2FC_light_vs_dark = -log2FoldChange
    )]
    x
  }), fill = TRUE)
}
