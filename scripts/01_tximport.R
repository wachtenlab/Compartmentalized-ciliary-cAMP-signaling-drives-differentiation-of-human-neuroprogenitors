#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tximport)
})
source("scripts/00_utils.R")

metadata <- read_sample_metadata()
quant_dir <- repo_path("data", "salmon_quant")
gtf_path <- repo_path("data", "reference", "gencode.v49.annotation.gtf.gz")
files <- file.path(quant_dir, metadata$sample, "quant.sf")
names(files) <- metadata$sample
require_files(c(gtf_path, files))

gtf <- fread(
  cmd = paste("gzip -dc", shQuote(gtf_path)),
  sep = "\t",
  header = FALSE,
  quote = "",
  data.table = TRUE
)
setnames(
  gtf,
  c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
)

extract_attribute <- function(x, key) {
  pattern <- paste0(".*", key, " \\\"([^\\\"]+)\\\".*")
  sub(pattern, "\\1", x)
}

tx2gene <- unique(gtf[
  feature == "transcript",
  .(
    transcript_id = extract_attribute(attribute, "transcript_id"),
    gene_id = extract_attribute(attribute, "gene_name")
  )
])

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreAfterBar = TRUE
)

output_dir <- ensure_directory(repo_path("results", "tximport"))
saveRDS(txi, file.path(output_dir, "tximport_gene_level.rds"))
fwrite(
  data.table(gene_id = rownames(txi$counts), as.data.table(txi$counts)),
  file.path(output_dir, "gene_counts.csv")
)

message("Imported ", length(files), " Salmon quantifications and ", nrow(txi$counts), " genes.")

