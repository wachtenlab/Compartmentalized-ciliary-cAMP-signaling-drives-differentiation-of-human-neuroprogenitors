#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("scripts/00_utils.R")

timepoints <- c("T0", "T8", "T24", "T72")
gene_path <- repo_path("data", "annotations", "figure3b_genes.tsv")
require_files(gene_path)
genes <- fread(gene_path)$gene
results <- read_light_vs_dark(timepoints)

plot_data <- results[gene %in% genes, .(
  timepoint,
  gene,
  log2FC_light_vs_dark,
  padj,
  significant = !is.na(padj) & padj < 0.05
)]
plot_data[, neglog10_padj := fifelse(significant, -log10(pmax(padj, 1e-300)), NA_real_)]
plot_data[, timepoint := factor(timepoint, levels = timepoints)]
plot_data[, gene := factor(gene, levels = rev(genes))]

figure_data <- copy(plot_data)
figure_data[, `:=`(timepoint = as.character(timepoint), gene = as.character(gene))]
fwrite(
  figure_data[order(match(gene, genes), match(timepoint, timepoints))],
  repo_path("data", "figure_data", "Figure3b_neurogenic_markers.tsv"),
  sep = "\t"
)

points <- plot_data[significant == TRUE]
colour_limit <- ceiling(max(abs(points$log2FC_light_vs_dark), na.rm = TRUE) * 2) / 2
size_limit <- ceiling(max(points$neglog10_padj, na.rm = TRUE) / 10) * 10

plot <- ggplot(plot_data, aes(timepoint, gene)) +
  geom_tile(fill = "white", colour = "black", linewidth = 0.35) +
  geom_point(
    data = points,
    aes(fill = log2FC_light_vs_dark, size = neglog10_padj),
    shape = 21,
    colour = "#303030",
    stroke = 0.35
  ) +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#C51B2D", midpoint = 0,
    limits = c(-colour_limit, colour_limit), name = "log2FC"
  ) +
  scale_size_continuous(range = c(1.8, 7.5), limits = c(0, size_limit), name = "-log10(padj)") +
  coord_fixed() +
  labs(x = "Timepoint (h)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.position = "right")

ensure_directory(repo_path("figures"))
ggsave(repo_path("figures", "Figure3b_neurogenic_markers.pdf"), plot, width = 5.2, height = 5.8)
ggsave(repo_path("figures", "Figure3b_neurogenic_markers.png"), plot, width = 5.2, height = 5.8, dpi = 400)
message("Figure 3b and its source-data table were written.")
