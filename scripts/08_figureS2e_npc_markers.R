#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("scripts/00_utils.R")

marker_path <- repo_path("data", "annotations", "neurogenic_progenitor_markers.tsv")
require_files(marker_path)
markers <- fread(marker_path, sep = "\t")
markers[, gene := toupper(gene)]
npc_order <- unique(markers[set == "Progenitor_maintenance", gene])
timepoints <- c("T0", "T8", "T24", "T72")

results <- read_light_vs_dark(timepoints)
npc <- results[gene %in% npc_order, .(timepoint, gene, log2FC_light_vs_dark, padj)]
npc[, significant := !is.na(padj) & padj < 0.05]
npc[, light_down := significant & log2FC_light_vs_dark < 0]
npc[, neglog10_padj := fifelse(significant, -log10(pmax(padj, 1e-300)), NA_real_)]
selection <- npc[, .(n_light_down = sum(light_down)), by = gene]
selected <- selection[n_light_down >= 3L, gene]
selected <- npc_order[npc_order %in% selected]

plot_data <- npc[gene %in% selected]
fwrite(
  plot_data[order(match(gene, selected), match(timepoint, timepoints))],
  repo_path("data", "figure_data", "FigureS2e_npc_markers.csv")
)
plot_data[, timepoint := factor(timepoint, levels = timepoints)]
plot_data[, gene := factor(gene, levels = rev(selected))]
points <- plot_data[significant == TRUE]
colour_limit <- ceiling(max(abs(points$log2FC_light_vs_dark), na.rm = TRUE) * 2) / 2
size_limit <- ceiling(max(points$neglog10_padj, na.rm = TRUE) / 5) * 5

plot <- ggplot(plot_data, aes(timepoint, gene)) +
  geom_tile(fill = "white", colour = "black", linewidth = 0.35) +
  geom_point(
    data = points,
    aes(fill = log2FC_light_vs_dark, size = neglog10_padj),
    shape = 21, colour = "black", stroke = 0.35
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
height <- max(6.5, length(selected) * 0.42 + 2.2)
ggsave(repo_path("figures", "FigureS2e_npc_markers.pdf"), plot, width = 5.6, height = height)
ggsave(repo_path("figures", "FigureS2e_npc_markers.png"), plot, width = 5.6, height = height, dpi = 400)
message("Supplementary Figure S2e and its source-data table were written.")
