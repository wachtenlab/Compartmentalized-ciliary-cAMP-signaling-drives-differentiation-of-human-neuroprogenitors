#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("scripts/00_utils.R")

count_path <- repo_path("data", "figure_data", "FigureS2b_deg_counts.tsv")
overlap_path <- repo_path("data", "figure_data", "FigureS2c_deg_overlap.tsv")
require_files(c(count_path, overlap_path))
counts <- fread(count_path, sep = "\t")
overlaps <- fread(overlap_path, sep = "\t")
conditions <- c("WT", "Light", "Dark")
timepoints <- c("8", "24", "72")

ellipse_coordinates <- function(cx, cy, rx, ry, angle, n = 240L) {
  theta <- seq(0, 2 * pi, length.out = n)
  rotation <- angle * pi / 180
  data.table(
    x = cx + rx * cos(theta) * cos(rotation) - ry * sin(theta) * sin(rotation),
    y = cy + rx * cos(theta) * sin(rotation) + ry * sin(theta) * cos(rotation)
  )
}

region_positions <- data.table(
  region = c("WT", "Light", "Dark", "WT&Light", "WT&Dark", "Light&Dark", "WT&Light&Dark"),
  x = c(-1.05, 1.05, 0, 0, -0.55, 0.55, 0),
  y = c(0.38, 0.38, -0.98, 0.55, -0.30, -0.30, 0.02)
)

for (tp in timepoints) {
  count_table <- counts[timepoint == paste0("T", tp)]
  overlap_table <- overlaps[timepoint == paste0("T", tp)]

  count_table[, condition := factor(condition, levels = conditions)]
  bar_plot <- ggplot(count_table, aes(condition, count, fill = direction)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.68, colour = "#4E5BBD") +
    geom_text(
      aes(label = count, colour = direction),
      position = position_dodge(width = 0.78), vjust = -0.35, size = 3.2, show.legend = FALSE
    ) +
    scale_fill_manual(values = c("Up" = "#AD1F24", "Down" = "white")) +
    scale_colour_manual(values = c("Up" = "#AD1F24", "Down" = "#4E5BBD")) +
    labs(x = NULL, y = "Number of DE genes", fill = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

  ellipse_data <- rbindlist(list(
    cbind(ellipse_coordinates(-0.45, 0.20, 1.20, 0.72, 18), set = "WT"),
    cbind(ellipse_coordinates(0.45, 0.20, 1.20, 0.72, -18), set = "Light"),
    cbind(ellipse_coordinates(0, -0.35, 1.10, 0.74, 0), set = "Dark")
  ))
  labels <- merge(overlap_table, region_positions, by = "region", all.x = TRUE)
  venn_plot <- ggplot() +
    geom_polygon(
      data = ellipse_data,
      aes(x, y, group = set, colour = set),
      fill = NA, linewidth = 0.55
    ) +
    geom_text(data = labels, aes(x, y, label = count), size = 3.4) +
    annotate("text", x = -1.15, y = 1.00, label = "WT", size = 3.5) +
    annotate("text", x = 1.15, y = 1.00, label = "Light", size = 3.5) +
    annotate("text", x = 0, y = -1.30, label = "Dark", size = 3.5) +
    scale_colour_manual(values = c("WT" = "#5A5A5A", "Light" = "#3A8CC1", "Dark" = "#9A9A9A")) +
    coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.45, 1.15), clip = "off") +
    labs(title = paste0("DE gene overlap (T", tp, " vs T0, padj < 0.05)")) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

  combined <- bar_plot / venn_plot + plot_layout(heights = c(0.9, 1.1))
  ensure_directory(repo_path("figures"))
  ggsave(
    repo_path("figures", paste0("FigureS2bc_T", tp, "_vs_T0.pdf")),
    combined, width = 6.2, height = 8.0
  )
}

message("Supplementary Figure S2b-c were written from their archived source-data tables.")
