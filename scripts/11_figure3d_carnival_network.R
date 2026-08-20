#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggraph)
  library(ggplot2)
  library(igraph)
  library(readr)
  library(scales)
  library(tibble)
})
source("scripts/00_utils.R")

result_nodes <- repo_path("results", "carnival", "all_model_nodes.csv")
result_edges <- repo_path("results", "carnival", "all_model_edges.csv")
source_nodes <- repo_path("data", "figure_data", "Figure3d_carnival_nodes.csv")
source_edges <- repo_path("data", "figure_data", "Figure3d_carnival_edges.csv")
nodes_path <- if (file.exists(result_nodes)) result_nodes else source_nodes
edges_path <- if (file.exists(result_edges)) result_edges else source_edges
require_files(c(nodes_path, edges_path))
nodes_all <- read_csv(nodes_path, show_col_types = FALSE)
edges_all <- read_csv(edges_path, show_col_types = FALSE)

consensus_edges <- edges_all |>
  filter(weight > 0) |>
  group_by(source, target, interaction) |>
  summarise(
    recurrence = n_distinct(paste(model, timepoint, sep = "__")),
    mean_weight = mean(weight),
    .groups = "drop"
  ) |>
  arrange(desc(recurrence), source, target)

consensus_nodes <- nodes_all |>
  filter(node_id %in% unique(c(consensus_edges$source, consensus_edges$target))) |>
  group_by(node_id) |>
  summarise(
    node_role = case_when(
      any(node_role == "CREB1 anchor") ~ "CREB1 anchor",
      any(node_role == "PRKACA anchor") ~ "PRKACA anchor",
      any(node_role == "Measured motif TF") ~ "Measured motif TF",
      TRUE ~ "Intermediate"
    ),
    mean_activity = mean(inferred_activity, na.rm = TRUE),
    max_abs_activity = max(abs(inferred_activity), na.rm = TRUE),
    .groups = "drop"
  )

vertices <- consensus_nodes |>
  transmute(
    name = node_id,
    NodeType = factor(case_when(
      node_role == "Measured motif TF" ~ "Measurement",
      node_role %in% c("CREB1 anchor", "PRKACA anchor") ~ "Hypothesis anchor",
      TRUE ~ "Intermediate"
    ), levels = c("Intermediate", "Measurement", "Hypothesis anchor")),
    AvgAct = mean_activity,
    NormAbsAct = rescale(max_abs_activity, to = c(0, 1)),
    label = recode(
      node_id,
      "PRKACA" = "PKA\n(PRKACA)",
      "PRKAA1_PRKAA2_PRKAB1_PRKAB2_PRKAG1_PRKAG2_PRKAG3" = "AMPK complex\n(PRKAA/B/G)",
      .default = node_id
    )
  )
graph <- graph_from_data_frame(consensus_edges, directed = TRUE, vertices = vertices)

distance_from <- function(node) {
  if (!(node %in% V(graph)$name)) return(rep(NA_real_, vcount(graph)))
  values <- as.numeric(distances(graph, v = node, to = V(graph), mode = "out")[1, ])
  values[is.infinite(values)] <- NA_real_
  values
}
layers <- tibble(
  node = V(graph)$name,
  perturbation = distance_from("Perturbation"),
  prkaca = distance_from("PRKACA"),
  creb1 = distance_from("CREB1")
) |>
  mutate(layer = case_when(
    node == "Perturbation" ~ 1L,
    node == "PRKACA" ~ 2L,
    node == "CREB1" ~ 3L,
    !is.na(creb1) ~ pmin(8L, 3L + as.integer(creb1)),
    !is.na(prkaca) ~ pmin(8L, 2L + as.integer(prkaca)),
    !is.na(perturbation) ~ pmin(8L, 1L + as.integer(perturbation)),
    TRUE ~ 8L
  ))

set.seed(436)
layout <- layout_with_sugiyama(
  graph,
  layers = layers$layer[match(V(graph)$name, layers$node)],
  hgap = 1.15,
  vgap = 1.10,
  maxiter = 250
)$layout
layout <- tibble(raw_x = layout[, 1], raw_y = layout[, 2], name = V(graph)$name) |>
  mutate(
    x = rescale(max(raw_y) - raw_y, to = c(-2.4, 2.4)),
    y = rescale(raw_x, to = c(-2.05, 2.05)),
    x = case_when(name == "Perturbation" ~ -2.25, name == "PRKACA" ~ -1.25,
                  name == "CREB1" ~ -0.25, TRUE ~ x),
    y = ifelse(name %in% c("Perturbation", "PRKACA", "CREB1"), 0, y)
  )
plot_layout <- create_layout(graph, layout = "manual", x = layout$x, y = layout$y)
activity_limit <- max(abs(vertices$AvgAct), na.rm = TRUE)

plot <- ggraph(plot_layout) +
  geom_edge_link(
    aes(
      edge_colour = factor(interaction),
      edge_linetype = factor(interaction),
      edge_width = recurrence
    ),
    arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"),
    end_cap = circle(2.5, "mm"), alpha = 0.72
  ) +
  geom_node_point(
    aes(fill = AvgAct, shape = NodeType, size = NormAbsAct),
    colour = "#252A32", stroke = 0.65
  ) +
  geom_node_label(
    aes(label = label), repel = TRUE, size = 2.5,
    label.size = 0.16, fill = alpha("white", 0.86), max.overlaps = Inf
  ) +
  scale_edge_colour_manual(values = c("-1" = "#2C7BB6", "1" = "#D73027"),
                           labels = c("-1" = "Inhibition", "1" = "Activation"),
                           name = "Interaction sign") +
  scale_edge_linetype_manual(values = c("-1" = "longdash", "1" = "solid"), guide = "none") +
  scale_edge_width_continuous(range = c(0.4, 2.5), breaks = 1:6,
                              name = "Recurrence across runs") +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
                       limits = c(-activity_limit, activity_limit), name = "AvgAct") +
  scale_shape_manual(values = c("Intermediate" = 21, "Measurement" = 22, "Hypothesis anchor" = 23),
                     name = "Node class") +
  scale_size_continuous(range = c(2.6, 7.2), guide = "none") +
  coord_equal(clip = "off") +
  theme_graph(base_size = 8, base_family = "Helvetica", background = "white") +
  theme(legend.position = "right", plot.margin = margin(8, 12, 8, 8))

output_dir <- ensure_directory(repo_path("results", "carnival"))
write_csv(consensus_nodes, file.path(output_dir, "Figure3d_consensus_nodes.csv"))
write_csv(consensus_edges, file.path(output_dir, "Figure3d_consensus_edges.csv"))
ensure_directory(repo_path("figures"))
ggsave(repo_path("figures", "Figure3d_carnival_network.pdf"), plot, width = 180, height = 160, units = "mm")
ggsave(repo_path("figures", "Figure3d_carnival_network.png"), plot, width = 180, height = 160,
       units = "mm", dpi = 600, bg = "white")
message("Figure 3d was written.")
