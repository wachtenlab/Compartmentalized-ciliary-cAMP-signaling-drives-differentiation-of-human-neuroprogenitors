#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
  library(tibble)
})
source("scripts/00_utils.R")

pkn_path <- repo_path("data", "carnival", "omnipath_signed_directed.csv")
measurement_path <- repo_path("data", "carnival", "tf_measurements.csv")
require_files(c(pkn_path, measurement_path))
maximum_steps <- 4L
timepoints <- c("T8", "T24", "T72")

pkn <- read_csv(pkn_path, show_col_types = FALSE) |>
  transmute(source, interaction = as.integer(interaction), target) |>
  filter(interaction %in% c(-1L, 1L), source != target) |>
  distinct()
measurements <- read_csv(measurement_path, show_col_types = FALSE) |>
  filter(timepoint %in% timepoints) |>
  mutate(target_sign = if_else(measurement > 0, 1L, -1L))
all_nodes <- sort(unique(c(pkn$source, pkn$target)))
source_nodes <- sort(unique(pkn$source))
measurements <- measurements |> filter(tf %in% all_nodes)

state_name <- function(node, sign) paste0(node, ifelse(sign > 0, "::POS", "::NEG"))
signed_edges <- bind_rows(
  pkn |> transmute(from = state_name(source, 1L), to = state_name(target, interaction)),
  pkn |> transmute(from = state_name(source, -1L), to = state_name(target, -interaction))
) |> distinct()
signed_graph <- graph_from_data_frame(
  signed_edges,
  directed = TRUE,
  vertices = tibble(name = c(state_name(all_nodes, 1L), state_name(all_nodes, -1L)))
)
candidate_states <- as.vector(rbind(state_name(source_nodes, 1L), state_name(source_nodes, -1L)))

distance_rows <- vector("list", nrow(measurements))
for (i in seq_len(nrow(measurements))) {
  row <- measurements[i, ]
  distance <- distances(
    signed_graph,
    v = c(state_name(row$tf, row$target_sign), state_name(row$tf, -row$target_sign)),
    to = candidate_states,
    mode = "in",
    weights = NA
  )
  distance_rows[[i]] <- tibble(
    timepoint = row$timepoint,
    candidate = rep(source_nodes, each = 2L),
    orientation = rep(c(1L, -1L), times = length(source_nodes)),
    desired_distance = as.numeric(distance[1, ]),
    opposite_distance = as.numeric(distance[2, ])
  )
}

distances_table <- bind_rows(distance_rows) |>
  mutate(
    desired = is.finite(desired_distance) & desired_distance >= 1 & desired_distance <= maximum_steps,
    opposite = is.finite(opposite_distance) & opposite_distance >= 1 & opposite_distance <= maximum_steps,
    relation = case_when(
      desired & !opposite ~ "support",
      desired & opposite ~ "ambiguous",
      !desired & opposite ~ "contradiction",
      TRUE ~ "unsupported"
    )
  )

ranking <- distances_table |>
  group_by(timepoint, candidate, orientation) |>
  summarise(
    n_targets = n(),
    signed_support = sum(relation == "support"),
    sign_ambiguous = sum(relation == "ambiguous"),
    signed_contradiction = sum(relation == "contradiction"),
    mean_support_distance = ifelse(any(desired), mean(desired_distance[desired]), Inf),
    .groups = "drop"
  ) |>
  mutate(signed_score = (signed_support + 0.5 * sign_ambiguous - signed_contradiction) / n_targets) |>
  group_by(timepoint, candidate) |>
  arrange(desc(signed_score), desc(signed_support), mean_support_distance, desc(orientation), .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  filter(signed_support + sign_ambiguous + signed_contradiction > 0) |>
  group_by(timepoint) |>
  arrange(desc(signed_score), desc(signed_support), sign_ambiguous, mean_support_distance, candidate,
          .by_group = TRUE) |>
  mutate(rank = row_number(), rank_percentile = 1 - (rank - 1) / pmax(n() - 1, 1)) |>
  ungroup()

overall <- ranking |>
  group_by(candidate) |>
  summarise(
    timepoints_ranked = n_distinct(timepoint),
    mean_rank_percentile = mean(rank_percentile),
    mean_signed_score = mean(signed_score),
    total_signed_support = sum(signed_support),
    median_rank = median(rank),
    .groups = "drop"
  ) |>
  arrange(desc(timepoints_ranked), desc(mean_rank_percentile), desc(mean_signed_score),
          desc(total_signed_support), median_rank, candidate) |>
  mutate(overall_rank = row_number())

comparators <- c(
  "CREB1", "PRKACA", "AKT1", "GSK3B", "MAPK1", "MAPK3", "MAPK14",
  "SRC", "CDK1", "CDK2", "CSNK2A1", "CSNK2A1_CSNK2B"
)
leading <- head(overall$candidate, 12L)
display <- unique(c(leading, comparators))
display <- display[display %in% ranking$candidate]
if (length(display) > 20L) {
  comparator_rows <- comparators[comparators %in% display]
  leading_rows <- overall$candidate[overall$candidate %in% display]
  display <- unique(c(head(leading_rows, 20L - length(comparator_rows)), comparator_rows))
}
display_order <- overall |> filter(candidate %in% display) |> arrange(overall_rank) |> pull(candidate)

plot_data <- ranking |>
  filter(candidate %in% display) |>
  complete(
    candidate = display,
    timepoint = timepoints,
    fill = list(signed_support = 0L, sign_ambiguous = 0L, signed_score = NA_real_, orientation = NA_integer_)
  ) |>
  mutate(
    candidate = factor(candidate, levels = rev(display_order)),
    timepoint = factor(timepoint, levels = timepoints),
    supported_weight = signed_support + 0.5 * sign_ambiguous,
    orientation_label = factor(orientation, levels = c(1L, -1L),
                               labels = c("Candidate activated", "Candidate inhibited"))
  )

source_data <- plot_data |>
  mutate(candidate = as.character(candidate), timepoint = as.character(timepoint))
write_csv(source_data, repo_path("data", "figure_data", "FigureS2f_upstream_candidates.csv"))

plot <- ggplot(plot_data, aes(timepoint, candidate)) +
  geom_point(
    aes(size = supported_weight, fill = signed_score, shape = orientation_label),
    colour = "#20252B", stroke = 0.5, alpha = 0.94
  ) +
  scale_shape_manual(values = c("Candidate activated" = 24, "Candidate inhibited" = 25), na.translate = FALSE) +
  scale_fill_gradient(low = "#DCECF2", high = "#176B87", name = "Signed score") +
  scale_size_continuous(range = c(2.0, 7.2), breaks = c(2, 5, 10, 15), name = "Supported TFs") +
  labs(x = NULL, y = NULL, shape = "Best orientation") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_line(colour = "#E5E7E9", linewidth = 0.3),
        panel.grid.minor = element_blank(), legend.position = "right")

ensure_directory(repo_path("figures"))
ggsave(repo_path("figures", "FigureS2f_upstream_candidates.pdf"), plot, width = 9.2, height = 7.2)
ggsave(repo_path("figures", "FigureS2f_upstream_candidates.png"), plot, width = 9.2, height = 7.2,
       dpi = 400, bg = "white")
message("Supplementary Figure S2f and its source-data table were written.")
