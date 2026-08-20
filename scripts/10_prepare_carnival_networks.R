#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(purrr)
  library(readr)
  library(tibble)
})
source("scripts/00_utils.R")

pkn_path <- repo_path("data", "carnival", "omnipath_focused.csv")
measurement_path <- repo_path("data", "carnival", "tf_measurements.csv")
require_files(c(pkn_path, measurement_path))

pkn <- read_csv(pkn_path, show_col_types = FALSE) |>
  transmute(source, interaction = as.integer(interaction), target) |>
  filter(interaction %in% c(-1L, 1L), source != target) |>
  distinct()
measurements <- read_csv(measurement_path, show_col_types = FALSE)
anchors <- c("CREB1", "PRKACA")
timepoints <- c("T8", "T24", "T72")
maximum_anchor_steps <- 4L
maximum_relay_steps <- 2L
output_dir <- ensure_directory(repo_path("results", "carnival", "pkn"))

focused_pkn <- pkn
write_csv(focused_pkn, file.path(output_dir, "pkn_focused.csv"))

focused_graph <- graph_from_data_frame(
  focused_pkn |> select(source, target, interaction),
  directed = TRUE
)

shortest_edge_path <- function(from, to, maximum_steps) {
  if (!all(c(from, to) %in% V(focused_graph)$name)) return(tibble())
  distance <- suppressWarnings(distances(focused_graph, v = from, to = to, mode = "out")[1, 1])
  if (!is.finite(distance) || distance < 1 || distance > maximum_steps) return(tibble())
  vertices <- V(focused_graph)$name[
    shortest_paths(focused_graph, from = from, to = to, mode = "out")$vpath[[1]]
  ]
  tibble(source = head(vertices, -1), target = tail(vertices, -1)) |>
    left_join(focused_pkn, by = c("source", "target")) |>
    select(source, interaction, target) |>
    distinct()
}

coverage_rows <- list()

for (timepoint in timepoints) {
  observed <- measurements |>
    filter(.data$timepoint == .env$timepoint) |>
    arrange(desc(abs(measurement)))
  present_tfs <- intersect(observed$tf, V(focused_graph)$name)
  paths <- list()

  for (anchor in anchors) {
    for (tf in observed$tf) {
      path <- shortest_edge_path(anchor, tf, maximum_anchor_steps)
      if (nrow(path) > 0L) paths[[length(paths) + 1L]] <- path
    }
  }
  for (source_tf in present_tfs) {
    for (target_tf in setdiff(present_tfs, source_tf)) {
      path <- shortest_edge_path(source_tf, target_tf, maximum_relay_steps)
      if (nrow(path) > 0L) paths[[length(paths) + 1L]] <- path
    }
  }

  if (length(paths) == 0L) stop("No model paths were found for ", timepoint, ".")
  timepoint_pkn <- bind_rows(paths) |>
    distinct(source, interaction, target) |>
    filter(source != target)
  nodes <- unique(c(timepoint_pkn$source, timepoint_pkn$target))
  coverage <- observed |>
    transmute(
      timepoint,
      tf,
      measurement,
      present_in_pkn = tf %in% nodes,
      reachable_from_CREB1 = map_lgl(tf, ~ nrow(shortest_edge_path("CREB1", .x, maximum_anchor_steps)) > 0L),
      reachable_from_PRKACA = map_lgl(tf, ~ nrow(shortest_edge_path("PRKACA", .x, maximum_anchor_steps)) > 0L)
    )
  coverage_rows[[timepoint]] <- coverage
  write_csv(timepoint_pkn, file.path(output_dir, paste0("pkn_", timepoint, ".csv")))
}

write_csv(bind_rows(coverage_rows), file.path(output_dir, "measurement_coverage.csv"))
message("Prepared CARNIVAL networks for ", paste(timepoints, collapse = ", "), ".")
