#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})
source("scripts/00_utils.R")

model_dir <- repo_path("results", "carnival", "models")
measurement_path <- repo_path("data", "carnival", "tf_measurements.csv")
require_files(measurement_path)
node_files <- list.files(model_dir, pattern = "_nodes\\.csv$", recursive = TRUE, full.names = TRUE)
edge_files <- list.files(model_dir, pattern = "_edges\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(node_files) != 9L || length(edge_files) != 9L) {
  stop("Expected node and edge tables for nine completed CARNIVAL models.")
}
measurements <- read_csv(measurement_path, show_col_types = FALSE)

file_metadata <- function(path, suffix) {
  stem <- str_remove(basename(path), suffix)
  tibble(
    path = path,
    timepoint = str_extract(stem, "^T[0-9]+"),
    model = str_remove(stem, "^T[0-9]+_")
  )
}

nodes <- map_dfr(node_files, function(path) {
  metadata <- file_metadata(path, "_nodes\\.csv$")
  read_csv(path, show_col_types = FALSE) |>
    transmute(
      node_id = Node,
      inferred_activity = AvgAct,
      timepoint = metadata$timepoint,
      model = metadata$model
    )
}) |>
  left_join(
    measurements |> select(timepoint, tf, motif_measurement = measurement),
    by = c("timepoint", "node_id" = "tf")
  ) |>
  mutate(
    node_role = case_when(
      node_id == "CREB1" ~ "CREB1 anchor",
      node_id == "PRKACA" ~ "PRKACA anchor",
      !is.na(motif_measurement) ~ "Measured motif TF",
      TRUE ~ "Intermediate"
    )
  )

edges <- map_dfr(edge_files, function(path) {
  metadata <- file_metadata(path, "_edges\\.csv$")
  read_csv(path, show_col_types = FALSE) |>
    transmute(
      source = Node1,
      interaction = as.numeric(Sign),
      target = Node2,
      weight = as.numeric(Weight),
      timepoint = metadata$timepoint,
      model = metadata$model
    )
})

output_dir <- ensure_directory(repo_path("results", "carnival"))
write_csv(nodes, file.path(output_dir, "all_model_nodes.csv"))
write_csv(edges, file.path(output_dir, "all_model_edges.csv"))
message("Exported combined CARNIVAL node and edge tables.")

