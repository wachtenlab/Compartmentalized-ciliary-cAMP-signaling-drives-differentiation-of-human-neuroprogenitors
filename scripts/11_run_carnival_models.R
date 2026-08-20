#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CARNIVAL)
  library(dplyr)
  library(readr)
  library(tidyr)
})
source("scripts/00_utils.R")

measurement_path <- repo_path("data", "carnival", "tf_measurements.csv")
pkn_dir <- repo_path("results", "carnival", "pkn")
require_files(c(
  measurement_path,
  file.path(pkn_dir, paste0("pkn_", c("T8", "T24", "T72"), ".csv"))
))
measurements <- read_csv(measurement_path, show_col_types = FALSE)
output_dir <- ensure_directory(repo_path("results", "carnival", "models"))

runs <- bind_rows(
  expand_grid(model = "CREB1", timepoint = c("T8", "T24", "T72"), perturbation = "CREB1"),
  expand_grid(model = "PRKACA", timepoint = c("T8", "T24", "T72"), perturbation = "PRKACA"),
  expand_grid(model = "inverse", timepoint = c("T8", "T24", "T72"), perturbation = NA_character_)
)

for (i in seq_len(nrow(runs))) {
  run <- runs[i, ]
  pkn <- read_csv(file.path(pkn_dir, paste0("pkn_", run$timepoint, ".csv")), show_col_types = FALSE) |>
    transmute(source, interaction = as.numeric(interaction), target)
  pkn_nodes <- unique(c(pkn$source, pkn$target))
  measurement_table <- measurements |>
    filter(timepoint == run$timepoint, tf %in% pkn_nodes) |>
    select(tf, measurement)
  measurement_vector <- setNames(measurement_table$measurement, measurement_table$tf)

  model_dir <- ensure_directory(file.path(output_dir, run$model))
  prefix <- paste0(run$timepoint, "_", run$model)
  options <- defaultLpSolveCarnivalOptions()
  options$outputFolder <- ensure_directory(file.path(model_dir, paste0(prefix, "_solver")))

  result <- if (run$model == "inverse") {
    runInverseCarnival(
      measurements = measurement_vector,
      priorKnowledgeNetwork = pkn,
      carnivalOptions = options
    )
  } else {
    runVanillaCarnival(
      perturbations = setNames(1, run$perturbation),
      measurements = measurement_vector,
      priorKnowledgeNetwork = pkn,
      carnivalOptions = options
    )
  }

  saveRDS(result, file.path(model_dir, paste0(prefix, ".rds")))
  write_csv(as_tibble(result$weightedSIF), file.path(model_dir, paste0(prefix, "_edges.csv")))
  write_csv(as_tibble(result$nodesAttributes), file.path(model_dir, paste0(prefix, "_nodes.csv")))
}

write_csv(runs, file.path(output_dir, "completed_runs.csv"))
message("Completed ", nrow(runs), " CARNIVAL models.")

