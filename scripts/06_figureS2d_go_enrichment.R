#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(data.table)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(patchwork)
})
source("scripts/00_utils.R")

term_path <- repo_path("data", "annotations", "figureS2d_go_terms.tsv")
require_files(term_path)
terms <- fread(term_path, sep = "\t")
timepoints <- c("T8", "T24", "T72")
alpha <- 0.05

symbol_to_go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  columns = c("GOALL", "ONTOLOGYALL"),
  keytype = "SYMBOL"
)
symbol_to_go <- unique(as.data.table(symbol_to_go)[
  ONTOLOGYALL == "BP" & GOALL %in% terms$go_id,
  .(gene = toupper(SYMBOL), go_id = GOALL)
])

enrichment_test <- function(target, universe, members) {
  target <- intersect(unique(target), universe)
  background <- setdiff(universe, target)
  table <- matrix(
    c(
      sum(target %in% members), sum(!(target %in% members)),
      sum(background %in% members), sum(!(background %in% members))
    ),
    nrow = 2L,
    byrow = TRUE
  )
  test <- fisher.test(table, alternative = "greater")
  list(hit_count = table[1, 1], odds_ratio = unname(test$estimate), pvalue = test$p.value)
}

results <- read_light_vs_dark(timepoints)
summary_rows <- list()
enrichment_rows <- list()

for (tp in timepoints) {
  table <- unique(results[timepoint == tp & !is.na(pvalue)], by = "gene")
  universe <- table$gene
  significant <- table[!is.na(padj) & padj < alpha]
  significant[, direction := ifelse(log2FC_light_vs_dark < 0, "Light-down", "Light-up")]

  for (i in seq_len(nrow(terms))) {
    go_value <- terms$go_id[[i]]
    members <- symbol_to_go[go_id == go_value, gene]
    hits <- significant[gene %in% members]
    tests <- list(
      "All DEGs" = significant$gene,
      "Light-down" = significant[direction == "Light-down", gene],
      "Light-up" = significant[direction == "Light-up", gene]
    )
    enrichment_rows[[length(enrichment_rows) + 1L]] <- rbindlist(lapply(names(tests), function(label) {
      x <- enrichment_test(tests[[label]], universe, members)
      data.table(
        timepoint = tp, go_id = go_value, tested_set = label,
        hit_count = x$hit_count, odds_ratio = x$odds_ratio, pvalue = x$pvalue
      )
    }))
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      timepoint = tp,
      go_id = go_value,
      n_significant_degs = nrow(significant),
      total_hit_count = nrow(hits),
      total_percentage_of_degs = 100 * nrow(hits) / nrow(significant),
      genes = paste(sort(hits$gene), collapse = ";")
    )
  }
}

enrichment <- rbindlist(enrichment_rows)
enrichment[, fdr := p.adjust(pvalue, method = "BH"), by = .(timepoint, tested_set)]
summary_table <- merge(rbindlist(summary_rows), terms, by = "go_id")
minimum_fdr <- enrichment[, .(minimum_fdr = min(fdr)), by = .(timepoint, go_id)]
summary_table <- merge(summary_table, minimum_fdr, by = c("timepoint", "go_id"))

fwrite(
  summary_table[order(timepoint, total_percentage_of_degs)],
  repo_path("data", "figure_data", "FigureS2d_go_term_summary.csv")
)
fwrite(enrichment, repo_path("results", "FigureS2d_go_enrichment.csv"))

maximum <- ceiling((max(summary_table$total_percentage_of_degs) + 0.4) * 2) / 2
make_panel <- function(tp) {
  panel <- summary_table[timepoint == tp][order(total_percentage_of_degs, term)]
  panel[, term := factor(term, levels = panel[order(-total_percentage_of_degs), as.character(term)])]
  ggplot(panel, aes(total_percentage_of_degs, term)) +
    geom_col(fill = "#F2AB2F", width = 0.78) +
    scale_x_continuous(limits = c(0, maximum), expand = expansion(mult = c(0, 0))) +
    labs(title = paste0("Light vs Dark (", tp, ")"), x = "Significant DEGs annotated to term (%)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
}

plot <- wrap_plots(lapply(timepoints, make_panel), nrow = 1)
ensure_directory(repo_path("figures"))
ggsave(repo_path("figures", "FigureS2d_go_terms.pdf"), plot, width = 15.5, height = 5.4)
ggsave(repo_path("figures", "FigureS2d_go_terms.png"), plot, width = 15.5, height = 5.4, dpi = 400)
message("Supplementary Figure S2d and its source-data table were written.")
