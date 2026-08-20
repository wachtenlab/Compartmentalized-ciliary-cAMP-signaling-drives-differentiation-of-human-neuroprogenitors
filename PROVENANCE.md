# Data provenance

The processed count matrix and sample metadata in this repository reproduce the final
75-sample DESeq2 contrast set used for Figure 3b and the downstream marker and GO
analyses.

The bar and overlap values displayed in Supplementary Figure S2b-c were generated
earlier from the complete experiment and retained as archived panel source data. The
exact displayed values are therefore deposited directly in:

- `data/figure_data/FigureS2b_deg_counts.tsv`
- `data/figure_data/FigureS2c_deg_overlap.tsv`

`scripts/05_figureS2bc_deg_counts_overlap.R` recreates the displayed panels from those
tables. This separation prevents a later refit of the DESeq2 object from silently
changing numbers already shown in the manuscript.

The CARNIVAL input measurements originate from the separately performed TF-motif
analysis. Only the three fields used by the network model (`timepoint`, `tf`, and
`measurement`) are retained in this repository.
