# RNA-seq analysis code for Frechen et al.

This repository contains the reproducible RNA-seq analysis used for the manuscript
"Compartmentalized ciliary cAMP signaling drives differentiation of human neuroprogenitors".
The workflow starts from Salmon `quant.sf` files, imports transcript abundance with
`tximport`, fits the DESeq2 interaction model, exports contrasts, and recreates the
RNA-seq-derived manuscript panels.

## Scope

Included manuscript panels:

- Figure 3b: neurogenic driver expression.
- Figure 3c: Transcription factor motif enrichment.
- Figure 3d: time-collapsed CARNIVAL consensus network.
- Supplementary Figure S2a: PCA.
- Supplementary Figure S2b: numbers of up- and downregulated genes.
- Supplementary Figure S2c: overlap of differentially expressed genes.
- Supplementary Figure S2d: selected GO biological-process summaries.
- Supplementary Figure S2e: neural progenitor marker expression.
- Supplementary Figure S2f: upstream candidate scan.

Not included:

- RNA-independent microscopy, biosensor, qPCR, and tissue-analysis panels.

See `FIGURE_MANIFEST.md` for the file-level mapping between scripts, inputs, and
manuscript panels.

## Repository layout

```text
data/
  annotations/   Curated gene and GO-term lists used by plotted panels
  carnival/      Signed TF measurements and the frozen OmniPath network
  counts/        Processed gene-level count matrix
  figure_data/   Plot-level source data for every included panel
  metadata/      Sample metadata
scripts/         Numbered analysis scripts
```


Create the environment with:

```bash
mamba env create -f environment.yml
mamba activate rnaseq-analysis
```

## Input data

The repository contains both the processed gene-level count matrix and the final
tximport object used for the DESeq2 model. To restart from Salmon output, place the
following files at the paths shown below:

```text
data/salmon_quant/<sample>/quant.sf
data/reference/gencode.v49.annotation.gtf.gz
```

The sample directory names must match `sample` in
`data/metadata/sample_metadata.tsv`. GENCODE release 49 was used.

Raw FASTQ files and Salmon quantification directories are not duplicated in this
code repository.

## Run order

Run scripts from the repository root:

```bash
Rscript scripts/01_tximport.R
Rscript scripts/02_deseq2.R
Rscript scripts/03_export_contrasts.R
Rscript scripts/04_RNASeq_analysis_preprocessing.ipynb
Rscript scripts/05_figure3b_neurogenic_markers.R
Rscript scripts/06_figureS2bc_deg_counts_overlap.R
Rscript scripts/07_figureS2d_go_enrichment.R
Rscript scripts/08_figureS2e_npc_markers.R
Rscript scripts/09_tf_motif_enrichment_analysis.Rmd
Rscript scripts/10_prepare_carnival_networks.R
Rscript scripts/11_run_carnival_models.R
Rscript scripts/12_export_carnival_results.R
Rscript scripts/13_figure3d_carnival_network.R
Rscript scripts/14_figureS2f_upstream_scan.R
```

`01_tximport.R` is optional when using the processed data included in this repository.
`02_deseq2.R` uses `data/counts/tximport_gene_level.rds` by default and otherwise
uses `data/counts/gene_counts.csv`.


## Statistical analysis

Gene-level abundance was summarized with `countsFromAbundance =
"lengthScaledTPM"`. DESeq2 used the design:

```r
~ timepoint + condition + timepoint:condition
```

WT and T0 were the reference levels. If a metadata table contains more than one
batch, `batch` is added as an additive covariate. Differentially expressed genes
were defined by Benjamini-Hochberg adjusted p-value below 0.05.

## Reproducibility notes

- The full signed, directed OmniPath network used in the analysis is frozen in
  `data/carnival/omnipath_signed_directed.csv`; the exact focused model network is
  retained in `data/carnival/omnipath_focused.csv`.
- Plot-level source tables are retained in `data/figure_data`, so every included
  panel can be inspected without rerunning model fitting.
