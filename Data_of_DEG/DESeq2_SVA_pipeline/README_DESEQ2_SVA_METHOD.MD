Senescence RNA-seq Pipeline: DESeq2 + VST + SVA + GO and Reactome

A compact end to end R workflow for RNA-seq normalization, PCA with batch visualization, differential expression with LFC shrinkage, curated DEG tables, and GO/Reactome enrichment. All outputs are written to the working directory with consistent names.

Contents

What this does

Inputs you must provide

Quick start

Required R packages

Design and methods

Outputs

Figures guide

Reproducibility

Suggested repo layout

Cite and acknowledge

What this does

Loads the count matrix and metadata, aligns samples, and builds a DESeq2 dataset

Normalizes counts and applies VST for PCA and heatmaps

Estimates surrogate variables with svaseq, then removes batch signal for visualization only with removeBatchEffect

Fits ~ Batch + Condition and produces DE with LFC shrinkage (apeglm preferred, ashr fallback)

Curates DEGs with strict filters and removes pseudogenes and LOC symbols

Generates volcano plots, Venn diagrams, sample-level heatmaps, and a 50-gene senescence panel

Runs GO:BP and Reactome enrichment on Entrez IDs with a tested-gene universe

Computes Acute vs Replicative overlaps for enriched terms and combines p values with Fisher’s method

Note: SVA and removeBatchEffect are used for visualization only. DE is run on the model space with Batch + Condition.

Inputs you must provide

Put these files next to the R script.

countmergematrix.csv
Rows are genes (HGNC symbols). Columns are sample IDs. Values are raw integer counts.

metadata_strict_even_batches_ordered.csv
Row names are sample IDs that match the count matrix column names. Required columns:

Condition with values in {control, replicative_senescence, acute_senescence}

Batch as a categorical label

No precomputed DE tables are required. The script performs DE and creates all downstream files.

Quick start

Run from R or a terminal in the folder that contains the inputs.

Rscript senescence_rnaseq_pipeline_portfolio.R


All outputs will appear in the working directory.

Required R packages

Install Bioconductor packages first, then CRAN packages.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "limma", "sva",
  "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "ReactomePA"
))

install.packages(c(
  "ggplot2", "pheatmap", "readr", "dplyr",
  "tidyr", "ggrepel", "ggvenn", "stringr"
))

# Optional shrinkers used by lfcShrink
BiocManager::install(c("apeglm","ashr"))


The script auto-detects apeglm. If it is not installed it falls back to ashr. If you want both to appear in sessionInfo(), add library(apeglm) and library(ashr) near the top.

Design and methods

Sample alignment. Intersect samples between counts and metadata. Enforce identical order.

Design. ~ Batch + Condition (both factors are used).

Filtering. Keep genes with counts ≥ 10 in at least two samples.

Normalization and transform. DESeq2 median of ratios. Then VST for PCA and heatmaps.

SVA for visualization. Run svaseq on the VST matrix. Use removeBatchEffect with Batch and SVs for PCA and heatmaps only.

Differential expression. Two contrasts: acute_senescence vs control and replicative_senescence vs control. LFC shrinkage uses apeglm if available else ashr. Wald p values and FDR are kept from the non-shrunken results.

Curated DEG filters. padj < 0.05, |log2FC| ≥ 1, baseMean ≥ 10. Remove rows with missing values, pseudogenes in Description, blank symbols, and symbols that start with LOC.

Visualization. Volcano plots with a small set of senescence markers labeled. Venn diagrams of DEG overlaps. Sample-level heatmaps for top 25 up and 25 down. Fixed 50-gene senescence panel heatmaps.

Enrichment. Map Symbols to Entrez. Use the tested-gene universe from the DESeq2 object. Run GO:BP and Reactome with BH FDR < 0.05.

Term overlap and meta p. For terms present in both contrasts, combine per-contrast p values with Fisher’s method and adjust again with BH.

Outputs

Normalization and PCA

normalized_counts_deseq2.csv

vst_matrix_before_batch.csv

vst_matrix_after_batch.csv

pca_before_batch.png

pca_after_batch.png

sva_surrogate_variables.csv

Differential expression and curated DEGs

DE_full_acute_vs_control_shrunk.csv

DE_full_replicative_vs_control_shrunk.csv

acute_up_degs.csv

acute_down_degs.csv

replicative_up_degs.csv

replicative_down_degs.csv

common_up_degs.csv

common_down_degs.csv

Plots

volcano_acute.jpg

volcano_replicative.jpg

venn_degs_up.png

venn_degs_down.png

heatmap_acute_top50_samples.png

heatmap_replicative_top50_samples.png

heatmap_senescence50_acute.png

heatmap_senescence50_replicative.png

Enrichment and overlaps

GO: go_bp_acute_up.csv, go_bp_acute_down.csv, go_bp_replicative_up.csv, go_bp_replicative_down.csv

Reactome: reactome_acute_up.csv, reactome_acute_down.csv, reactome_replicative_up.csv, reactome_replicative_down.csv

Fisher overlap tables:

go_bp_overlap_up_with_fisher.csv, go_bp_overlap_down_with_fisher.csv

reactome_overlap_up_with_fisher.csv, reactome_overlap_down_with_fisher.csv

Venns on IDs:

go_bp_venn_up_ids.png, go_bp_venn_down_ids.png

reactome_venn_up_ids.png, reactome_venn_down_ids.png

Top 20 overlap bar charts:

go_bp_overlap_up_top20.png, go_bp_overlap_down_top20.png

reactome_overlap_up_top20.png, reactome_overlap_down_top20.png

Figures guide

PCA before vs after shows how SVA removal of batch for display sharpens the condition pattern.

Volcano shows Up, Down, Significant, and NS regions with classical markers annotated.

Top 50 heatmaps use the 25 most significant up and 25 most significant down genes per contrast.

Senescence marker heatmaps show a fixed panel of 50 genes for a quick check.

GO and Reactome Venns show term overlaps between Acute and Replicative.

Top overlap charts rank shared terms by Fisher-combined evidence.

Reproducibility

The script sets set.seed(20250825) before SVA and at the top of the file.

The full environment is saved to sessionInfo.txt.

For a frozen setup, use renv if you prefer:

install.packages("renv")
renv::init()
renv::snapshot()

Suggested repo layout
.
├── data/                      # optional, put inputs here
├── scripts/
│   └── senescence_rnaseq_pipeline_portfolio.R
├── results/
│   ├── tables/
│   └── plots/
└── README.md


You can move outputs after the run:

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(pattern = "\\.png$|\\.jpg$"), "results/plots", overwrite = TRUE)
file.copy(list.files(pattern = "\\.csv$"), "results/tables", overwrite = TRUE)

Cite and acknowledge

Please cite DESeq2, limma, sva, clusterProfiler, ReactomePA, org.Hs.eg.db, ggplot2, and pheatmap. Include your own study citation. Add a license to the repo.

Owner: Mohd Shahzaib
