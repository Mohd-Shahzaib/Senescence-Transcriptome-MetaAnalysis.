# Senescence RNA‑seq Pipeline: DESeq2 + VST + SVA + GO and Reactome

A compact, end‑to‑end R workflow for RNA‑seq normalization, batch visualization with SVA, DEG curation, and downstream enrichment. Outputs are named consistently for a clean portfolio and are written to the working directory by default.

---

## Contents

* [What this does](#what-this-does)
* [Inputs you must provide](#inputs-you-must-provide)
* [Quick start](#quick-start)
* [Required R packages](#required-r-packages)
* [Design and methods](#design-and-methods)
* [Outputs](#outputs)
* [Figures guide](#figures-guide)
* [Reproducibility tips](#reproducibility-tips)
* [Suggested repo layout](#suggested-repo-layout)
* [Cite and acknowledge](#cite-and-acknowledge)

---

## What this does

1. Loads your count matrix and metadata, aligns samples, and prepares a DESeq2 dataset.
2. Normalizes counts and applies VST for PCA and heatmaps.
3. Estimates surrogate variables with `svaseq`, then removes batch signal for visualization using `removeBatchEffect`.
4. Reads precomputed DGE tables for Acute vs Control and Replicative vs Control, filters DEGs using strict rules, and exports clean lists.
5. Generates volcano plots and compact heatmaps for top DEGs and a 50‑gene senescence panel.
6. Runs GO Biological Process and Reactome enrichment for Up and Down sets separately.
7. Computes overlaps between Acute and Replicative enrichments, combines p values with Fisher’s method, adjusts with BH, and saves ranked tables and bar charts.

> Note: SVA and `removeBatchEffect` are used for visualization only. Statistical testing should be performed on the original DESeq2 model space that accounts for design terms.

---

## Inputs you must provide

Place these files in the same folder as the R script.

1. **`countmergematrix.csv`**
   Rows are genes. Columns are sample IDs. Values are raw integer counts.

2. **`metadata_strict_even_batches_ordered.csv`**
   Row names are sample IDs matching the count matrix column names. Required columns:

   * `Condition` with values in `{control, replicative_senescence, acute_senescence}`
   * `Batch` as a categorical label such as `B1`, `B2`, `B3`

3. **`DGE_Acute_vs_Control_Annotated_CORRECTED_EXP.csv`**
   1 row per gene. Required columns: `Symbol`, `Description`, `padj`, `log2FoldChange`, `GeneID` (Entrez). Additional columns are fine.

4. **`DGE_Replicative_vs_Control_Annotated_CORRECTED_EXP.csv`**
   Same columns as in the Acute file.

> The enrichment steps expect `GeneID` to contain human Entrez IDs. If your IDs are different, map to Entrez first.

---

## Quick start

Run the pipeline with base R or from a terminal.

```bash
Rscript senescence_rnaseq_pipeline_portfolio.R
```

All outputs will appear in the working directory. Paths in the script are kept simple on purpose, so you can drop this into any project.

---

## Required R packages

Install Bioconductor packages first, then CRAN packages.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "sva", "org.Hs.eg.db", "ReactomePA", "clusterProfiler"
))

install.packages(c(
  "limma", "ggplot2", "pheatmap", "readr", "dplyr", "ggrepel",
  "VennDiagram", "ggvenn", "stringr", "tidyverse"
))
```

Tested with modern R (4.x). Use the latest stable release for best results.

---

## Design and methods

* **Sample alignment**: samples are intersected between counts and metadata, order is enforced, and consistency is checked.
* **Design**: flexible formula that includes `Batch` and `Condition` if both vary. Otherwise falls back to the available factor.
* **Filtering**: genes with at least 10 counts in at least 2 samples are kept.
* **Normalization and transform**: size factors via DESeq2, then VST for PCA and heatmaps.
* **SVA**: surrogate variables from `svaseq(assay(vsd), mod, mod0)`. Batch signal is removed for PCA only using `removeBatchEffect` with both `batch` and `covariates`.
* **DGE curation**: reads precomputed DGE tables, removes pseudogenes and `LOC*`, keeps `padj < 0.05` and `|log2FC| ≥ 1`. Volcano plots cap extreme fold changes for readability.
* **Overlap of DEGs**: intersection of Symbols between Acute and Replicative for Up and for Down, saved to CSV.
* **GO:BP enrichment**: `clusterProfiler::enrichGO` with `OrgDb = org.Hs.eg.db`, `keyType = ENTREZID`, BH adjusted.
* **Reactome enrichment**: `ReactomePA::enrichPathway` with organism `human` and BH adjusted.
* **Acute vs Replicative overlap on terms**: inner join on `ID` and `Description`, Fisher’s combined p, BH adjusted, then top 20 terms are plotted.

---

## Outputs

**Normalization and PCA**

* `pca_before_batch_correction.png`
* `pca_after_batch_correction.png`
* `normalized_counts_deseq2.csv`
* `vst_matrix_before_batch.csv`
* `vst_matrix_after_batch.csv`
* `sva_surrogate_variables.csv`

**DEG tables**

* `acute_up_degs.csv`
* `acute_down_degs.csv`
* `replicative_up_degs.csv`
* `replicative_down_degs.csv`
* `common_up_degs.csv`
* `common_down_degs.csv`

**Venn diagrams**

* `venn_degs_up.png`
* `venn_degs_down.png`

**Volcano plots**

* `volcano_acute.jpg`
* `volcano_replicative.jpg`

**Heatmaps**

* `heatmap_acute_top50_degs.png`
* `heatmap_replicative_top50_degs.png`
* `heatmap_senescence50_acute.png`
* `heatmap_senescence50_replicative.png`

**GO:BP enrichment**

* `go_bp_acute_up.csv`
* `go_bp_acute_down.csv`
* `go_bp_replicative_up.csv`
* `go_bp_replicative_down.csv`
* `go_bp_venn_up.png`
* `go_bp_venn_down.png`
* `go_bp_overlap_up_with_fisher.csv`
* `go_bp_overlap_down_with_fisher.csv`
* `go_bp_overlap_up_top20.png`
* `go_bp_overlap_down_top20.png`

**Reactome enrichment**

* `reactome_acute_up.csv`
* `reactome_acute_down.csv`
* `reactome_replicative_up.csv`
* `reactome_replicative_down.csv`
* `reactome_venn_up.png`
* `reactome_venn_down.png`
* `reactome_overlap_up_with_fisher.csv`
* `reactome_overlap_down_with_fisher.csv`
* `reactome_overlap_up_top20.png`
* `reactome_overlap_down_top20.png`

---

## Figures guide

* **PCA before and after** shows how much batch and condition structure dominate. After SVA removal of batch for visualization the clusters by condition should sharpen.
* **Volcano** highlights Up, Down, Significant, and NS. Key senescence markers are annotated.
* **Top 50 heatmaps** are split into 25 strongest Up and 25 strongest Down by adjusted p value.
* **Senescence gene set heatmaps** show a fixed 50‑gene panel for quick visual checks.
* **GO and Reactome Venns** reveal term overlap between Acute and Replicative signatures. Bar charts rank the most consistent terms using Fisher‑combined evidence.

---

## Reproducibility tips

* The script sets `set.seed(20250825)`.
* Save the package versions used by exporting `sessionInfo()` to a text file.

```r
sink("sessionInfo.txt"); sessionInfo(); sink()
```

* If you prefer a locked environment, initialize `renv` in the repo and snapshot.

```r
install.packages("renv"); renv::init(); renv::snapshot()
```

---

## Suggested repo layout

The script writes outputs to the working directory. If you want a tidy repo, you can create folders and move files after the run.

```
.
├── data/                      # place input CSVs here (optional)
├── scripts/
│   └── senescence_rnaseq_pipeline_portfolio.R
├── results/
│   ├── tables/
│   └── plots/
└── README.md
```

Example post‑run organization in R:

```r
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(pattern = "\\\.png$|\\\.jpg$"), "results/plots", overwrite = TRUE)
file.copy(list.files(pattern = "\\\.csv$"), "results/tables", overwrite = TRUE)
```

---

## Cite and acknowledge

If you use this repository in a talk or manuscript, please cite the core R packages you rely on in this workflow: DESeq2, limma, sva, clusterProfiler, ReactomePA, org.Hs.eg.db, ggplot2, and pheatmap. Add your own study citation here as well.

Owner: Mohd Shahzaib
