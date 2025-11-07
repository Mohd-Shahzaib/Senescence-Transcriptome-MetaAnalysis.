<p align="center"> <img src="https://img.shields.io/badge/RNA--seq-DESeq2-blue" /> <img src="https://img.shields.io/badge/Batch-SVA%20%2B%20limma-green" /> <img src="https://img.shields.io/badge/Enrichment-GO%20%2B%20Reactome-orange" /> <img src="https://img.shields.io/badge/Reproducible-Seeded%20%2B%20sessionInfo-brightgreen" /> </p> <h1 align="center">Senescence RNA-seq Pipeline</h1> <p align="center">DESeq2 + VST + SVA for visualization + LFC shrinkage + GO and Reactome</p> <p align="center"> <a href="#quick-run">Quick Run</a> • <a href="#highlights">Highlights</a> • <a href="#pipeline">Pipeline</a> • <a href="#what-you-get">What You Get</a> • <a href="#method-guardrails">Method Guardrails</a> • <a href="#reproducibility">Reproducibility</a> </p>
Quick Run
Rscript senescence_rnaseq_pipeline_portfolio.R


Place countmergematrix.csv and metadata_strict_even_batches_ordered.csv next to the script. Outputs appear in the working directory.

Highlights

✅ Correct design: model is ~ Batch + Condition for DE

✅ SVA for display only: use svaseq and removeBatchEffect to show batch reduction in PCA, never for testing

✅ Stable effect sizes: LFC shrinkage uses apeglm if available, else ashr; p values and FDR from non-shrunken results

✅ Strict curation: keep padj < 0.05, |log2FC| ≥ 1, baseMean ≥ 10; drop pseudogenes and LOC*

✅ Defined universe: enrichment uses all tested genes from the DESeq2 object

Pipeline
Counts + Metadata
        │
        ▼
DESeq2 size factors → VST
        │
        ├─ PCA (before)
        ├─ SVA on VST → removeBatchEffect (viz only)
        └─ PCA (after)
        │
        ▼
DE: ~ Batch + Condition
        │
        ├─ LFC shrinkage (apeglm → ashr)
        ├─ Curated DEGs
        ├─ Volcano • Venn • Heatmaps
        └─ GO:BP • Reactome with universe
              └─ Overlap Acute vs Replicative + Fisher combined p

What You Get

Tables

normalized_counts_deseq2.csv, vst_matrix_before_batch.csv, vst_matrix_after_batch.csv, sva_surrogate_variables.csv

DE_full_acute_vs_control_shrunk.csv, DE_full_replicative_vs_control_shrunk.csv

acute_up_degs.csv, acute_down_degs.csv, replicative_up_degs.csv, replicative_down_degs.csv

common_up_degs.csv, common_down_degs.csv

GO: go_bp_*.csv and Reactome: reactome_*.csv

Overlaps with Fisher: *_overlap_*_with_fisher.csv

Figures

PCA: pca_before_batch.png, pca_after_batch.png

Volcano: volcano_acute.jpg, volcano_replicative.jpg

Heatmaps: heatmap_acute_top50_samples.png, heatmap_replicative_top50_samples.png

Senescence panel: heatmap_senescence50_acute.png, heatmap_senescence50_replicative.png

Venns on IDs: go_bp_venn_*_ids.png, reactome_venn_*_ids.png

Top shared terms: *_overlap_*_top20.png

Mini Gallery
<p align="center"> <img src="docs/thumbs/pca_before.png" width="31%" /> <img src="docs/thumbs/pca_after.png" width="31%" /> <img src="docs/thumbs/volcano.png" width="31%" /> </p> <p align="center"> <img src="docs/thumbs/top50_heatmap.png" width="31%" /> <img src="docs/thumbs/senescence50.png" width="31%" /> <img src="docs/thumbs/go_overlap.png" width="31%" /> </p>

Tip: create docs/thumbs/ with 600 px wide copies of your plots so the README stays lightweight.

Method Guardrails

No testing on batch-removed matrices. Batch removal is for display only

LFC shrinkage stabilizes effect sizes; inference uses the non-shrunken Wald results

Universe is the set of tested genes from the DESeq2 object

Filters remove pseudogenes and LOC* symbols and require basic expression and FDR thresholds

Reproducibility

Seed is fixed: set.seed(20250825)

Full environment is saved to sessionInfo.txt

Consider a locked setup:

install.packages("renv")
renv::init()
renv::snapshot()

<details> <summary><b>Inputs</b> (click to expand)</summary>

countmergematrix.csv
Rows are genes. Columns are sample IDs. Raw integer counts.

metadata_strict_even_batches_ordered.csv
Row names are sample IDs. Columns:
Condition ∈ {control, replicative_senescence, acute_senescence}
Batch as a categorical label

</details> <details> <summary><b>Install</b> (click to expand)</summary>
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","limma","sva","AnnotationDbi","org.Hs.eg.db","clusterProfiler","ReactomePA"))
install.packages(c("ggplot2","pheatmap","readr","dplyr","tidyr","ggrepel","ggvenn","stringr"))
BiocManager::install(c("apeglm","ashr"))   # optional shrinkers

</details>
