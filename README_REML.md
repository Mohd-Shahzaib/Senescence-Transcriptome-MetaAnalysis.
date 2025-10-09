# REML Meta-Analysis Pipeline for Senescence Transcriptomics

### Author
**Mohd Shahzaib**  
PhD Candidate, University of Campania “Luigi Vanvitelli”  
© 2025 Mohd Shahzaib. All rights reserved.

---

## Overview

This branch hosts the standalone implementation of **Approach 1: REML-based meta-analysis** from the dual validation framework for transcriptomic integration of **cellular senescence datasets**.

It performs a **random-effects model (REML)** meta-analysis across multiple senescence datasets, combining differential expression results from `DESeq2` into a single unified transcriptional signature.

---

## Workflow Summary

1. **Input Files**
   - `countmergematrix.csv` — merged gene expression matrix  
   - `metadata.csv` — sample annotation with Batch and Condition  
   - `Human_gene_annotation.csv` — GeneID–Symbol mapping file  

2. **Steps**
   - DESeq2 per-dataset normalization  
   - REML meta-analysis using `metafor::rma()`  
   - I² heterogeneity computation  
   - Pseudogene filtering  
   - GO:BP and Reactome enrichment  
   - Volcano + Heatmap visualization  

3. **Outputs**
   - `Meta_DEG_Senescence_Combined_Metafor.csv`  
   - `Significant_Upregulated_logFC1.csv`, `Significant_Downregulated_logFC1.csv`  
   - `GO_BP_*` and `Reactome_*` enrichment tables  
   - `Volcano_Meta_DEGs.png`, `Heatmap_Top50.jpg`  
   - `session_info.txt`

---

## Dependencies

R version ≥ 4.0  
Required packages:  
`DESeq2`, `dplyr`, `metafor`, `clusterProfiler`, `ReactomePA`, `EnhancedVolcano`, `pheatmap`, `ggplot2`, `org.Hs.eg.db`, `patchwork`, `stringr`.

To install:
```R
install.packages(c("dplyr","tibble","readr","ggplot2","stringr"))
BiocManager::install(c("DESeq2","metafor","clusterProfiler","org.Hs.eg.db","ReactomePA","EnhancedVolcano","pheatmap","patchwork"))
