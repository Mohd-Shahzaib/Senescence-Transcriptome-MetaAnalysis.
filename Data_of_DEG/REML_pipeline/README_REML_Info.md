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


########################################################################
# pipeline_meta_senescence.R
# Complete meta-analysis + overlap integration pipeline
# Author: Mohd Shahzaib
# Version: Public Release (2025)
#
# Description:
# This pipeline integrates transcriptomic data from multiple senescence models
# using two complementary strategies:
#   (1) REML-based meta-analysis (Approach 1)
#   (2) SVA-corrected DESeq2 overlap (Approach 2)
#
# Both workflows have been unified in this reproducible R pipeline.
# The code is identical in logic and parameterization to that used in the PhD thesis.
#

#   - No manual file re-inputs between stages
#   - Outputs include all DEG tables, GO/Reactome enrichment, and QC plots
########################################################################
#   (1) Following is the approach 1 REML-based meta-analysis (Approach 1) in R. 
# ============================================================
# SECTION 1: INITIAL SETUP AND CONFIGURATION
# ============================================================
counts_file <- "countmergematrix.csv"                  # Merged raw counts (genes x samples)
meta_file   <- "metadata.csv"  # Metadata with Condition & Batch
annot_file  <- "Human_gene_annotation.csv"             # Annotation: GeneID, Symbol, Description, EnsemblGeneID
outdir      <- "."                                     # Output directory (default: current)
print_checks <- TRUE                                   # If TRUE, prints key checkpoints

# ============================================================
# SECTION 2: LOAD REQUIRED LIBRARIES
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(metafor)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(EnhancedVolcano)
  library(pheatmap)
  library(patchwork)
  library(stringr)
})

# ============================================================
# SECTION 3: INPUT VALIDATION
# ============================================================
if (!file.exists(counts_file)) stop("counts_file not found: ", counts_file)
if (!file.exists(meta_file))   stop("meta_file not found: ", meta_file)
if (!file.exists(annot_file))  stop("annot_file not found: ", annot_file)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ============================================================
# SECTION 4: DATA LOADING AND PREPROCESSING
# ============================================================
raw_counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta       <- read.csv(meta_file, row.names = 1, check.names = FALSE)
annot      <- read.csv(annot_file, check.names = FALSE)

# Align samples between count and metadata
samples_to_keep <- intersect(colnames(raw_counts), rownames(meta))
raw_counts <- raw_counts[, samples_to_keep, drop = FALSE]
meta       <- meta[samples_to_keep, , drop = FALSE]
stopifnot(identical(colnames(raw_counts), rownames(meta)))

# Annotation subset (clean unique GeneIDs)
annot_sub <- annot %>%
  select(GeneID, Symbol, Description, EnsemblGeneID) %>%
  distinct(GeneID, .keep_all = TRUE)
annot_sub$GeneID <- as.character(annot_sub$GeneID)

# ============================================================
# SECTION 5: DIFFERENTIAL EXPRESSION ANALYSIS (PER DATASET)
# ============================================================
run_dge_per_dataset <- function(batch_name, raw_counts, meta) {
  cat("Processing batch:", batch_name, "\n")
  sub_samples <- rownames(meta)[meta$Batch == batch_name]
  if (length(sub_samples) < 2) return(NULL)

  sub_counts <- raw_counts[, sub_samples, drop = FALSE]
  sub_meta   <- meta[sub_samples, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(round(sub_counts), sub_meta, design = ~ Condition)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)

  out_list <- list()
  if ("replicative_senescence" %in% sub_meta$Condition) {
    res_RS <- results(dds, contrast = c("Condition", "replicative_senescence", "control"))
    res_RS <- as.data.frame(res_RS)
    res_RS$GeneID <- rownames(res_RS)
    res_RS$Batch <- batch_name
    res_RS$ConditionType <- "replicative_senescence"
    out_list[["RS"]] <- res_RS
  }
  if ("acute_senescence" %in% sub_meta$Condition) {
    res_AS <- results(dds, contrast = c("Condition", "acute_senescence", "control"))
    res_AS <- as.data.frame(res_AS)
    res_AS$GeneID <- rownames(res_AS)
    res_AS$Batch <- batch_name
    res_AS$ConditionType <- "acute_senescence"
    out_list[["AS"]] <- res_AS
  }
  bind_rows(out_list)
}

all_batches <- unique(meta$Batch)
dge_all <- map_dfr(all_batches, run_dge_per_dataset, raw_counts = raw_counts, meta = meta)
cat("Total unique genes across all dataset DGEs:", length(unique(dge_all$GeneID)), "\n")

# ============================================================
# SECTION 6: RANDOM-EFFECTS META-ANALYSIS (REML)
# ============================================================
meta_prep <- dge_all %>% filter(!is.na(log2FoldChange), !is.na(lfcSE), lfcSE > 0)
valid_genes <- meta_prep %>%
  group_by(ConditionType, GeneID) %>%
  summarise(n_valid = sum(!is.na(log2FoldChange)), .groups = "drop") %>%
  filter(n_valid >= 2) %>%
  pull(GeneID)

meta_filtered <- meta_prep %>% filter(GeneID %in% valid_genes)

meta_out <- meta_filtered %>%
  group_by(ConditionType, GeneID) %>%
  group_modify(~{
    yi <- .x$log2FoldChange; sei <- .x$lfcSE
    if (length(yi) >= 2) {
      res <- tryCatch(rma(yi = yi, sei = sei, method = "REML"), error = function(e) NULL)
      if (!is.null(res))
        tibble(meta_log2FC = res$b, meta_pval = res$pval, meta_ci_lb = res$ci.lb, meta_ci_ub = res$ci.ub)
      else
        tibble(meta_log2FC = NA, meta_pval = NA, meta_ci_lb = NA, meta_ci_ub = NA)
    } else tibble(meta_log2FC = NA, meta_pval = NA, meta_ci_lb = NA, meta_ci_ub = NA)
  }) %>%
  ungroup() %>%
  group_by(ConditionType) %>%
  mutate(meta_padj = p.adjust(meta_pval, "BH")) %>%
  ungroup()

meta_annot <- left_join(meta_out, annot_sub, by = "GeneID")
write.csv(meta_annot, file.path(outdir, "Meta_DEG_file.csv"), row.names = FALSE)

# ============================================================
# SECTION 7: COMBINED META-INTEGRATION (RS + AS)
# ============================================================
meta_rs <- meta_annot %>% filter(ConditionType == "replicative_senescence")
meta_as <- meta_annot %>% filter(ConditionType == "acute_senescence")
meta_pair <- inner_join(meta_rs, meta_as, by = "GeneID", suffix = c("_RS", "_AS"))

combined_results <- pmap_dfr(meta_pair, function(...) {
  row <- list(...)
  yi  <- c(row$meta_log2FC_RS, row$meta_log2FC_AS)
  sei <- c((row$meta_ci_ub_RS - row$meta_ci_lb_RS)/(2*1.96),
           (row$meta_ci_ub_AS - row$meta_ci_lb_AS)/(2*1.96))
  if (any(is.na(yi)) || any(is.na(sei)) || any(sei <= 0)) return(tibble(GeneID = row$GeneID, I2 = NA))
  res <- tryCatch(rma(yi = yi, sei = sei, method = "REML"), error = function(e) NULL)
  if (!is.null(res))
    tibble(GeneID = row$GeneID, combined_log2FC = res$b, combined_pval = res$pval, I2 = res$I2)
})
combined_results <- combined_results %>% mutate(combined_padj = p.adjust(combined_pval, "BH"))
write.csv(combined_results, file.path(outdir, "Meta_DEG_Senescence_Combined_Metafor.csv"), row.names = FALSE)

# ============================================================
# SECTION 8: DEG FILTERING AND PSEUDOGENE REMOVAL
# ============================================================
filter_low_quality <- function(df) {
  df %>%
    filter(!is.na(Symbol),
           !grepl("pseudogene|uncharacterized|LOC", Description, ignore.case = TRUE),
           !grepl("^LOC", Symbol, ignore.case = TRUE),
           abs(combined_log2FC) >= 1)
}

sig_up   <- combined_results %>% filter(combined_padj < 0.05, combined_log2FC >= 1)
sig_down <- combined_results %>% filter(combined_padj < 0.05, combined_log2FC <= -1)
write.csv(sig_up, "Significant_Upregulated_logFC1.csv", row.names = FALSE)
write.csv(sig_down, "Significant_Downregulated_logFC1.csv", row.names = FALSE)

clean_up   <- filter_low_quality(sig_up)
clean_down <- filter_low_quality(sig_down)
write.csv(clean_up, "Upregulated_DEGs.csv", row.names = FALSE)
write.csv(clean_down, "Downregulated_DEGs.csv", row.names = FALSE)

# ============================================================
# SECTION 9: FUNCTIONAL ENRICHMENT (GO + REACTOME)
# ============================================================
if (nrow(clean_up) > 0) {
  up_entrez <- bitr(clean_up$Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  go_up <- enrichGO(gene = up_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(go_up), "GO_BP_Upregulated_logfc1.csv", row.names = FALSE)
  reactome_up <- enrichPathway(gene = up_entrez, organism = "human", readable = TRUE)
  write.csv(as.data.frame(reactome_up), "Reactome_Upregulated.csv", row.names = FALSE)
}

if (nrow(clean_down) > 0) {
  down_entrez <- bitr(clean_down$Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  go_down <- enrichGO(gene = down_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(go_down), "GO_BP_Downregulated_logfc1.csv", row.names = FALSE)
  reactome_down <- enrichPathway(gene = down_entrez, organism = "human", readable = TRUE)
  write.csv(as.data.frame(reactome_down), "Reactome_Downregulated.csv", row.names = FALSE)
}

# ============================================================
# SECTION 10: VISUALIZATIONS (VOLCANO + HEATMAP)
# ============================================================
deg <- combined_results
senescence_markers <- c("CDKN2A","CDKN1A","TP53","RB1","LMNB1","CCNB1")  # optional: highlight canonical senescence markers for readability in plots

EnhancedVolcano(
  deg,
  lab = deg$Symbol,
  x = "combined_log2FC",
  y = "combined_padj",
  xlab = bquote(~Log[2]~ "Fold Change"),
  ylab = bquote(~-Log[10]~ "Adjusted P"),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  col = c("grey70", "#4DAF4A", "#377EB8", "#E41A1C"),
  colAlpha = 0.85,
  labSize = 5,
  drawConnectors = TRUE,
  widthConnectors = 0.4,
  colConnectors = "black",
  selectLab = senescence_markers,
  title = "Meta-analysis Volcano Plot"
)
ggsave("Volcano_Meta_DEGs.png", width = 9, height = 7, dpi = 1200)  # choose any file name or output

top_degs <- deg %>% arrange(combined_padj) %>% slice_head(n = 50)
mat <- matrix(top_degs$combined_log2FC, ncol = 1)
rownames(mat) <- top_degs$Symbol
png("Heatmap_Top50.jpg", width = 1200, height = 1600, res = 300)
pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 50 DEGs (Combined Meta-Analysis)",
         fontsize = 10, fontsize_row = 8)
dev.off()

# ============================================================
# SECTION 11: SESSION INFORMATION FOR REPRODUCIBILITY
# ============================================================
sink(file.path(outdir, "session_info.txt"))
sessionInfo()
sink()

# ============================================================
# SECTION 12: FINAL MESSAGE
# ============================================================
cat("✅ Meta-analysis + overlap pipeline finished successfully.\n")
cat("All outputs saved to:", normalizePath(outdir), "\n")
cat("Generated files include:\n")
cat(" - Meta_DEG_results.csv\n")
cat(" - Meta_DEG_Senescence_Combined_Metafor.csv\n")
cat(" - * DEG lists\n")
cat(" - GO and Reactome enrichment tables\n")
cat(" - Volcano and Heatmap plots\n")
cat(" - session_info.txt (environment record)\n")
