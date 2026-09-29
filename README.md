# Senescence Systems Biology

### Reproducible transcriptomic meta-analysis and network-informed discovery in cellular senescence

[![Language: R](https://img.shields.io/badge/Language-R-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Focus: Cellular Senescence](https://img.shields.io/badge/Focus-Cellular%20Senescence-B03A2E)](#scientific-rationale)
[![Methods: REML | DESeq2 | SVA | proteomics ](https://img.shields.io/badge/Methods-REML%20%7C%20DESeq2%20%7C%20SVA-1F6F8B)](#implemented-workflows)
[![License: MIT](https://img.shields.io/badge/License-MIT-2EA44F)](LICENSE)

> **A reproducible R resource for identifying robust senescence-associated transcriptional programmes across heterogeneous datasets and prioritizing candidates for pathway, network, proteomic, and multi-omics validation.**

---

## Scientific rationale

Cellular senescence is a context-dependent cell state rather than a single universal gene signature. Its molecular features vary with cell type, senescence trigger, tissue environment, donor background, experimental design, and technical platform. Consequently, individual transcriptomic studies may capture both generalisable senescence biology and study-specific variation.

This repository provides transparent workflows to distinguish robust cross-study signals from heterogeneous effects. It combines random-effects meta-analysis with batch-aware RNA-seq modelling to prioritize reproducible candidate genes and molecular programmes for systems-level interpretation and experimental validation.

### Central question

> **Which transcriptional programmes are consistently associated with cellular senescence across independent studies after accounting for biological and technical heterogeneity?**

---

## Implemented workflows

| Research layer | Workflow | Core methods | Scientific purpose |
|---|---|---|---|
| **Transcriptomics** | Cross-study meta-analysis | Random-effects model (REML) | Synthesizes gene-level evidence across independent datasets while explicitly modelling between-study heterogeneity |
| **Transcriptomics** | Batch-aware RNA-seq analysis | DESeq2, variance-stabilizing transformation, surrogate-variable analysis (SVA) | Identifies condition-associated transcriptional programmes while reducing unwanted technical variation |
| **Mitochondrial peptide biology** | Comparative Humanin physicochemical profiling | R, `Peptides`, sequence-derived descriptors, hierarchical clustering, reference-centred heatmap | Quantifies cross-species variation in peptide physicochemical features linked to stability, charge, hydrophobicity, and oxidative-stress susceptibility |
| **Structural bioinformatics** | Humanin–receptor protein–peptide docking | HDOCKlite, Python batch automation, top-ranked complex generation | Generates comparative structural hypotheses for Humanin interactions with selected candidate receptors |
| **Structural bioinformatics** | Docked-complex interface energetics | PRODIGY, Python parsing, predicted binding affinity, \(K_d\), interface contacts | Extracts comparative predicted interface-affinity and contact metrics from docking-derived protein–peptide complexes |
| **Reproducibility** | Environment and analysis capture | R session information, scripted workflows | Documents computational dependencies and supports transparent reruns |

---

## Research use

The workflows in this repository are designed to support:

- Identification of robust senescence-associated genes across independent datasets.
- Assessment of consistency and heterogeneity in gene-level effects.
- Transcriptome-guided pathway and gene-set interpretation.
- Prioritization of candidate regulators for protein–protein interaction and network analysis.
- Downstream comparison with proteomic, single-cell, spatial, and other complementary molecular datasets.
- Generation of testable hypotheses in aging, cellular senescence, regenerative medicine, and age-related disease.

---

## Scientific interpretation

This resource supports **evidence synthesis and candidate prioritization**. Its outputs should be interpreted as transcriptomic associations, not as direct evidence of:

- Causal regulatory activity.
- Universal specificity for cellular senescence.
- Protein-level abundance or functional activity.
- Direct therapeutic relevance.

Robust biological conclusions require orthogonal validation through independent cohorts, proteomics, genetic perturbation, functional assays, or single-cell/spatial analyses.

---

## Related research programme

This repository is part of a broader interdisciplinary research programme integrating aging biology, cellular senescence, mitochondrial regulation, transcriptomics, proteomics, network biology, structural bioinformatics, and experimentally grounded stem-cell models.

### Systems biology, transcriptomics, and multi-omics

- **Shahzaib M. et al.** *The interactome era: Integrating RNA-seq, proteomics, and network biology to decode cellular senescence.* **Ageing Research Reviews** (2025).

- **Shahzaib M. et al.** *A conserved regulatory architecture stabilizes cellular senescence across distinct triggers in human fibroblasts.* **GeroScience** (2026).

- **Shahzaib M. et al.** *Integrative Meta-Analysis of Transcriptomic Networks Reveals Core Signatures and Master Regulators of Cellular Senescence.* Computational workflow and research project associated with this repository.

### Mitochondrial peptides and structural bioinformatics

- **Shahzaib M. et al.** *Humanin as an evolutionarily tuned mitochondrial peptide: Insights from mammalian oxidative stress diversity.* **Free Radical Biology and Medicine** (2026).

- **Shahzaib M. et al.** *Network Topology and Interactomic Analysis Reveal the Regulatory Framework of the Humanin Protein Family (MTRNR2Lx Class).* **Biomolecules** (2026).

### Experimental senescence and mitochondrial phenotyping

- **Samiminemati A. et al.** *Mesenchymal Stromal Cell Isolation and Induction of Acute and Replicative Senescence.* In: **Stem Cells and Aging: Methods and Protocols** (2024), pp. 203–215.

- **Samiminemati A. et al.** *Methods to Detect and Compare Cellular and Mitochondrial Changes in Senescent and Healthy Mesenchymal Stem Cells.* In: **Stem Cells and Aging: Methods and Protocols** (2024), pp. 95–124.

Full publication record: [Google Scholar](https://scholar.google.com/citations?user=DKht8zEAAAAJ&hl=en)
---

## Reproducibility

Each workflow directory contains the analysis code, workflow-specific documentation, and R session information. Before reproducing or adapting a pipeline, review the corresponding README, confirm the input-data and metadata structure, and install R packages compatible with the documented environment.

---

## Citation

If you use or adapt this repository, please cite:

> Shahzaib M. (2025). *Integrative Meta-Analysis of Transcriptomic Networks Reveals Core Signatures and Master Regulators of Cellular Senescence.*

**Repository:**  
[https://github.com/Mohd-Shahzaib/Senescence-Transcriptome-MetaAnalysis.](https://github.com/Mohd-Shahzaib/Senescence-Transcriptome-MetaAnalysis.)

---

## License

Released under the [MIT License](LICENSE).

© 2025 Mohd Shahzaib
