# Senescence Transcriptomic Meta-Analysis.
R pipeline for meta-analysis  transcriptomic integration of cellular senescence datasets.
# Senescence Transcriptome Meta-Analysis

This repository hosts the full pipeline and analysis scripts used to perform
a large-scale, cross-study meta-analysis of cellular senescence transcriptomics.

## Contents
- **pipeline_meta_senescence.R** — unified REML + SVA pipeline
- **session_info.txt** — environment details for reproducibility
  
## 📁 Repository Structure

The following tree shows the internal organization of the repository, including both REML and DESeq2 + SVA analysis pipelines:

```text
Senescence-Transcriptome-MetaAnalysis/
│
├── 📁 Data_of_DEG/
│   ├── 📁 REML_pipeline/                 # Random-Effects Meta-Analysis (REML model)
│   │   ├── R Code REML METHOD FOR META-ANALYSIS
│   │   ├── REML PIPELINE environment session info for R
│   │   ├── README_REML_Info.md
│   │   └── README.md
│   │
│   └── 📁 DESeq2_SVA_pipeline/           # DESeq2 + SVA differential expression analysis
│       ├── R Code RNA-seq Batch Correction (DESeq2 + VST + SVA)
│       ├── DESeq2 and SVA batch correction session info for the R environment
│       ├── README_DESEQ2_SVA_METHOD.md
│       └── README.md
│
├── LICENSE
├── README.md
└── .github/workflows/
    └── codeql.yml

```


## Citation
If you use this pipeline or derived results, please cite:
> Mohd Shahzaib (2025). *Integrative Meta-Analysis of Transcriptomic Networks Reveals Core Signatures and Master Regulators of Cellular Senescence*.

> GitHub: [https://github.com/Mohd-Shahzaib/Senescence-Transcriptome-MetaAnalysis./](https://github.com/Mohd-Shahzaib/Senescence-Transcriptome-MetaAnalysis./)

## License
MIT License © 2025 Mohd Shahzaib
