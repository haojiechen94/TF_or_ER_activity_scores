![R](https://img.shields.io/badge/R-%3E%3D4.4-blue)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-lightgrey)
![License](https://img.shields.io/badge/license-GPL--3-green)

# TregOmic

## Joint Inference of Regulon Activity and Regulatory Potential from Multi-omic Data

**TregOmic** is a computational framework for simultaneously inferring **regulon activities (RA)** and the **regulatory potentials (RP)** of transcriptional regulators by integrating prior information from regulator–target regulatory element and regulator–target gene associations.

TregOmic supports the analysis of **ATAC-seq**, **ChIP-seq**, and **RNA-seq** data, enabling systematic characterization of regulatory programs across biological samples.

<p align="center">
  <img src="https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/Supplementary_figure1.jpg" width="900">
</p>

---

## Overview

TregOmic provides a framework for simultaneously inferring regulon activities and the regulatory potentials of transcriptional regulators by integrating prior information from regulator–target regulatory element and regulator–target gene associations.

Beyond regulon inference, TregOmic supports a range of downstream integrative multi-omic analyses, including:

- Identification of core regulators associated with sample phenotypes.
- Discovery of putative cis-acting genomic alterations associated with regulon activities.
- Discovery of putative trans-acting genomic alterations associated with regulon activities.
- Characterization of post-translational modification sites that may influence regulator activity.
- Integration of genomic, epigenomic, transcriptomic, and proteomic data to investigate regulatory mechanisms underlying cellular phenotypes.

---

## Key Features

### Regulatory Activity Inference

- Quantification of transcriptional regulator activities across samples.
- Simultaneous estimation of regulator-specific regulatory potentials.
- Integration of transcriptional regulator activities and multi-omic data.

### Most variable Feature Identification

- Identification of most variable regulatory elements using ATAC-seq or ChIP-seq data.
- Identification of most variable genes using RNA-seq data.

### Multi-omic Integration

- Mutation–regulon association analysis.
- Cis-effect mutation discovery.
- Trans-effect mutation discovery.
- Post-translational modification analysis.
- Phenotype-associated regulator discovery.

### Visualization

- PCA visualization.
- t-SNE visualization.
- Regulon activity comparison across sample groups.
- Regulatory network exploration.

---

## Installation

### Prerequisites

TregOmic depends on several CRAN and Bioconductor packages.

Install Bioconductor packages:

```r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

BiocManager::install(c(
    "MAnorm2",
    "motifmatchr",
    "GenomicRanges",
    "IRanges",
    "TFBSTools",
    "JASPAR2018",
    "BSgenome",
    "BSgenome.Hsapiens.UCSC.hg19",
    "DESeq2",
    "SummarizedExperiment"
))
```

Install CRAN packages:

```r
install.packages(c(
    "rstan",
    "Rtsne",
    "pcaMethods",
    "RColorBrewer",
    "scales",
    "ggpubr",
    "tidyr"
))
```

### Install TregOmic

Download the source package and install:

```r
install.packages(
    "path/to/TregOmic_0.1.0.tar.gz",
    repos = NULL,
    type = "source"
)
```

---

## Preparing Input Data

### ATAC-seq / ChIP-seq

Raw sequencing data (FASTQ format) can be processed using the Epigenetic Analysis Platform (EAP):

https://github.com/haojiechen94/EAP

After processing, TregOmic requires:

- Proximal peak count matrix
- Distal peak count matrix
- Sample metadata table

Alternatively, users may prepare input files following the MAnorm2 format:

https://github.com/tushiqi/MAnorm2

---

## Example Dataset

Example datasets are available at: Zenodo

---

## Example Workflow

### Step 1. Identify Most Variable Regulatory Elements

```r
library(TregOmic)

HyperChIP_res <- HyperChIP_ATAC_seq(
    "./data/proximal_peak_regions_2000bp.txt",
    "./data/distal_peak_regions_2000bp.txt",
    "./data/GAC_cellines_H3K27ac_ChIP_seq_metadata.csv",
    categorical_variable = "tissue_type",
    top_number_of_PCs = 2,
    perplexity = 0,
    filtered_chromosomes = c("chrX", "chrY", "chrM"),
    fdr_cutoff = 0.01
)
```

### Step 2. Infer Regulon Activity and Regulatory Potential

```r
TregOmic_res <- TregOmic_ATAC_seq(
    HyperChIP_res,
    c(
        "GATA6",
        "HNF4A",
        "HNF4G",
        "TEAD1",
        "TEAD2",
        "TEAD3",
        "TEAD4"
    ),
    fdr_cutoff = 0.01,
    tol_rel_obj = 0.001,
    z_transform = 1,
    genome =
        BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    JASPAR =
        JASPAR2018::JASPAR2018
)
```

### Step 3. Visualize Regulon Activity

```r
boxplot_RA(
    TregOmic_res,
    TR = "HNF4A",
    categorical_variable = "tissue_type"
)
```

---

## Main Output Objects

| Object | Description |
|----------|-------------|
| RA matrix | Regulon activity matrix across samples |
| RP matrix | Regulatory potential matrix |
| Prior matrix | Prior regulator-target association matrix |
| Most variable peak regions | Significant most variable regulatory elements |
| Most variable genes | Top-ranked most vairable genes |
| PCA/t-SNE results | Low-dimensional sample embeddings |

---

## Citation

If you use TregOmic in your research, please cite:

> Chen HJ et al. TregOmic: Joint inference of regulon activity and regulatory potential from multi-omic data. Manuscript in preparation.

---

## Contact

**Haojie Chen**

Email: chenhaojiecompbio@gmail.com

GitHub: https://github.com/haojiechen94/TF_or_ER_activity_scores
