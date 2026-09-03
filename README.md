<div align="center">

# TregOmic

### Joint inference of transcriptional regulator activity and regulatory potential from epigenomic and transcriptomic profiles

**From regulatory activity inference to multi-layer dissection of transcriptional control**

[![R](https://img.shields.io/badge/R-%E2%89%A54.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)](#installation)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](#license)
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](#installation)

[Overview](#overview) •
[Concept](#conceptual-framework) •
[Installation](#installation) •
[Quick start](#quick-start) •
[Multi-omic analyses](#multi-omic-regulatory-analysis) •
[Input & output](#input-and-output) •
[Citation](#citation)

</div>

---

## Overview

**TregOmic** is an R framework for quantitative analysis of transcriptional regulatory programs from **ATAC-seq, ChIP-seq, and RNA-seq** profiles.

Rather than using target-gene or motif enrichment alone as a surrogate for regulator activity, TregOmic jointly estimates:

- **Regulon activity (RA)** — the activity of a transcriptional regulator across biological samples.
- **Regulatory potential (RP)** — the inferred contribution of regulator-associated regulatory elements or target genes to the observed molecular profile.

The framework combines molecular measurements with prior regulator–target relationships and uses a probabilistic model to infer sample-specific regulatory activities together with regulator-specific regulatory potentials.

TregOmic further connects inferred regulon activity to **genomic, proteomic, and post-translational regulatory layers**, enabling systematic investigation of the mechanisms that shape transcriptional regulator activity.

> **In short:** TregOmic moves from *“Which regulators are active?”* to *“Why is a regulator active in a given biological context?”*

---

## Conceptual framework

```text
                             PRIOR KNOWLEDGE
                  ┌─────────────────────────────────┐
                  │ TF–regulatory element links    │
                  │ TF–target gene relationships   │
                  └───────────────┬─────────────────┘
                                  │
                                  ▼
 MOLECULAR PROFILES        ┌───────────────────────────┐
 ┌──────────────────┐      │        TregOmic           │
 │ ATAC-seq         │─────▶│                           │
 │ ChIP-seq         │─────▶│  Joint probabilistic     │
 │ RNA-seq          │─────▶│  inference of RA and RP  │
 └──────────────────┘      └─────────────┬─────────────┘
                                        │
                   ┌────────────────────┴────────────────────┐
                   ▼                                         ▼
          REGULON ACTIVITY                           REGULATORY POTENTIAL
           across samples                        of targets / regulatory elements
                   │
                   ▼
       MULTI-OMIC REGULATORY DISSECTION
 ┌───────────────────────────────────────────────────────────────┐
 │ phenotype associations                                       │
 │ cis- and trans-acting somatic mutations                      │
 │ protein abundance                                             │
 │ phosphorylation / other PTM-associated regulation            │
 │ multi-layer regulatory architecture                          │
 └───────────────────────────────────────────────────────────────┘
```

The current package directly implements workflows for epigenomic and transcriptomic profiles. Single-cell or spatial transcriptomic data can be analyzed after constructing appropriate **cluster-, region-, or pseudobulk-level expression matrices** for the RNA-based workflow.

---

## Why TregOmic?

Transcription-factor abundance is often an incomplete proxy for transcriptional activity. Regulatory output can be altered by chromatin accessibility, enhancer usage, mutations, protein abundance, phosphorylation, interacting factors, and other regulatory processes.

TregOmic is designed around three questions:

1. **Which transcriptional regulators are active?**
2. **Which regulatory elements or target genes contribute to their activity?**
3. **Which molecular alterations may explain changes in regulator activity?**

This design makes TregOmic suitable for studies of:

- lineage-specific transcriptional programs;
- cancer progression and regulatory reprogramming;
- phenotype-associated core transcriptional regulators;
- mutation-driven cis- and trans-regulatory effects;
- phosphorylation-associated regulation of transcription factors;
- multi-layer genomic–proteomic control of regulator activity.

---

## Key capabilities

| Module | Main function | Purpose |
|---|---|---|
| Variable epigenomic feature analysis | `HyperChIP_ATAC_seq()` | Identify highly variable regulatory elements from ATAC-seq or ChIP-seq profiles |
| Epigenomic regulon inference | `TregOmic_ATAC_seq()` | Jointly infer regulon activity and regulatory potential from ATAC/ChIP-seq |
| Variable transcriptomic feature analysis | `DESeq2_RNA_seq()` | Identify highly variable genes and prepare RNA-seq data |
| Transcriptomic regulon inference | `TregOmic_RNA_seq()` | Jointly infer regulon activity and regulatory potential from RNA-seq |
| Matrix extraction | `transform_to_matrix()` | Convert TregOmic results into a regulon activity matrix |
| Phenotype visualization | `boxplot_RA()` | Compare regulator activity across sample groups |
| Mutation association | `Mutation_affect_regulon()` | Identify putative cis- and trans-acting mutations associated with regulon activity |
| PTM association | `PTM_affect_regulon()` | Test whether regulator PTMs explain activity beyond protein abundance |
| Mutation visualization | `boxplot_RA2()` / `hbar_plot()` | Visualize mutation-associated activity changes and trans effects |
| Multi-layer modeling | `Multi_layer_regulation_model()` | Integrate mutations, protein abundance, and phosphosites using Elastic Net |

---

## Installation

### 1. Install Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
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

### 2. Install CRAN dependencies

```r
install.packages(c(
    "rstan",
    "Rtsne",
    "pcaMethods",
    "RColorBrewer",
    "scales",
    "ggpubr",
    "tidyr",
    "patchwork",
    "future.apply",
    "tidyverse",
    "glmnet"
))
```

### 3. Install TregOmic

Clone this repository:

```bash
git clone https://github.com/haojiechen94/TF_or_ER_activity_scores.git
cd TF_or_ER_activity_scores
```

Install the source package in R:

```r
install.packages(
    "TregOmic_0.1.0.tar.gz",
    repos = NULL,
    type = "source"
)
```

Then load the package:

```r
library(TregOmic)
```

---

# Quick start

## Workflow A — ATAC-seq / ChIP-seq

### Step 1. Identify variable regulatory elements

```r
library(TregOmic)

HyperChIP_res <- HyperChIP_ATAC_seq(
    input_proximal = "./data/proximal_peak_regions_2000bp.txt",
    input_distal   = "./data/distal_peak_regions_2000bp.txt",
    metadata       = "./data/sample_metadata.csv",
    categorical_variable = "tissue_type",
    top_number_of_PCs = 2,
    perplexity = 0,
    filtered_chromosomes = c("chrX", "chrY", "chrM"),
    fdr_cutoff = 0.01
)
```

This step identifies highly variable regulatory regions while accounting for the mean–variance relationship in epigenomic count data.

### Step 2. Infer regulon activity and regulatory potential

```r
TregOmic_res <- TregOmic_ATAC_seq(
    HyperChIP_res,
    genes_annotated_with_peaks = c(
        "GATA6", "HNF4A", "HNF4G",
        "TEAD1", "TEAD2", "TEAD3", "TEAD4",
        "RUNX2"
    ),
    fdr_cutoff = 0.01,
    tol_rel_obj = 0.001,
    z_transform = 1,
    genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    JASPAR = JASPAR2018::JASPAR2018
)
```

### Step 3. Compare regulator activity across phenotypes

```r
boxplot_RA(
    TregOmic_res,
    TR = "HNF4A",
    categorical_variable = "tissue_type"
)
```

<p align="center">
  <img src="https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/1.png" width="400">
</p>


```r
boxplot_RA(
    TregOmic_res,
    TR = "RUNX2",
    categorical_variable = "tissue_type"
)
```

<p align="center">
  <img src="https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/2.png" width="400">
</p>

### Step 4. Export a regulon activity matrix

```r
RA_matrix <- transform_to_matrix(TregOmic_res)
```

The resulting matrix can be used for clustering, phenotype association, survival analysis, multi-omic integration, or downstream visualization.

---

## Workflow B — RNA-seq

### Step 1. Prepare RNA-seq data and identify variable genes

```r
DESeq2_res <- DESeq2_RNA_seq(
    input_count_table = "./data/raw_read_counts.txt",
    metadata = "./data/sample_info.txt",
    categorical_variable = "Stage",
    top_number_of_PCs = 2,
    perplexity = 0,
    gene_expressed_in_min_number_of_samples = 3,
    top_number_of_HVGs = 10000
)
```

### Step 2. Infer transcriptomic regulon activity

```r
TregOmic_RNA_res <- TregOmic_RNA_seq(
    DESeq2_res,
    regulatory_network_path = "./data/human_net.txt",
    iter = 10000,
    output_samples = 300,
    tol_rel_obj = 0.0005,
    z_transform = 1,
    selected_TRs = c(
        "NKX2-1",
        "CEBPB",
        "CEBPD",
        "RUNX2"
    )
)
```

### Step 3. Visualize phenotype-associated activity

```r
TregOmic_res_zscore_10000$metadata$Stage<-factor(TregOmic_res_zscore_10000$metadata$Stage,
                                                 levels = c('Stage_I','Stage_II','Stage_III','Stage_IV'))
boxplot_RA(
    TregOmic_RNA_res,
    TR = "RUNX2",
    categorical_variable = "Stage"
)
```
<p align="center">
  <img src="https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/3.png" width="400">
</p>
---

# Multi-omic regulatory analysis

A major goal of TregOmic is to use regulon activity as a quantitative molecular phenotype and determine which upstream regulatory layers explain its variation.

---

## 1. Mutation–regulon activity associations

`Mutation_affect_regulon()` identifies recurrent somatic mutations associated with regulator activity.

The analysis distinguishes:

- **cis effects** — mutation of a regulator gene associated with activity of the same regulator;
- **trans effects** — mutation of one gene associated with activity of another regulator.

A multivariable regression model can incorporate clinical or molecular covariates.

$$ RA \sim \beta_0 + \beta_{mutation} \cdot Muation+\beta_sex \cdot Sex + \beta_{age} \cdot Age + \beta_{BMI} \cdot BMI + \beta_{TMB} \cdot TMB + \epsilon $$

```r
mutation_res <- Mutation_affect_regulon(
    mutation_matrix = "./data/LSCC_mutations.txt",
    RA_matrix = "./data/LSCC_TR_activity.txt",
    metadata = "./data/LSCC_meta.txt",
    covariates = c("Age", "Sex", "BMI", "TMB"),
    recurrent_gene_cutoff = 0.05
)
```

Inspect cis effects:

```r
head(mutation_res$cis)
```

Visualize the activity of a regulator between mutant and wild-type samples:

```r
boxplot_RA2(
    mutation_matrix = "./data/LSCC_mutations.txt",
    RA_matrix = "./data/LSCC_TR_activity.txt",
    gene = "NFE2L2",
    TR = "NFE2L2"
)
```
<p align="center">
  <img src="https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/4.png" width="400">
</p>

Visualize genes ranked by the number of significant trans-regulatory associations:

```r
hbar_plot(
    mutation_res$trans,
    pval_cutoff = 0.001,
    topn = 10,
    color = "purple"
)
```
<p align="center">
  <img src="https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/5.png" width="400">
</p>
---

## 2. PTM-associated regulation of transcriptional activity

Protein abundance alone does not necessarily determine transcription-factor activity. TregOmic therefore provides a model-comparison strategy for identifying regulator PTMs associated with additional variation in regulon activity.

For each regulator, two nested models are compared:

$$ RA \sim PE $$

and

$$
RA \sim PE + PTM_1 + PTM_2 + \cdots + PTM_n
$$

where:

- **RA** = regulon activity;
- **PE** = protein expression;
- **PTM** = regulator-specific modification-site abundance.

The full and reduced models are compared by ANOVA, followed by multiple-testing correction across regulators.

```r
ptm_res <- PTM_affect_regulon(
    RA_path = "./data/LUAD_TR_activity.txt",
    proteomics_path = "./data/LUAD_proteomics.txt",
    phosphomics_path = "./data/LUAD_phosphomics.txt"
)
```

Inspect regulators whose phosphosite profiles explain additional activity variation:

```r
head(ptm_res$phos_vs_pro)
```

Identify phosphosites correlated with the activity of a specific regulator:

```r
sites <- Get_correlated_sites(
    "HNF4A",
    ptm_res$Phos,
    ptm_res$RA
)
```

Visualize protein abundance, phosphosite abundance, and regulon activity:

```r
scatter_plot_Pro_vs_Phos(
    "HNF4A",
    sites$site[1],
    ptm_res$Pro,
    ptm_res$Phos,
    ptm_res$RA
)
```

---

## 3. Multi-layer regulatory architecture

`Multi_layer_regulation_model()` integrates several candidate regulatory layers into a regulator-specific model.

For each eligible transcriptional regulator, the model combines:

- regulator protein abundance;
- recurrent somatic mutations;
- regulator-specific phosphosite abundance.

Regulon activity is bounded between 0 and 1 and is logit-transformed before fitting an **Elastic Net** regression model.

```r
multi_layer_res <- Multi_layer_regulation_model(
    RA_path = "./data/LSCC_TR_activity.txt",
    proteomics_path = "./data/LSCC_proteomics.txt",
    phosphomics_path = "./data/LSCC_phosphomics.txt",
    mutation_path = "./data/LSCC_mutations.txt",
    recurrent_gene_cutoff = 0.10
)
```

For each analyzed regulator, the returned table contains selected features and their coefficients, ranked by absolute coefficient magnitude.

```r
multi_layer_res[["NFE2L2"]]
```

This module is intended to help distinguish regulators whose activities are dominated by protein abundance from those with substantial genetic or post-translational control.

---

# Single-cell and spatial applications

TregOmic currently operates on sample-by-feature matrices rather than directly on individual cells or spots.

Single-cell RNA-seq or spatial transcriptomic data can therefore be used by constructing biologically meaningful **pseudobulk profiles**, for example:

```text
single cells / spatial spots
          │
          ├── cluster
          ├── cell state
          ├── tissue region
          ├── histological annotation
          └── patient × compartment
          │
          ▼
 aggregated gene-expression matrix
          │
          ▼
     DESeq2_RNA_seq()
          │
          ▼
     TregOmic_RNA_seq()
          │
          ▼
 regulon activity across regions / states
```

This strategy enables regulator activity to be compared among cell states, genotypes, tissue regions, pathological compartments, or other aggregated biological units while retaining a statistically stable sample-level representation.

---

# Input and output

## ATAC-seq / ChIP-seq input

TregOmic expects count matrices for proximal and distal regulatory regions together with sample metadata.

Typical input:

```text
proximal_peak_regions_2000bp.txt
distal_peak_regions_2000bp.txt
sample_metadata.csv
```

Raw ATAC-seq or ChIP-seq reads can first be processed using the **Epigenomic Analysis Platform (EAP)**:

https://github.com/haojiechen94/EAP

MAnorm2-compatible count matrices can also be used as input:

https://github.com/tushiqi/MAnorm2

---

## RNA-seq input

The RNA workflow requires:

1. a raw gene-count matrix;
2. sample metadata;
3. a regulator–target network (e.g., DoRothEA).

A typical regulator–target network has the form:

```text
source      target      mor
HNF4A       APOA1       1
HNF4A       CYP2C9      1
TP53        BCL2       -1
```

where `mor` represents the direction of regulation:

- `1` — positive regulation;
- `-1` — negative regulation.

---

## Main TregOmic outputs

The central output of `TregOmic_ATAC_seq()` and `TregOmic_RNA_seq()` contains regulator-specific inference results and sample metadata.

Conceptually, the output includes:

| Output | Description |
|---|---|
| **RA** | Sample-specific regulon activity |
| **RP** | Target- or element-specific regulatory potential |
| **Prior information** | Regulator–target constraints used by the model |
| **Metadata** | Sample phenotype information |
| **Variable features** | Variable regulatory elements or genes identified during preprocessing |

Use:

```r
transform_to_matrix(TregOmic_res)
```

to obtain a conventional regulator × sample activity matrix for downstream analysis.

---

# Statistical model

TregOmic represents the observed molecular profile as a combination of regulator activities and regulatory potentials:

$$
X \approx A W^{T}
$$

where:

- \(X\) is the sample-by-feature molecular profile;
- \(A\) is the sample-by-regulator **regulon activity** matrix;
- \(W\) is the feature-by-regulator **regulatory potential** matrix.

Prior regulator–target relationships constrain the direction and distribution of regulatory potentials. For transcriptomic data, prior relationships may encode activating and repressing effects, whereas epigenomic regulatory associations provide prior support for regulator–regulatory-element connections.

The model is fitted using **Stan**, allowing regulator activity and regulatory potential to be inferred jointly rather than sequentially.

---

# Recommended analysis strategy

A typical TregOmic study can be organized into four analytical layers:

```text
Layer 1 — Molecular heterogeneity
          identify variable peaks / genes

Layer 2 — Regulatory state
          infer RA and RP

Layer 3 — Biological association
          relate RA to phenotype, progression, subtype, or spatial state

Layer 4 — Regulatory mechanism
          mutation → protein → PTM → regulon activity
```

This separation is useful because it distinguishes **observed molecular variation**, **inferred regulatory state**, and **candidate upstream mechanisms**.

---

# Reproducibility

For reproducible analysis, we recommend reporting:

- TregOmic version;
- R and Stan versions;
- genome build and motif database;
- regulator–target network source;
- variable-feature selection thresholds;
- number of model iterations and retained posterior samples;
- transformation strategy (`z_transform`);
- recurrent mutation threshold;
- covariates used in mutation association models.

A fixed random seed is used in core model-fitting procedures where applicable.

---

# Example data

Curated example datasets used by the package workflows are hosted on Zenodo (https://zenodo.org/records/20754798).

Examples include regulon activity, proteomics, phosphoproteomics, mutation, and metadata matrices used to demonstrate downstream multi-omic analyses.


---

# Repository structure

```text
TF_or_ER_activity_scores/
├── data/                     # Example input data
├── images/                   # Figures and workflow illustrations
├── TregOmic_0.1.0.tar.gz     # R source package
└── README.md
```

---

# Function reference

### Core inference

```r
HyperChIP_ATAC_seq()
TregOmic_ATAC_seq()

DESeq2_RNA_seq()
TregOmic_RNA_seq()
```

### Result extraction and visualization

```r
transform_to_matrix()
boxplot_RA()
```

### Mutation analysis

```r
Mutation_affect_regulon()
boxplot_RA2()
hbar_plot()
```

### Post-translational regulation

```r
PTM_affect_regulon()
Get_correlated_sites()
scatter_plot_Pro_vs_Phos()
```

### Multi-layer integration

```r
Multi_layer_regulation_model()
```

---

# Citation

If you use **TregOmic** in your research, please cite:

> **Chen HJ et al.** *TregOmic: Joint inference of regulon activity and regulatory potential from multi-omic data.* Manuscript in preparation.

A formal citation will be added after publication.

If the preprocessing workflow from EAP is used, please also cite the corresponding EAP publication.

---

# Contributing

Issues, bug reports, and suggestions are welcome.

When reporting a problem, please include:

- a minimal reproducible example;
- R version;
- TregOmic version;
- operating system;
- relevant warning or error message.

Please submit issues through the GitHub repository.

---

# License

TregOmic is distributed under the **MIT License**.

---

# Contact

**Haojie Chen**  
Computational biology / cancer epigenomics

- GitHub: https://github.com/haojiechen94
- Repository: https://github.com/haojiechen94/TF_or_ER_activity_scores
- Email: chenhaojiecompbio@gmail.com

---

<div align="center">

### TregOmic

**Quantifying regulatory state. Dissecting regulatory mechanism.**

</div>
