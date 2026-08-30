# TregOmic

**Joint inference of regulon activity and regulatory potential from multi-omic data**

![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)
![R](https://img.shields.io/badge/language-R-276DC3.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

TregOmic is an R package for inferring **regulon activity (RA)** and **regulatory potential (RP)** of transcriptional regulators from epigenomic and transcriptomic data. It integrates prior regulator–target information with sample-level molecular profiles and provides downstream analyses for linking regulator activity to genomic alterations and post-translational modifications.

TregOmic currently supports:

- **ATAC-seq / ChIP-seq**: identification of most variable regulatory regions followed by RA/RP inference using transcription factor motif information.
- **RNA-seq**: identification of most variable genes followed by RA/RP inference using a signed regulator–target network.
- **Mutation–regulon analysis**: identification of putative cis- and trans-acting mutations associated with regulon activity.
- **PTM–regulon analysis**: evaluation of post-translational modification signals associated with regulator activity beyond protein abundance.
- **Single-cell applications**: application to cluster-level or other pseudobulk RNA-seq profiles derived from scRNA-seq data.

---

## Overview

TregOmic is designed around two complementary quantities:

- **Regulon activity (RA):** the inferred activity of a transcriptional regulator in each sample.
- **Regulatory potential (RP):** the inferred contribution of regulatory elements or target genes to regulator-associated variation.

A typical workflow is:

```text
ATAC/ChIP-seq                           RNA-seq
     |                                    |
     v                                    v
HyperChIP_ATAC_seq()                DESeq2_RNA_seq()
     |                                    |
     v                                    v
Most variable regions              Most variable genes
     |                                    |
     +--------------+---------------------+
                    |
                    v
          TregOmic_ATAC_seq() /
          TregOmic_RNA_seq()
                    |
                    v
       Regulon activity + Regulatory potential
                    |
          +---------+----------+
          |                    |
          v                    v
Mutation_affect_regulon()   PTM_affect_regulon()
          |                    |
          v                    v
 cis/trans genomic effects   PTM-associated regulation
```

---

## Installation

### 1. Install Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c(
    "MAnorm2",
    "motifmatchr",
    "GenomicRanges",
    "IRanges",
    "BSgenome",
    "BSgenome.Hsapiens.UCSC.hg19",
    "TFBSTools",
    "JASPAR2018",
    "DESeq2",
    "SummarizedExperiment",
    "pcaMethods"
))
```

### 2. Install CRAN dependencies

```r
install.packages(c(
    "rstan",
    "RColorBrewer",
    "Rtsne",
    "scales",
    "ggplot2",
    "ggpubr",
    "tidyr",
    "patchwork",
    "future",
    "future.apply"
))
```

### 3. Install TregOmic

Download `TregOmic_0.1.0.tar.gz` from this repository and run:

```r
install.packages(
    "TregOmic_0.1.0.tar.gz",
    repos = NULL,
    type = "source"
)

library(TregOmic)
```

> TregOmic uses `rstan` for Bayesian model fitting. A working R/C++ toolchain may therefore be required depending on your operating system and R installation.

---

## Main functions

| Function | Description |
|---|---|
| `HyperChIP_ATAC_seq()` | Identifies most variable regulatory regions from ATAC-seq or ChIP-seq data. |
| `TregOmic_ATAC_seq()` | Infers regulon activity and regulatory potential from epigenomic data. |
| `DESeq2_RNA_seq()` | Identifies most variable genes and prepares RNA-seq data for TregOmic. |
| `TregOmic_RNA_seq()` | Infers regulon activity and regulatory potential from transcriptomic data. |
| `boxplot_RA()` | Visualizes regulon activity across sample groups. |
| `transform_to_matrix()` | Converts TregOmic output into RA, RP, and prior-information matrices. |
| `Mutation_affect_regulon()` | Identifies mutations associated with regulon activity through cis- or trans-effects. |
| `hbar_plot()` | Ranks mutated genes by the number of significantly associated regulon activities. |
| `boxplot_RA2()` | Compares regulon activity between mutant and wild-type samples. |
| `PTM_affect_regulon()` | Identifies regulators whose PTM profiles explain RA beyond protein abundance. |
| `download_data()` | Downloads curated example datasets used in downstream analyses. |

---

# Basic usage

## 1. ATAC-seq / ChIP-seq workflow

### Input data

`HyperChIP_ATAC_seq()` requires:

1. A proximal peak table.
2. A distal peak table.
3. A comma-separated metadata table.

The proximal and distal peak files should contain genomic coordinates in the first three columns, followed by sample-level read-count and occupancy information. Such files can be generated using the [Epigenetic Analysis Platform (EAP)](https://github.com/haojiechen94/EAP).

Example data are also available in the [`data`](https://github.com/haojiechen94/TF_or_ER_activity_scores/tree/master/data) directory.

### Step 1. Identify most variable regulatory regions

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

`HyperChIP_ATAC_seq()` performs HyperChIP-based variability analysis and returns the identified variable regions together with PCA, t-SNE, and sample metadata.

### Step 2. Infer regulon activity and regulatory potential

```r
TregOmic_res <- TregOmic_ATAC_seq(
    HyperChIP_res,
    c(
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

TregOmic scans the selected variable regions for transcription factor motifs and uses the resulting prior matrix during Bayesian inference of RA and RP.

### Step 3. Visualize regulon activity

```r
boxplot_RA(
    TregOmic_res,
    TR = "RUNX2",
    categorical_variable = "tissue_type"
)
```

### Step 4. Export RA and RP matrices

```r
TregOmic_matrix <- transform_to_matrix(TregOmic_res)

RA_matrix <- TregOmic_matrix$Regulon_activity
RP_matrix <- TregOmic_matrix$Rgulatory_potential
prior_matrix <- TregOmic_matrix$Prior_information
```

---

## 2. RNA-seq workflow

### Input data

`DESeq2_RNA_seq()` requires:

- A raw RNA-seq count matrix with genes in rows and samples in columns.
- A comma-separated metadata table with samples in rows.
- A categorical variable used for PCA/t-SNE visualization.

### Step 1. Identify most variable genes

```r
DESeq2_res <- DESeq2_RNA_seq(
    "./data/LUAD_raw_read_counts.txt",
    "./data/LUAD_sample_info.csv",
    categorical_variable = "Stage",
    top_number_of_PCs = 2,
    perplexity = 0,
    gene_expressed_in_min_number_of_samples = 3,
    top_number_of_HVGs = 10000
)
```

### Step 2. Infer regulon activity and regulatory potential

TregOmic uses a signed regulator–target network for transcriptomic analysis. The network should contain regulator (`source`), target gene (`target`), and mode-of-regulation (`mor`) information.

```r
TregOmic_res <- TregOmic_RNA_seq(
    DESeq2_res,
    regulatory_network_path = "./data/human_net.txt",
    iter = 10000,
    output_samples = 300,
    tol_rel_obj = 0.0005,
    z_transform = 1,
    selected_TRs = c("NKX2-1", "CEBPB", "CEBPD", "RUNX2")
)
```

### Step 3. Compare regulon activity across phenotypes

```r
boxplot_RA(
    TregOmic_res,
    TR = "RUNX2",
    categorical_variable = "Stage"
)
```

---

# Advanced applications

## 3. PTM-associated regulation of regulon activity

TregOmic can integrate regulon activity, proteomics, and phosphoproteomics/PTM-omics data to test whether PTM measurements explain additional variation in regulator activity beyond total protein abundance.

For each regulator, the package compares two nested models:

```text
RA ~ Protein
RA ~ Protein + PTM1 + PTM2 + ... + PTMn
```

using an ANOVA model comparison.

### Download the LUAD example dataset

```r
dir.create("./data/LUAD", recursive = TRUE, showWarnings = FALSE)

download_data(
    "https://zenodo.org/records/20754798/files/LUAD_TR_activity.txt",
    "./data/LUAD/LUAD_TR_activity.txt"
)

download_data(
    "https://zenodo.org/records/20754798/files/LUAD_proteomics.txt",
    "./data/LUAD/LUAD_proteomics.txt"
)

download_data(
    "https://zenodo.org/records/20754798/files/LUAD_phosphomics.txt",
    "./data/LUAD/LUAD_phosphomics.txt"
)
```

### Run PTM–regulon analysis

```r
PTM_res <- PTM_affect_regulon(
    "./data/LUAD/LUAD_TR_activity.txt",
    "./data/LUAD/LUAD_proteomics.txt",
    "./data/LUAD/LUAD_phosphomics.txt"
)

head(PTM_res$phos_vs_pro)
```

The returned object currently contains:

```text
PTM_res$phos_vs_pro   regulator-level ANOVA results
PTM_res$RA            regulon activity matrix
PTM_res$Pro           protein abundance matrix
PTM_res$Phos          phosphoproteomics/PTM matrix
```

For phosphoproteomic input, row names should begin with the corresponding gene symbol followed by `|`, because TregOmic uses the first field of the row name to map PTM sites to regulators.

---

## 4. Mutation–regulon association analysis

`Mutation_affect_regulon()` identifies recurrent mutations associated with regulator activity using multivariable linear regression. Regulon activities are restricted away from exactly 0 and 1 and then logit-transformed before model fitting.

The analysis contains two components:

- **cis-effect:** mutation of a transcriptional regulator is tested against the activity of the same regulator.
- **trans-effect:** each recurrently mutated gene is tested against the activities of all regulators.

Optional clinical or molecular covariates can be included in both models.

### Input format

- **Mutation matrix:** genes in rows, samples in columns; mutation status coded numerically, typically `0`/`1`.
- **RA matrix:** transcriptional regulators in rows, samples in columns.
- **Metadata:** samples in rows and covariates in columns; the current function expects a tab-delimited file.

### Download the LSCC example dataset

```r
dir.create("./data/LSCC", recursive = TRUE, showWarnings = FALSE)

download_data(
    "https://zenodo.org/records/22177302/files/LSCC_meta.txt",
    "./data/LSCC/LSCC_meta.txt"
)

download_data(
    "https://zenodo.org/records/22177302/files/LSCC_mutations.txt",
    "./data/LSCC/LSCC_mutations.txt"
)

download_data(
    "https://zenodo.org/records/22177302/files/LSCC_TR_activity.txt",
    "./data/LSCC/LSCC_TR_activity.txt"
)
```

### Run the association analysis

```r
mutation_res <- Mutation_affect_regulon(
    "./data/LSCC/LSCC_mutations.txt",
    "./data/LSCC/LSCC_TR_activity.txt",
    "./data/LSCC/LSCC_meta.txt",
    covariates = c("Age", "Sex", "BMI", "TMB"),
    recurrent_gene_cutoff = 0.05
)
```

Inspect cis-acting associations:

```r
head(mutation_res$cis)
```

Inspect trans-acting associations:

```r
head(mutation_res$trans)
```

Each result contains the estimated mutation coefficient (`beta`), nominal `pval`, and Benjamini–Hochberg adjusted `padj`.

### Compare mutant and wild-type samples

```r
boxplot_RA2(
    "./data/LSCC/LSCC_mutations.txt",
    "./data/LSCC/LSCC_TR_activity.txt",
    TR = "NFE2L2",
    gene = "NFE2L2"
)
```

### Rank recurrent mutations by trans-effects

```r
hbar_plot(
    mutation_res$trans,
    pval_cutoff = 0.001,
    topn = 10,
    color = "purple"
)
```

The bar height represents the number of transcriptional regulators whose activities are significantly associated with mutation of each gene at the specified nominal P-value cutoff.

> In version 0.1.0, trans-effect analysis uses `future::multisession` with 10 workers. Make sure sufficient CPU and memory resources are available.

---

## 5. Application to single-cell RNA-seq data

TregOmic operates on sample-level matrices. For scRNA-seq data, one practical strategy is to aggregate cells into biologically meaningful pseudobulk profiles, such as cell clusters, cell types, or sample-by-cell-type groups, and then apply the RNA-seq workflow.

The example test workflow performs the following steps:

1. Process scRNA-seq data using Seurat.
2. Define cell clusters.
3. Aggregate raw counts within selected clusters to generate cluster-level pseudobulk profiles.
4. Create metadata describing the pseudobulk profiles.
5. Run `DESeq2_RNA_seq()` followed by `TregOmic_RNA_seq()`.

A minimal pseudobulk example is:

```r
library(Seurat)

# After clustering a Seurat object named `obj`:
obj <- JoinLayers(obj)

clusters <- unique(obj$seurat_clusters)

pseudo_counts <- sapply(clusters, function(cl) {
    cells <- colnames(obj)[obj$seurat_clusters == cl]
    rowSums(obj@assays$RNA$counts[, cells, drop = FALSE])
})

colnames(pseudo_counts) <- paste0("cluster_", clusters)
```

The resulting pseudobulk matrix can then be used as input to `DESeq2_RNA_seq()` together with the corresponding metadata.

---

## Output structure

### `TregOmic_ATAC_seq()` and `TregOmic_RNA_seq()`

Both functions return a list containing:

```text
RA_and_RP   regulator-specific inference results
metadata    sample metadata
```

For each regulator in `RA_and_RP`:

```text
TR_activity     inferred regulon activity across samples
weights         inferred regulatory potential
prior_matrix    prior regulator–feature association information
```

Use:

```r
matrix_res <- transform_to_matrix(TregOmic_res)
```

to combine regulator-specific results into matrices for downstream analyses.

---

## Notes on input consistency

Sample identifiers must be consistent across all input matrices and metadata tables. Downstream multi-omic functions automatically retain overlapping samples, but harmonizing sample IDs before analysis is strongly recommended.

For the current version:

- `HyperChIP_ATAC_seq()` and `DESeq2_RNA_seq()` read metadata as **comma-separated** files.
- `Mutation_affect_regulon()` reads metadata as a **tab-delimited** file.
- Regulon activity values used in mutation association analysis are expected to lie between 0 and 1.

---

## Example datasets

Two curated datasets are used in the advanced examples:

- **LUAD proteogenomic data for PTM analysis:** [Zenodo record 20754798](https://zenodo.org/records/20754798)
- **LSCC mutation/regulon data for mutation association analysis:** [Zenodo record 22177302](https://zenodo.org/records/22177302)

Additional epigenomic example files are available in the repository [`data`](https://github.com/haojiechen94/TF_or_ER_activity_scores/tree/master/data) directory.

---

## Citation

If you use TregOmic in your research, please cite:

> Chen HJ et al. **TregOmic: Joint inference of regulon activity and regulatory potential from multi-omic data.** Manuscript in preparation.

Please update this section with the final publication information when available.

---

## Related resources

- **EAP — Epigenetic Analysis Platform:** https://github.com/haojiechen94/EAP
- **TregOmic repository:** https://github.com/haojiechen94/TF_or_ER_activity_scores

---

## Contact

**Haojie Chen**  
Email: chenhaojiecompbio@gmail.com

For bug reports, feature requests, or questions about TregOmic, please use the GitHub Issues page or contact the maintainer.

---

## License

TregOmic is distributed under the MIT License.
