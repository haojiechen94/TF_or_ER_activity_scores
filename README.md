# HyperBMF: jointly infers TF/ER activity and its regulatory network in ATAC/ChIP-seq data
Estimate TF or ER activity based on bulk ATAC/ChIP-seq profiles
![image](https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/workflow.png)

Raw ChIP/ATAC-seq data (in FASTQ format) can be uploaded and processed by EAP (epignomic analysis platform, details see https://github.com/haojiechen94/EAP). After completion of the analyses, two hypervariable analysis results (in proximal regions and distal regions) can be used as input for HyperBMF, or users could prepare the input files according to the Guideline in MAnorm2 (https://github.com/tushiqi/MAnorm2).

Script HyperBMF.R demeonstrates a simple example of transcription regulator activity estimation and its regulation network inference.
Data used for running this example is available in directory [data](https://github.com/haojiechen94/TF_or_ER_activity_scores/tree/master/data) (Ho, Shamaine Wei Ting et al.).

Requirements:
MAnorm2; rstan; motifmatchr; GenomicRanges; SummarizedExperiment; BSgenome; TFBSTools; JASPAR2018; BSgenome.Hsapiens.UCSC.hg19

R packages: Using install.packages("[HyperBMF](https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/HyperBMF_0.0.0.9000.tar.gz)") to install HyperBMF.

Webtools and docker tools are developing.

Citation:<br>
[1] Chen, Haojie et al. “HyperChIP: identification of hypervariable signals across ChIP-seq or ATAC-seq samples.” Genome biology vol. 23,1 62. 28 Feb. 2022, doi:10.1186/s13059-022-02627-9<br>
[2] https://www.cs.helsinki.fi/u/sakaya/tutorial/<br>
[3] Gao, Shang et al. “A Bayesian inference transcription factor activity model for the analysis of single-cell transcriptomes.” Genome research, vol. 31,7 1296–1311. 30 Jun. 2021, doi:10.1101/gr.265595.120<br>
[4] https://mc-stan.org/docs/stan-users-guide/index.html<br>
[5] Ho, Shamaine Wei Ting et al. “Regulatory enhancer profiling of mesenchymal-type gastric cancer reveals subtype-specific epigenomic landscapes and targetable vulnerabilities.” Gut vol. 72,2 (2023): 226-241. doi:10.1136/gutjnl-2021-326483


<p align="center">
  <a href="#">
     <img src="https://api.visitorbadge.io/api/visitors?path=https://github.com/haojiechen94/TF_or_ER_activity_scores" />
   </a>
</p>
