# TF_or_ER_activity_scores
Estimate TF or ER activity based on bulk ATAC/ChIP-seq profiles
![image](https://github.com/haojiechen94/TF_or_ER_activity_scores/blob/master/images/HyperBMF.png)

Raw ChIP/ATAC-seq data (in FASTQ format) can be uploaded and processing by EAP (epignomic analysis platform, details see https://github.com/haojiechen94/EAP). After completion of the analyses, two hypervariable analysis results (in proximal regions and distal regions) can be used as input for HyperBMF.
Or users could prepare the input files according to the Guidelines in MAnorm2 (https://github.com/tushiqi/MAnorm2).

Script HyperBMF.R demeonstrate a simple example of transcription regulator activity estimation and its regulation network inference.
Data for run this example is available in directory data (Ho, Shamaine Wei Ting et al.).

Requirements:
rstan;motifmatchr;GenomicRanges;SummarizedExperiment;BSgenome;TFBSTools;JASPAR2018;BSgenome.Hsapiens.UCSC.hg19

Webtools and docker version are developing.

Citation:<br>
[1] Chen, Haojie et al. “HyperChIP: identification of hypervariable signals across ChIP-seq or ATAC-seq samples.” Genome biology vol. 23,1 62. 28 Feb. 2022, doi:10.1186/s13059-022-02627-9<br>
[2] Gao, Shang et al. “A Bayesian inference transcription factor activity model for the analysis of single-cell transcriptomes.” Genome research, vol. 31,7 1296–1311. 30 Jun. 2021, doi:10.1101/gr.265595.120<br>
[3] https://mc-stan.org/docs/stan-users-guide/index.html<br>
[4] Ho, Shamaine Wei Ting et al. “Regulatory enhancer profiling of mesenchymal-type gastric cancer reveals subtype-specific epigenomic landscapes and targetable vulnerabilities.” Gut vol. 72,2 (2023): 226-241. doi:10.1136/gutjnl-2021-326483
