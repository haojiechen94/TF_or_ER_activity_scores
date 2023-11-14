library(rstan)
library(motifmatchr)
library(GenomicRanges)
library(SummarizedExperiment)
library(BSgenome)
library(TFBSTools)
library(JASPAR2018)
library(BSgenome.Hsapiens.UCSC.hg19)

getJasparMotifs <- function(species = "Homo sapiens", 
                              collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}



#Users can use EAP to generate these data, details see https://github.com/haojiechen94/EAP
load('../Proximal_hypervariable_analysis.RData')
load('../Distal_hypervariable_analysis.RData')


num_of_samples<-0
for(i in colnames(proximal_result)){
    if(grepl('.read_cnt',i)){
        num_of_samples<-num_of_samples+1
    }
}

#define hypervariable peak regions (HVRs)
HVRs<-rbind(proximal_result[proximal_result$fdr<0.01,],distal_result[distal_result$fdr<0.01,])

#get HVRs genomic positions
HVRs_coordinates<-GRanges(seqnames=c(as.character(HVRs$chrom)),
                          ranges=IRanges(start=c(HVRs$start),
                                         end=c(HVRs$end)))

PFMatrixList<-c(getJasparMotifs(species="Homo sapiens",collection="CORE"),
                getJasparMotifs(species="Homo sapiens",collection="UNVALIDATED"))


#TF motifs enrichment
motif_pos<-matchMotifs(PFMatrixList,
                       HVRs_coordinates,
                       genome=BSgenome.Hsapiens.UCSC.hg19, 
                       out='matches')

match_pos<-motifMatches(motif_pos)
match_pos<-as.matrix(match_pos)
names<-list()
for(i in c(1:length(PFMatrixList))){
    names<-append(names,name(PFMatrixList[[i]]))
}

colnames(match_pos)<-as.character(names)

chipseq_weight <- matrix(1, nrow = dim(match_pos)[1], ncol = dim(match_pos)[2])
for(i in 1:length(colnames(match_pos))){
    chipseq_weight[,i] <- ifelse(match_pos[,i], 1, 0)
}
Mask_matrix<-chipseq_weight

data<-HVRs[,colnames(HVRs)[4:(4+num_of_samples-1)]]

X <- t(as.matrix(data))
N <- dim(X)[1]
M <- dim(X)[2]
TF <- dim(match_pos)[2]
data_to_model <- list(N = N, D = D, K = K, X = X, P = Mask_matrix)

#applying HyperBMF 
set.seed(100)
TF_activity_model<-"
data {
    int<lower=0> N; //Number of samples
    int<lower=0> M; //Number of HVRs
    int<lower=0> TF; //Number of TFs
    matrix[M, N] X; //Input matrix
    matrix[M, TF] P; //Prior matrix
}

parameters {
    matrix[M, TF] W; //Weight matrix
    matrix[TF, K] A; //TF activity matrix
    real<lower=0> tau; //Noise term
    vector<lower=0>[TF] alpha; //ARD prior  
}

transformed parameters {
    matrix<lower=0>[M, TF] t_alpha;
    real<lower=0> t_tau;
    for(wr in 1:M){
        for(wc in 1:TF){
            t_alpha[wr, wc] = P[wr, wc] == 1 ? inv(sqrt(alpha[wc])) : 0.01;
        }
    }
    t_tau = inv(sqrt(tau));
}

model {
    tau ~ gamma(1,1);
    to_vector(A) ~ beta(0.5,0.5);
    alpha ~ gamma(1e-3,1e-3);
    for(r in 1:M){
        for(c in 1:TF){
            W[r,c] ~ normal(0, t_alpha[r, c]);
        }
    }
    to_vector(X) ~ normal(to_vector(W*A), t_tau);
}

"


infer_TF_activity_and_regulatory_network <- stan_model(model_code = TF_activity_model)
fit.vb <- vb(infer_TF_activity_and_regulatory_network, 
	data = data_to_model, algorithm = "meanfield",iter = 8000, output_samples = 300, tol_rel_obj = 0.005)


result_matrix <- apply(extract(fit.vb,"Z")[[1]], c(2,3), mean)
TF_activity<-t(result_matrix)
colnames(TF_activity)<-colnames(data)
rownames(TF_activity)<-colnames(match_pos)
write.table(TF_activity,'../TF_activity.txt',sep='\t',quote=F)


result_matrix <- apply(extract(fit.vb,"W")[[1]], c(2,3), mean)
weights<-result_matrix
colnames(weights)<-colnames(match_pos)
write.table(cbind(HVRs[,c('chrom','start','end')],weights),'../weights.txt',sep='\t')

