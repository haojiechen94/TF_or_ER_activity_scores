genes_regulatory_model<-"
data {
    int<lower=0> N; //Number of samples
    int<lower=0> M; //Number of HVRs
    int<lower=0> TF; //Number of TFs
    int<lower=0> K; //Number of latent components
    matrix[M, N] X; //Input matrix
    matrix[M, TF] P; //Prior matrix
}

parameters {
    matrix[M, TF] W; //Weight matrix
    matrix[TF, K] L; //Latent matrix
    matrix[K, N] C; //Composition matrix
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
    to_vector(L) ~ normal(0,1);
    to_vector(C) ~ normal(0,1);
    alpha ~ gamma(1e-3,1e-3);
    for(r in 1:M){
        for(c in 1:TF){
            W[r,c] ~ normal(0, t_alpha[r, c]);
        }
    }
    to_vector(X) ~ normal(to_vector(W*L*C), t_tau);
}

"