////////////////////////////////////////////////////////////
// This Stan program fits an AR1 model of GPP as a linear
// function of covariates with fixed observation error
////////////////////////////////////////////////////////////


data {
    int<lower=1> S;             // number of sites
    int<lower=0> N;             // number of observations at each site
    int<lower=1> K;             // number of covariates
    // int<lower=1,upper=S> ss[N]; // site for each observation
    matrix[N,K] X[S];           // model array with covariates
    matrix[N,S] P;              // productivity (response)
    matrix[N,S] P_sd;           // known sd of observations
    // int miss_dat[N,S];          // matrix of 0/1 indicating missing data
    int start_idx[S];           // starting index of each site
    int end_idx[S];             // ending index of each site
    int new_year;               // index where 2021 starts
}

parameters {
    real<lower=0,upper=1> phi;  // ar1 coefficient
    vector[K+1] gamma;          // population level coefficients
    vector[S] beta;             // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sigma;        // standard deviation of process error
    matrix[N,S] mu;             // underlying mean of process
}


model {
    //priors
    gamma ~ normal(0,5);
    phi ~ beta(1,1);
    tau ~ cauchy(0, 2.5);
    sigma ~ normal(0,1);

    for(s in 1:S){
        beta[s] ~ normal(gamma[1], tau);
    }

    for(s in 1:S){
        for(n in start_idx[s]:end_idx[s]){

            if(n == start_idx[s]|| n == new_year){
                //restart the AR process on each new time series
                mu[n,s] ~ normal(P[n,s], sigma);

            } else {
                mu[n,s] ~ normal(beta[s] + X[s][n,] * gamma[2:K+1] + phi * mu[n-1,s], sigma);
            }
            P[n,s] ~ normal(mu[n,s], P_sd[n,s]);
        }
    }

}

generated quantities {

    matrix[N,S] y_tilde;
    matrix[N,S] log_lik;

    for(s in 1:S){
        for(n in start_idx[s]:end_idx[s]) {
            y_tilde[n,s] = normal_rng(mu[n,s], P_sd[n,s]);
            log_lik[n,s] = normal_lpdf(P[n,s] | mu[n,s], P_sd[n,s]);
        }
    }
}
