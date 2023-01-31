////////////////////////////////////////////////////////////
// This Stan program fits an AR1 model of GPP as a linear
// function of covariates with fixed observation error
////////////////////////////////////////////////////////////


data {
    int<lower=0> N;             // total number of observations
    int<lower=1> K;             // number of covariates
    int<lower=1> S;             // number of sites
    int<lower=1,upper=S> ss[N]; // site for each observation
    matrix[N,K] X;              // model matrix with covariates
    vector[N] P;                // productivity (response)
    vector[N] P_sd;             // sd on GPP
    int new_ts[N];              // vector of 0/1 indicating new site years
}

parameters {
    real<lower=0,upper=1> phi;  // ar1 coefficient
    vector[K+1] gamma;          // population level coefficients
    vector[S] beta;             // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sigma;        // process error
    vector[N] mu;               // true process mean
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

    for(n in 1:N){
        if(new_ts[n] == 1){
            mu[n] ~ normal(P[n], sigma);
        }
        else{
            mu[n] ~ normal(beta[ss[n]] + X[n,] * gamma[2:K+1] + phi * mu[n-1], sigma);
        }
    }

    P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] mu_tilde; //linear predictor
    vector[N] y_tilde;

    for(n in 1:N){
        if(new_ts[n] == 1){
            mu_tilde[n] = normal_rng(P[n], sigma);
        }
        else{
            mu_tilde[n] = normal_rng(beta[ss[n]] + X[n,] * gamma[2:K+1] + phi * mu_tilde[n-1], sigma);
        }

        y_tilde[n] = normal_rng(mu_tilde[n], P_sd[n]);
    }

}
