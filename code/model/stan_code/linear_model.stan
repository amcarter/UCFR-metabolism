////////////////////////////////////////////////////////////////
// This Stan program fits a hierarchical linear regression model
// of GPP as a function of covariates with fixed observation error
////////////////////////////////////////////////////////////////


data {
    int<lower=0> N;             // total number of observations
    int<lower=1> K;             // number of covariates
    int<lower=1> S;             // number of sites
    int<lower=1,upper=S> ss[N]; // site for each observation
    matrix[N,K] X;              // model matrix with covariates
    vector[N] P;                // productivity (response)
    vector[N] P_sd;             // sd on GPP
}

parameters {
    vector[K+1] gamma;          // population level coefficients
    vector[S] beta;             // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sigma;        // process error
    vector[N] mu;               // true process mean
}

model {
    //priors
    gamma ~ normal(0,5);
    tau ~ cauchy(0, 2.5);
    sigma ~ normal(0,1);

    for(s in 1:S){
        beta[s] ~ normal(gamma[1], tau);
    }

    mu ~ normal(beta[ss] + X * gamma[2:K+1], sigma);

    P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] mu_tilde; //linear predictor
    vector[N] y_tilde;

    for(n in 1:N){
        mu_tilde[n] = normal_rng(beta[ss[n]] + X[n,] * gamma[2:K+1], sigma);
        y_tilde[n] = normal_rng(mu_tilde[n], P_sd[n]);
    }

}
