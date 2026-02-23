data {
    int<lower=2> N;  // number of observations
    int<lower=1> S;  // number of sites
    int<lower=0> K;  // number of covariates
    matrix[N,K] X;   // covariate matrix
    int<lower=1,upper=S> ss[N];
    vector[N] P;     // observations
    int<lower=0,upper=1> new_ts[N]; // vector of 0/1 indicating new site years
}
parameters {
    vector[K+1] gamma;
    vector[S] beta;
    real<lower=0> tau;
    real<lower=0> sigma;         // error scale
    real<lower=0,upper=1> phi;   // AR1 coefficient
    real<lower=0,upper=1> theta; // MA1 coefficient
}
transformed parameters {
    vector[N] mu;         // mean
    vector[N] nu;         // prediction for time t
    vector[N] epsilon;    // error terms
    mu = beta[ss] + X * gamma[2:(K+1)];
    for(t in 1:N){
        if(new_ts[t] == 1){
            nu[t] = mu[t] + phi * mu[t];
            epsilon[t] = P[t] - nu[t];
        } else {
            nu[t] = mu[t] + phi * P[t-1] + theta * epsilon[t-1];
            epsilon[t] = P[t] - nu[t];
        }
    }
}
model {
    gamma ~ normal(0, 5);
    theta ~ cauchy(0, 2.5);
    phi ~ cauchy(0, 2.5);
    sigma ~ cauchy(0, 2.5);
    beta ~ normal(gamma[1], tau);
    epsilon ~ normal(0,sigma);
}
generated quantities{
    vector[N] y_hat;
    vector[N] log_lik;
    for(t in 1:N){
        y_hat[t] = normal_rng(nu[t], sigma);
        log_lik[t] = normal_lpdf(P[t] | nu[t], sigma);
    }
}
