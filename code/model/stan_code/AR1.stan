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
}
transformed parameters {
    vector[N] mu;         // mean
    vector[N] nu;         // prediction for time t
    mu = beta[ss];
    if(K>0){
        mu += X * gamma[2:(K+1)];
    }
    nu = P - mu;

}
model {
    gamma ~ normal(0, 5);
    phi ~ cauchy(0, 2.5);
    sigma ~ cauchy(0, 2.5);
    beta ~ normal(gamma[1], tau);

    for(t in 1:N){
        if(new_ts[t] == 1){
            P[t] ~ normal(mu[t], sigma);

        } else {
            P[t] ~ normal(mu[t] + phi * nu[t-1], sigma);
        }
    }
}

generated quantities{
    vector[N] y_hat;
    vector[N] log_lik;
    for(t in 1:N){
        if(new_ts[t] == 1){
            y_hat[t] = normal_rng(mu[t], sigma);

        } else {
            y_hat[t] = normal_rng(mu[t] + phi * nu[t-1], sigma);
        }

        log_lik[t] = normal_lpdf(P[t] | nu[t], sigma);
    }
}
