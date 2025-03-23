data {
    int<lower=2> N;  // number of observations
    int<lower=1> S;  // number of sites
    int<lower=0> K;  // number of covariates
    matrix[N,K] X;   // covariate matrix
    int<lower=1,upper=S> ss[N];
    vector[N] P;     // observations
    // vector[N] y_sd;
    int<lower=0,upper=1> new_ts[N]; // vector of 0/1 indicating new site years
}
parameters {
    vector[K+1] gamma;
    vector[S] beta;
    real<lower=0> tau;
    real<lower=0> sigma;  // error scale
    // real<lower=0> sigma_obs;  // error scale
    real<lower=0,upper=1> theta;           // lag coefficient
    // vector[N] y_s;
}
transformed parameters {
    vector[N] mu;         // mean
    vector[N] epsilon;    // error terms
    mu = beta[ss] + X * gamma[2:(K+1)];
    for(t in 1:N){
        if(new_ts[t] == 1){
            epsilon[t] = P[t] - mu[t];
        } else {
            epsilon[t] = ( P[t] - mu[t]
                           - theta * epsilon[t - 1] );
        }
    }
}
model {
    gamma ~ normal(0, 5);
    theta ~ cauchy(0, 2.5);
    sigma ~ cauchy(0, 2.5);
    // sigma_obs ~ cauchy(0, 2.5);
    beta ~ normal(gamma[1], tau);

    for (t in 1:N){
        // if(new_ts[t] == 1){
        //     y_s[t] ~ normal(y[t], sigma);
        // }
        if(new_ts[t] == 0){
            P[t] ~ normal(mu[t]
                            + theta * epsilon[t - 1],
                            sigma);
        }
    }

    // y ~ normal(y_s, y_sd);
}
generated quantities{
    vector[N] y_hat;
    vector[N] log_lik;

    for(t in 1:N){
        if(new_ts[t] == 1){
            y_hat[t] = normal_rng(P[t], sigma);
            log_lik[t] = normal_lpdf(P[t] | mu[t], sigma);
        } else{
            y_hat[t] = normal_rng(mu[t] + theta * epsilon[t-1], sigma);
            log_lik[t] = normal_lpdf(P[t] | mu[t] + theta * epsilon[t-1], sigma);
        }
    }
}
