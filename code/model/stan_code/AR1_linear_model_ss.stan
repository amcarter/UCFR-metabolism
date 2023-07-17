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
    vector[N] P_sd;             // known sd of observations
    int new_ts[N];              // vector of 0/1 indicating new site years
}

parameters {
    real<lower=0,upper=1> phi;  // ar1 coefficient
    vector[K+1] gamma;          // population level coefficients
    vector[S] beta;             // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sigma;        // standard deviation of process error
    vector[N] P_state;          // underlying true GPP value
}

transformed parameters {
    vector[N] mu;               // underlying mean of process
    // vector[N] y2;               // underlying mean of process

    mu = beta[ss] + X * gamma[2:(K+1)];
    // P_state = P - mu;
}

model {
    //priors
    gamma ~ normal(0,5);
    phi ~ beta(1,1);
    // P_sd ~ cauchy(0,5);
    tau ~ cauchy(0, 2.5);
    sigma ~ cauchy(0,5);
    P_state ~ gamma(1,2);

    for(s in 1:S){
        beta[s] ~ normal(gamma[1], tau);
    }

    for(n in 1:N){
        if(new_ts[n] == 1){
            // restart the AR process on each new time series
            P_state[n] ~ normal(P[n], sigma);
        }
        else{
            // P_state[n] ~ normal(mu[n] + phi * (P[n-1] - mu[n-1]), sigma);
            P_state[n] ~ normal(mu[n] + phi * (P_state[n-1] - mu[n-1]), sigma);
        }

        // Likelihood
        P[n] ~ normal(P_state[n], P_sd[n]);
    }

}

generated quantities {

    vector[N] y_hat;
    vector[N] log_lik;


    for(n in 1:N){
    //     if(new_ts[n] == 1){
    //         // restart the AR process on each new time series
    //         y_hat[n] = normal_rng(P[n], sigma);
    //     }
    //     else{
    //         y_hat[n] = normal_rng(mu[n] + phi * y_hat[n-1], sigma);
    //     }
    //     y_tilde[n] = normal_rng(y_hat[n], P_sd[n]);
        y_hat[n] = normal_rng(P_state[n], P_sd[n]);
        log_lik[n] = normal_lpdf(P[n] | mu[n], P_sd[n]);
    }

}

