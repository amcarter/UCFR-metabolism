////////////////////////////////////////////////////////////
// This Stan program fits an linear model of multiple time series that
// all share an underlying AR1 process but diverge from that as a function
// of their covariates
////////////////////////////////////////////////////////////


data {
    int<lower=1> S;             // number of sites
    int<lower=0> N;             // number of observations at each site
    int<lower=1> K;             // number of covariates
    // int<lower=1,upper=S> ss[N]; // site for each observation
    vector[S] X[N,K];           // model array with covariates
    vector[S] y[N];             // productivity (response)
    // vector[S] P_sd[N];          // known sd of observations
    // int miss_dat[N,S];          // matrix of 0/1 indicating missing datai
    // int start_idx[S];           // starting index of each site
    // int end_idx[S];             // ending index of each site
    // int new_year;               // index where 2021 starts
}

parameters {
    real<lower=0,upper=1> phi;  // ar1 coefficient
    row_vector[K+1] gamma;          // population level coefficients
    vector[S] beta;             // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sdp;          // standard deviation of process error
    real<lower=0> sdo;          // standard deviation of observation error
    vector[N] epsilon;          // process error at each timestep
    // vector[S] P_state[N];       // underlying state of process
}

transformed parameters{
    vector[N] theta;       // latent ar process
    vector[S] mu[N];       // means S for num sites and N for length of time series
    cov_matrix[S] Sigma;

    theta[1] = sdp / sqrt(1 - phi);
    for(t in 2:N){
        theta[t] = theta[t - 1] * phi + epsilon[t];
    }
    for(t in 1:N){
        // create means from covariates stored in an array with dimensions S x P x N
        for(s in 1:S){
            mu[t,s] = beta[s] + theta[t];//rep_vector(theta[t], S);
            for(k in 2:K){
                mu[t,s] += gamma[k] * X[t,k,s];
            }
        }
    }

    Sigma = diag_matrix(rep_vector(sdo, S));

}

model{

    for(t in 1:N){
        epsilon[t] ~ normal(0, sdp);

        // y[t] is a vector of length S and Id is SxS identity matrix
        y[t] ~ multi_normal(mu[t], Sigma);
    }

    // for(s in 1:S){
    //     for(n in start_idx[s]:end_idx[s]){
    //
    //         if(n == start_idx[s]|| n == new_year){
    //             //restart the AR process on each new time series
    //             mu[n,s] ~ normal(P[n,s], sigma);
    //
    //         } else {
    //             mu[n,s] ~ normal(beta[s] + X[s][n,] * gamma[2:K+1] + phi * mu[n-1,s], sigma);
    //         }
    //         P[n,s] ~ normal(mu[n,s], P_sd[n,s]);
    //     }
    // }

    //priors
    gamma ~ normal(0,5);
    phi ~ beta(1,1);
    tau ~ cauchy(0, 2.5);
    sdo ~ normal(0,1);
    sdp ~ normal(0,1);
    for(s in 1:S){
        beta[s] ~ normal(gamma[1], tau);
    }

}

generated quantities {

    vector[S] y_hat[N];

    for(t in 1:N){
        y_hat[t] = multi_normal_rng(mu[t], Sigma);
    }
}
