////////////////////////////////////////////////////////////
// This Stan program fits an linear model of multiple time series that
// all share an underlying AR1 process but diverge from that as a function
// of their covariates and of a random intercept
////////////////////////////////////////////////////////////

data {
    int<lower=1> S;             // number of sites
    int<lower=1> N;             // number of observations at each site
    int<lower=0,upper=N*S> N_mis; // number of missing observations
    int<lower=1,upper=N*S> ii_mis[N_mis]; //indices of missing obs by column
    int<lower=1> K;             // number of covariates
    matrix[N,K] X[S];           // model array with covariates
    vector[S] y[N];             // response
    int<lower=1,upper=N> new_year; // index when a new year starts
}

parameters {
    vector[N_mis] y_mis;
    real<lower=0,upper=1> phi;  // ar1 coefficient
    vector[K+1] beta;           // population level coefficients
    vector[S] beta0;            // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sdp;          // standard deviation of process error
    real<lower=0> sdo;          // standard deviation of observation error
    vector[N] epsilon;          // process error at each timestep
}

transformed parameters{
    vector[S] y_comp[N];
    vector[N] theta;       // latent ar process
    vector[S] mu[N];       // means S for num sites and N for length of time series
    cov_matrix[S] Sigma;
    Sigma = diag_matrix(rep_vector(sdo, S));

    y_comp = y;

    for(i in 1:N_mis){
        int ii_t;
        ii_t = ii_mis[i]%N;
        if(ii_t==0){
            ii_t = N;
        }
        y_comp[ii_t,(ii_mis[i]%/%N+1)] = y_mis[i];
    }

    // initialize the AR first timestep
    theta[1] = sdp / sqrt(1 - phi);

    //complete the AR process
    for(t in 2:N){
        if(t == new_year){
            theta[t] = sdp / sqrt(1 - phi);
        }else{
            theta[t] = theta[t - 1] * phi + epsilon[t];
        }
    }

    // create means from covariates stored in an array with dimensions S x N X K
    for(t in 1:N){
        for(s in 1:S){
            mu[t,s] = beta0[s] + X[s,t] * beta[2:(K+1)] + theta[t];
        }
    }
}


model{
    for(t in 1:N){
        epsilon[t] ~ normal(0, sdp);
        y_comp[t] ~ multi_normal(mu[t], Sigma);
    }

    //priors
    beta ~ normal(0,5);
    phi ~ beta(1,1);
    tau ~ cauchy(0, 2.5);
    sdo ~ normal(0,1);
    sdp ~ normal(0,1);
    for(s in 1:S){
        beta0[s] ~ normal(beta[1], tau);
    }

}

generated quantities {

    vector[S] y_hat[N];

    for(t in 1:N){
        y_hat[t] = multi_normal_rng(mu[t], Sigma);
    }
}
