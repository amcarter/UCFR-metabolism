data {
    int<lower=0> N;            // total number of observations
    int<lower=1> K;            // number of covariates
    vector[N] light;           // covariate
    vector[N] biomass;         // covariate
    vector[N] P;               // productivity (response)
    vector[N] P_sd;            // sd on GPP
    int new_ts [N];            // vector of 0/1 indicating new time sequences
}

parameters {
    real<lower=0,upper=1> phi; // ar1 coefficient
    vector[K] gamma;           // population level coefficients
    real<lower=0> sigma;       // process error
    vector[N] mu;              // underlying state
}

model {
    //priors
    gamma ~ normal(0,5);
    phi ~ beta(1,1);
    sigma ~ normal(0,1);

    // initial values
    mu[1] ~ normal(P[1], sigma);
    P[1] ~ normal(mu[1], P_sd[1]);

    for(n in 2:N){ // reinitialize model at the beginning of each site year
        if(new_ts[n] == 1){
            mu[n] ~ normal(P[n], sigma);
        }

        else{
            mu[n] ~ normal(gamma[1] + gamma[2]* light[n] + gamma[3]*biomass[n] + phi*mu[n-1], sigma);
        }

        P[n] ~ normal(mu[n], P_sd[n]);
    }
}

// generated quantities {
//     vector[N] y_tilde;
//     vector[N] mu; //linear predictor
//
//     mu[1] = normal_rng(P[1], sigma);
//     P[1] = normal_rng(mu[1], P_sd[1]);
//
//     for(n in 2:N){
//         mu[n] = normal_rng(gamma[1] + gamma[2]* light[n] + gamma[3]*biomass[n] + phi*mu[n-1], sigma);
//         P[n] = normal_rng(mu[n], P_sd[n]);
//     }
//
// }
//
