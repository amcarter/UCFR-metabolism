data {
    int<lower=0> N;             // total number of observations
    int<lower=1> K;             // number of covariates
    vector[N] light;            // covariate
    vector[N] biomass;          // covariate
    vector[N] K600;             // covariate
    vector[N] P;                // productivity (response)
    vector[N] P_sd;             // sd on GPP
    int new_ts [N];             // vector of 0/1 indicating new site years
}

parameters {
    real<lower=0,upper=1> phi; // ar1 coefficient
    vector[K] gamma;           // population level coefficients
    real<lower=0> sigma;       // process error
    vector[N] mu;              // linear predictor
}

model {
    //priors
    gamma ~ normal(0,5);
    phi ~ beta(1,1);
    sigma ~ normal(0,1);

    mu[1] ~ normal(P[1], sigma);

    for(n in 2:N){
        if(new_ts[n] == 1){
            mu[n] ~ normal(P[n], sigma);
        }

        else{
            mu[n] ~ normal(gamma[1] + gamma[2]* light[n] + gamma[3]*biomass[n] + gamma[4]*K600[n] + phi*mu[n-1], sigma);
        }
  }

  P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] mu_tilde; //linear predictor
    vector[N] y_tilde;

    mu_tilde[1] = normal_rng(P[1], sigma);
    y_tilde[1] = normal_rng(mu_tilde[1], P_sd[1]);

    for(n in 2:N){
        if(new_ts[n] == 1){
            mu_tilde[n] = normal_rng(P[n], sigma);
        }
        else{
            mu_tilde[n] = normal_rng(gamma[1] + gamma[2]* light[n] + gamma[3]*biomass[n] + gamma[4] * K600[n] + phi*mu_tilde[n-1], sigma);
        }

        y_tilde[n] = normal_rng(mu_tilde[n], P_sd[n]);
    }

}
