data {
  int<lower=0> N;             // total number of observations
  int<lower=1> K;             // number of covariates
  vector[N] light;            // covariate
  vector[N] biomass;          // covariate
  vector[N] P;                // productivity (response)
  vector[N] P_sd;             // sd on GPP
}

parameters {
  vector[K] gamma;     // population level coefficients

}

model {

  vector[N] mu; //linear predictor

  //priors
  gamma ~ normal(0,5);

  for(n in 1:N){
      mu[n] = gamma[1] + gamma[2]* light[n] + gamma[3]*biomass[n];
      P[n] ~ normal(mu[n], P_sd[n]);
  }

  //likelihood

}

generated quantities {
    vector[N] y_tilde;
    vector[N] mu; //linear predictor

    for(n in 1:N){
        mu[n] = gamma[1] + gamma[2]* light[n] + gamma[3]*biomass[n];
        y_tilde[n] = normal_rng(mu[n], P_sd[n]);
    }

}

