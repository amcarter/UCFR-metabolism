// STAN code to run linear models of GPP using modeled Biomass data


data {
  int<lower=0> N;
  vector[N] P;       // productivity (response)
  // vector[N] site;    // grouping variable
  vector[N] light;   // covariate
  vector[N] biomass; // covariate
  // vector[N] epil_mean;
}

parameters {
  vector[3] beta;       // intercept and covariate parameters
  real<lower=0> sigma;  // sd of innovations
}

model {

  //priors
  beta ~ normal(0,5);
  sigma~ cauchy(0,1);

  //likelihood
  P ~ normal(beta[1] + beta[2] * light + beta[3] * biomass, sigma);

}

