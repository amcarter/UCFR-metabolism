// STAN code to run linear models of GPP using modeled Biomass data


data {
  int<lower=0> N;
  vector[N] P;
  vector[N] light;
  vector[N] biomass;
  // vector[N] epil_mean;
}

parameters {
  vector[3] beta;
  real<lower=0> sigma;
}

model {

  //priors
  beta ~ normal(0,5);
  sigma~ cauchy(0,1);

  //likelihood
  P ~ normal(beta[1] + beta[2] * light + beta[3] * biomass, sigma);

}

