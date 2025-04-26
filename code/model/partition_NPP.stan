//
// This Stan program defines a model partitioning daily biomass production (NPP)
// across filamentous and epilithic algal functional groups based on variation
// in algal biomass and light using a Monod-type formulation for light saturation
// and a carrying capacity for filamentous biomass to capture the dynamics of self-shading.

data {
  int<lower=0> N;
  vector[N] NPP; // daily NPP (GPP - AR)
  vector[N] B_f;  // filamentous algal biomass
  vector[N] B_e;  // epilithic algal biomass
}

parameters {
  real mu_f;    // filamentous algal growth rate
  real mu_e;    // epilithic algal growth rate
  real<lower=0> sigma; // standard deviation of residual error
}

model {
  # Model daily NPP as function of algal biomass with normally distributed error
  NPP ~ normal(mu_f * B_f + mu_e * B_e, sigma);
}

