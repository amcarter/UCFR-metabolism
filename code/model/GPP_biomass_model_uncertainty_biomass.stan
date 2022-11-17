data {
  int<lower=0> N;             // total number of observations
  int<lower=1> K;             // number of covariates
  int<lower=1> S;             // number of groups (sites)
  int<lower=1,upper=S> ss[N]; // site for each observation
  matrix[N,K] X;              // covariate matrix
  vector[N] P;                // productivity (response)
}

parameters {
  vector[K+1] gamma;     // population level coefficients
  real<lower=0> tau;     //standard deviation of intercepts

  vector[S] beta;        // vector of site level intercepts
  real<lower=0> sigma;   // sd of individual observations
}

model {

  vector[N] biomass_draw ~ normal(biomass, biomass_se);

  vector[N] mu; //linear predictor

  //priors
  gamma ~ normal(0,5);
  tau ~ cauchy(0,2.5);
  sigma ~ gamma(2,1);

  for(s in 1:S){
      beta[s] ~ normal(gamma[1],tau);
  }

  for(n in 1:N){
      mu[n] = beta[ss[n]] + biomass_draw[n] * gamma[2] + X[n]*gamma[2:(K+1)];
  }

  //likelihood
  P ~ normal(mu, sigma);

}

