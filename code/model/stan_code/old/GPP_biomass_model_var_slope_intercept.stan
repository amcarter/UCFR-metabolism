data {
  int<lower=0> N;             // total number of observations
  int<lower=1> K;             // number of covariates + intercept
  int<lower=1> S;             // number of groups (sites)
  int<lower=1,upper=S> ss[N]; // site for each observation
  matrix[N,K] X;              // covariate matrix
  vector[N] P;                // productivity (response)
}

parameters {
  vector[K] gamma;       // population level coefficientss
  vector<lower=0>[K] tau;//standard deviation of regression coefficients

  vector[K] beta[S];     //matrix of site level reg coefficients
  real<lower=0> sigma;   // sd of individual observations
}

model {

  vector[N] mu; //linear predictor

  //priors
  gamma ~ normal(0,5);
  tau ~ cauchy(0,2.5);
  sigma ~ gamma(2,1);

  for(s in 1:S){
      beta[s] ~ normal(gamma,tau);
  }

  for(n in 1:N){
      mu[n] = X[n]*beta[ss[n]];
  }

  //likelihood
  P ~ normal(mu, sigma);

}

