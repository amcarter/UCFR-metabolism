//
data {
  int<lower=1> N;  // total number of observations
  vector[N] NPP;  // response variable
  // data for noise-free variables
  int<lower=1> Mme;  // number of groups
  vector[N] Bf_mod;  // noisy values
  vector<lower=0>[N] Bf_se;  // measurement noise
  vector[N] Be_mod;  // noisy values
  vector<lower=0>[N] Be_se;  // measurement noise
  int<lower=1> NCme;  // number of latent correlations
}
transformed data {
}
parameters {
  real<lower=0> mu_f;
  real<lower=0> mu_e;
  real<lower=0> sigma;  // dispersion parameter
  // parameters for noise free variables
  vector[Mme] meanme;  // latent means
  vector<lower=0>[Mme] sdme;  // latent SDs
  matrix[Mme, N] zme;  // standardized latent values
  cholesky_factor_corr[Mme] Lme;  // cholesky factor of the latent correlation matrix
}
transformed parameters {
  matrix[N, Mme] Xme1;  // actual latent values
  // using separate vectors increases efficiency
  vector[N] B_f;
  // using separate vectors increases efficiency
  vector[N] B_e;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual latent values
  Xme1 = rep_matrix(transpose(meanme), N) + transpose(diag_pre_multiply(sdme, Lme) * zme);
  B_f = Xme1[, 1];
  B_e = Xme1[, 2];
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(Lme | 1);
}
model {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += mu_f * B_f[n] + mu_e * B_e[n];
    }
    target += normal_lpdf(NPP | mu, sigma);

  // priors including constants
  target += lprior;
  target += normal_lpdf(Bf_mod | B_f, Bf_se);
  target += normal_lpdf(Be_mod | B_e, Be_se);
  target += std_normal_lpdf(to_vector(zme));
}
generated quantities {
  // obtain latent correlation matrix
  corr_matrix[Mme] Corme = multiply_lower_tri_self_transpose(Lme);
  vector<lower=-1,upper=1>[NCme] corme;
  // extract upper diagonal of correlation matrix
  for (k in 1:Mme) {
    for (j in 1:(k - 1)) {
      corme[choose(k - 1, 2) + j] = Corme[j, k];
    }
  }
}
