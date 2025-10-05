//
// This Stan program defines a model partitioning daily biomass production (NPP)
// across filamentous and epilithic algal functional groups based on variation
// in algal biomass and light using a Monod-type formulation for light saturation
// and a carrying capacity for filamentous biomass to capture the dynamics of self-shading.

data {
    int<lower=1> N;
    // int<lower=0> N_site;  // numnber of sites
    vector[N] NPP;        // daily NPP (GPP - AR)

    // measurement error model for algal biomass
    int<lower=1> Mme;  // number of measurement error groups
    vector[N] B_f_mean;  // modeled filamentous algal biomass
    vector<lower=0>[N] B_f_se;      // se of modeled filamentous algal biomass
    vector[N] B_e;        // epilithic algal biomass
    vector[N] light;      // daily cumulative light (J/m2)
    // int site[N];          // vector of sites
}

parameters {
    real<lower=0> mu_f;     // filamentous algal growth rate
    // real<lower=0> sigma_f;  // sd of filamentous algal growth rates across sites
    real<lower=0> mu_e;     // epilithic algal growth rate
    // real<lower=0> sigma_e;  // sd of epilithic algal growth rates across sites
    //real<lower=0> K_Bf;    // Density dependence of filamentous biomass
    real<lower=0> K_I;      // half saturation for light
    real<lower=0> sigma;    // standard deviation of residual error
    vector[Mme] meanme; // latent means
    vector<lower=0>[Mme] sdme; // latend sds
    matrix[Mme, N] zme // standardized latent values
    cholesky_factor_corr[Mme] Lme; // cholesky factor of the latent correlation matrix
    // vector<lower=0>[N_site] mu_fj;    // filamentous algal growth rate at each site
    // vector<lower=0>[N_site] mu_ej;    // epilithic algal growth rate at each site

    # true values of algal biomass:
    // real B_f[N];    // filamentous algal biomass
}

transformed parameters {
    matrix[N, Mme] Xme; // actual latent values
    vector[N] B_f;

    // compute actual latent values
    Xme = rep_matrix(transpose(meanme), N) + transpose(diag_pre_multiply(sdme, Lme) * zme);
    B_f = Xme[, 1];

}

model {

    // initialize linear predictor:
    vector[N] mu = rep_vector(0.0, N);

    // prior distributions
    mu_f ~ normal(0, 1);
    sigma_f ~ normal(0, 0.1);
    mu_fj ~ gamma(mu_f, sigma_f);
    mu_e ~ normal(0, 1);
    sigma_e ~ normal(0, 0.1);
    mu_ej ~ gamma(mu_e, sigma_e);
    //K_Bf ~ normal(0.2, 0.1);
    K_I ~ normal(5, 1);
    sigma ~ normal(0, 2);
    //B_f ~ normal(40, 50);
    B_f_mean ~ normal(B_f, B_f_se)

    // Model daily NPP as function of algal biomass with normally distributed error
    //NPP ~ normal((mu_f * B_f .* (1/(1 + K_Bf * B_f)) + mu_e * B_e) .* (light./(K_I + light)), sigma);
    for(i in 1:N){
        NPP[i] ~ normal((mu_fj[site[i]] * B_f[i] + mu_ej[site[i]] * B_e[i]) .* (light[i]./(K_I + light[i])), sigma);
    }
}

generated quantities {
    vector[N] NPP_post;

    for(i in 1:N){
//        NPP_post[i] = normal_rng((mu_f*B_f[i] * (1/(1 + K_Bf * B_f[i])) + mu_e*B_e[i]) * (light[i]/(K_I + light[i])), sigma);
        NPP_post[i] = normal_rng((mu_fj[site[i]]*B_f[i] + mu_ej[site[i]]*B_e[i]) * (light[i]/(K_I + light[i])), sigma);
    }
}
//istributed// with mean 'mu' and standard deviation 'sigma'.// // Learn more about model development with Stan at: // //    http://mc-stan.org/users/interfaces/rstan.html //    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started // // The input data is a vector 'y' of length 'N'.// The parameters accepted by the model. Our model// accepts two parameters 'mu' and 'sigma'.// The model to be estimated. We model the output// 'y' to be normally distributed with mean 'mu'// and standard deviation 'sigma'. y ~ normal(mu, sigma);


