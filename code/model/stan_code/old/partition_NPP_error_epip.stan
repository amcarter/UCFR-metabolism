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
    int<lower=1> Mme;          // number of covariates with measurement error (2)
    int<lower=1> NCme;         // number of latent correlations (1)
    vector[N] B_f_mean;        // modeled filamentous algal biomass
    vector<lower=0>[N] B_f_se; // se of modeled filamentous algal biomass
    vector[N] B_ep_mean;        // modeled epilithic algal biomass
    vector<lower=0>[N] B_ep_se; // se of modeled epilithic algal biomass
    vector[N] B_e_mean;        // modeled epilithic algal biomass
    vector<lower=0>[N] B_e_se; // se of modeled epilithic algal biomass
    vector[N] light;           // daily cumulative light (J/m2)
    // int site[N];               // vector of sites
}

parameters {
    real<lower=0> mu_f;      // filamentous algal growth rate
    // real<lower=0> sigma_f;  // sd of filamentous algal growth rates across sites
    real<lower=0> mu_ep;      // epilithic algal growth rate
    real<lower=0> mu_e;      // epilithic algal growth rate
    // real<lower=0> sigma_e;  // sd of epilithic algal growth rates across sites
    //real<lower=0> K_Bf;    // Density dependence of filamentous biomass
    real<lower=0> K_I;       // half saturation for light
    real<lower=0> sigma;     // standard deviation of residual error

    // measurement error latent variables:
    vector[Mme] meanme;        // latent means
    vector<lower=0>[Mme] sdme; // latend sds
    matrix[Mme, N] zme;        // standardized latent values
    cholesky_factor_corr[Mme] Lme; // cholesky factor of the latent correlation matrix
    // vector<lower=0>[N_site] mu_fj;    // filamentous algal growth rate at each site
    // vector<lower=0>[N_site] mu_ej;    // epilithic algal growth rate at each site

    // true values of algal biomass:
    // real B_f[N];    // filamentous algal biomass
}

transformed parameters {
    matrix[N, Mme] Xme; // actual latent values for algal biomass
    // using separate vectors increases efficiency
    vector[N] B_f;      // true filamentous biomass
    vector[N] B_ep;      // true epilithic biomass
    vector[N] B_e;      // true epilithic biomass
    real lprior = 0;    // prior contributions to the log posterior

    // compute actual latent values
    Xme = rep_matrix(transpose(meanme), N) + transpose(diag_pre_multiply(sdme, Lme) * zme);
    B_f = Xme[, 1];
    B_ep = Xme[, 2];
    B_e = Xme[, 3];

    // prior distributions:
    lprior += normal_lpdf(mu_f | 0, 1);
    lprior += normal_lpdf(mu_ep | 0, 1);
    lprior += normal_lpdf(mu_e | 0, 1);
    lprior += normal_lpdf(K_I | 5, 1);
    lprior += student_t_lpdf(sigma | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    lprior += normal_lpdf(meanme[1] | 50, 10);
    lprior += normal_lpdf(meanme[1] | 30, 10);
    lprior += student_t_lpdf(sdme | 3, 0, 2.5) - 2 * student_t_lccdf(0 | 3, 0, 2.5);
    lprior += lkj_corr_cholesky_lpdf(Lme | 1);
}

model {

    // initialize linear predictor:
    vector[N] mu = rep_vector(0.0, N);
    for(n in 1:N){
        // add more terms to the linear predictor
        // Model daily NPP as function of algal biomass with normally distributed error
        mu[n] += (mu_f * B_f[n] + mu_ep * B_ep[n] + mu_e * B_e[n]) .* (light[n]./(K_I + light[n]));
    }
    //NPP ~ normal((mu_f * B_f .* (1/(1 + K_Bf * B_f)) + mu_e * B_e) .* (light./(K_I + light)), sigma);

    //likelihood:
    target += normal_lpdf(NPP | mu, sigma);

    // prior distributions
    target += lprior;
    target += normal_lpdf(B_f_mean | B_f, B_f_se);
    target += normal_lpdf(B_ep_mean | B_ep, B_ep_se);
    target += normal_lpdf(B_e_mean | B_e, B_e_se);
    target += std_normal_lpdf(to_vector(zme));
}

generated quantities {
    // obtain latent correlation matrix
    corr_matrix[Mme] Corme = multiply_lower_tri_self_transpose(Lme);
    vector<lower=-1,upper=1>[NCme] corme;
    //exact upper diagonal of correlation matrix
    for (k in 1:Mme) {
        for (j in 1:(k-1)) {
            corme[choose(k-1, 2) + j] = Corme[j,k];
        }
    }

    vector[N] NPP_post;

    for(i in 1:N){
//        NPP_post[i] = normal_rng((mu_f*B_f[i] * (1/(1 + K_Bf * B_f[i])) + mu_e*B_e[i]) * (light[i]/(K_I + light[i])), sigma);
        NPP_post[i] = normal_rng((mu_f*B_f[i] + mu_ep*B_ep[i] + mu_e*B_e[i]) * (light[i]/(K_I + light[i])), sigma);
    }
}


