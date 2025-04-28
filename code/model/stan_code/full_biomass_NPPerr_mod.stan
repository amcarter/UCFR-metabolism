data {
    int<lower=1> N;  // total number of observations
    vector[N] epil_meas;  // measured epiithon algae
    vector[N] fila_meas;  // measured filamentous algae
    vector[N] NPP_est;  // estimated NPP
    real<lower=0> NPP_se;  // standard error on NPP estimates
    vector[N] light;  // MODIS derived light
    int<lower=0> Nmi_epil;  // number of missings
    int<lower=0> Nmi_fila;  // number of missings
    int<lower=0> Nmi_NPP;  // number of missings
    array[Nmi_epil] int<lower=1> Jmi_epil;  // positions of missings
    array[Nmi_fila] int<lower=1> Jmi_fila;  // positions of missings
    array[Nmi_NPP] int<lower=1> Jmi_NPP;  // positions of missings
    int<lower=1> K; // number of population-level effects
    // data for splines
    int Ks;  // number of linear effects
    matrix[N, Ks] Xs;  // design matrix for the linear effects
    // data for spline 1
    int nb_1;  // number of bases
    array[nb_1] int knots_1;  // number of knots
    // basis function matrices
    matrix[N, knots_1[1]] Zs_1_1;
    // data for spline 2
    int nb_2;  // number of bases
    array[nb_2] int knots_2;  // number of knots
    // basis function matrices
    matrix[N, knots_2[1]] Zs_2_1;
    matrix[N, knots_2[2]] Zs_2_2;
    matrix[N, knots_2[3]] Zs_2_3;
}
transformed data {
}
parameters {
    vector<lower=0>[Nmi_epil] Ymi_epil;  // estimated missing epils
    vector<lower=0>[Nmi_fila] Ymi_fila;  // estimated missing filas
    vector[Nmi_NPP] Ymi_NPP;  // estimated missing NPPs
    vector[N] NPP;  // true NPPs
    vector<lower=0>[K] b;  // regression coefficients for algal biomass
    // real<lower=0> K_I;    // saturation coefficient for light
    real<lower=0> sigma; // dispersion parameter
    real Intercept_epil;  // temporary intercept for centered predictors
    real Intercept_fila;  // temporary intercept for centered predictors
    vector[Ks] bs_epil;  // unpenalized spline coefficients
    vector[Ks] bs_fila;  // unpenalized spline coefficients
    // parameters for spline 1
    // standardized penalized spline coefficients
    vector[knots_1[1]] zs_1_1_epil;
    vector[knots_1[1]] zs_1_1_fila;
    vector<lower=0>[nb_1] sds_1_epil;  // SDs of penalized spline coefficients
    vector<lower=0>[nb_1] sds_1_fila;  // SDs of penalized spline coefficients
    // parameters for spline 2
    // standardized penalized spline coefficients
    vector[knots_2[1]] zs_2_1_epil;
    vector[knots_2[1]] zs_2_1_fila;
    // standardized penalized spline coefficients
    vector[knots_2[2]] zs_2_2_epil;
    vector[knots_2[2]] zs_2_2_fila;
    // standardized penalized spline coefficients
    vector[knots_2[3]] zs_2_3_epil;
    vector[knots_2[3]] zs_2_3_fila;
    vector<lower=0>[nb_2] sds_2_epil;  // SDs of penalized spline coefficients
    vector<lower=0>[nb_2] sds_2_fila;  // SDs of penalized spline coefficients
    real<lower=0> shape_epil;  // shape parameter
    real<lower=0> shape_fila;  // shape parameter
}
transformed parameters {
    // penalized spline coefficients
    vector[knots_1[1]] s_1_1_epil;
    vector[knots_1[1]] s_1_1_fila;
    // penalized spline coefficients
    vector[knots_2[1]] s_2_1_epil;
    vector[knots_2[1]] s_2_1_fila;
    // penalized spline coefficients
    vector[knots_2[2]] s_2_2_epil;
    vector[knots_2[2]] s_2_2_fila;
    // penalized spline coefficients
    vector[knots_2[3]] s_2_3_epil;
    vector[knots_2[3]] s_2_3_fila;
    real lprior = 0;  // prior contributions to the log posterior
    // compute penalized spline coefficients
    s_1_1_epil = sds_1_epil[1] * zs_1_1_epil;
    s_1_1_fila = sds_1_fila[1] * zs_1_1_fila;
    // compute penalized spline coefficients
    s_2_1_epil = sds_2_epil[1] * zs_2_1_epil;
    s_2_1_fila = sds_2_fila[1] * zs_2_1_fila;
    // compute penalized spline coefficients
    s_2_2_epil = sds_2_epil[2] * zs_2_2_epil;
    s_2_2_fila = sds_2_fila[2] * zs_2_2_fila;
    // compute penalized spline coefficients
    s_2_3_epil = sds_2_epil[3] * zs_2_3_epil;
    s_2_3_fila = sds_2_fila[3] * zs_2_3_fila;
    lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    lprior += student_t_lpdf(Intercept_epil | 3, 0, 2.5);
    lprior += student_t_lpdf(Intercept_fila | 3, 0, 2.5);
    lprior += student_t_lpdf(sds_1_epil | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    lprior += student_t_lpdf(sds_1_fila | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    lprior += student_t_lpdf(sds_2_epil | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    lprior += student_t_lpdf(sds_2_fila | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    lprior += gamma_lpdf(shape_epil | 0.01, 0.01);
    lprior += gamma_lpdf(shape_fila | 0.01, 0.01);
}
model {
    // likelihood including constants
    // vector combining observed and missing responses
    vector[N] Yl_epil = epil_meas;
    vector[N] Yl_fila = fila_meas;
    // initialize linear predictor term
    vector[N] mu_epil = rep_vector(0.0, N);
    vector[N] mu_fila = rep_vector(0.0, N);
    Yl_epil[Jmi_epil] = Ymi_epil;
    Yl_fila[Jmi_fila] = Ymi_fila;
    mu_epil += Intercept_epil + Xs * bs_epil + Zs_1_1 * s_1_1_epil +
    Zs_2_1 * s_2_1_epil + Zs_2_2 * s_2_2_epil + Zs_2_3 * s_2_3_epil;
    mu_fila += Intercept_fila + Xs * bs_fila + Zs_1_1 * s_1_1_fila +
    Zs_2_1 * s_2_1_fila + Zs_2_2 * s_2_2_fila + Zs_2_3 * s_2_3_fila;
    mu_epil = exp(mu_epil);
    target += gamma_lpdf(Yl_epil | shape_epil, shape_epil ./ mu_epil);
    target += gamma_lpdf(Yl_fila | shape_fila, shape_fila ./ mu_fila);
    // Model for NPP as a function of biomass:
    vector[N] Yl_NPP = NPP_est;
    Yl_NPP[Jmi_NPP] = Ymi_NPP;
    vector[N] mu;
    mu = ((b[1])*Yl_epil + (b[2])*Yl_fila ) .* (light ./ (0.5 + light));
    target += normal_lpdf(NPP | mu, sigma);
    target += normal_lpdf(Yl_NPP | NPP, NPP_se);
    // priors including constants
    target += lprior;
    target += std_normal_lpdf(zs_1_1_epil);
    target += std_normal_lpdf(zs_1_1_fila);
    target += std_normal_lpdf(zs_2_1_epil);
    target += std_normal_lpdf(zs_2_1_fila);
    target += std_normal_lpdf(zs_2_2_epil);
    target += std_normal_lpdf(zs_2_2_fila);
    target += std_normal_lpdf(zs_2_3_epil);
    target += std_normal_lpdf(zs_2_3_fila);
}
generated quantities {
    // actual population-level intercept
    real b_Intercept_epil = Intercept_epil;
    real b_Intercept_fila = Intercept_fila;
    // vector combining observed and missing responses
    vector[N] Yl_epil = epil_meas;
    vector[N] Yl_fila = fila_meas;
    vector[N] y_rep;
    vector[N] y_rep2;
    Yl_epil[Jmi_epil] = Ymi_epil;
    Yl_fila[Jmi_fila] = Ymi_fila;

    for (n in 1:N) {
        y_rep2[n] = normal_rng(((b[1]) * Yl_epil[n] + (b[2]) * Yl_fila[n])*(light[n] / (0.5 + light[n])), sigma);
        y_rep[n] = normal_rng(NPP[n], NPP_se);
    }
}
