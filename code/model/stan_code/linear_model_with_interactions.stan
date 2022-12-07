data {
    int<lower=0> N;             // total number of observations
    int<lower=1> S;             // number of sites
    int<lower=1,upper=S> ss[N]; // site for each observation
    vector[N] light;            // covariate
    vector[N] K600;             // covariate
    vector[N] biomass;          // covariate
    vector[N] P;                // productivity (response)
    vector[N] P_sd;             // sd on GPP
}

parameters {
    vector[5] gamma;          // population level coefficients
    vector[S] beta;             // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sigma;        // process error
    vector[N] mu;               // linear predictor
}

model {
    //priors
    gamma ~ normal(0,5);
    tau ~ cauchy(0, 2.5);
    sigma ~ normal(0,1);

    for(s in 1:S){
        beta[s] ~ normal(gamma[1], tau);
    }

    mu[1] ~ normal(P[1], sigma);

    for(n in 2:N){
        mu[n] ~ normal(beta[ss[n]] + gamma[2]* light[n] + gamma[3] * K600[n] + gamma[4]*biomass[n]+ gamma[5] * biomass[n]* light[n], sigma);
    }

    P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] mu_tilde; //linear predictor
    vector[N] y_tilde;

    mu_tilde[1] = normal_rng(P[1], sigma);
    y_tilde[1] = normal_rng(mu_tilde[1], P_sd[1]);

    for(n in 2:N){
        mu_tilde[n] = normal_rng(beta[ss[n]] + gamma[2]*light[n] + gamma[3] * K600[n] + gamma[4]*biomass[n] + gamma[5]*biomass[n]* light[n], sigma);
        y_tilde[n] = normal_rng(mu_tilde[n], P_sd[n]);
    }

}
