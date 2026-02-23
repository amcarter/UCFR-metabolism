data {
    int<lower=0> N;             // total number of observations
    int<lower=1> K;             // number of biomass categories
    int<lower=1> S;             // number of sites
    int<lower=1,upper=S> ss[N]; // site for each observation
    vector[N] K600;             // covariate
    vector[N] light;            // covariate
    matrix[N,K] biomass;        // covariate
    vector[N] P;                // productivity (response)
    vector[N] P_sd;             // sd on GPP
    int new_ts[N];              // vector of 0/1 indicating new site years
}

transformed data{
    vector[N] log_P;

    log_P = log(P);
}

parameters {
    real<lower=0,upper=1> phi;  // ar1 coefficient
    vector[3+K] gamma;          // population level coefficients
    vector[S] beta;             // site level intercepts
    real<lower=0> tau;          // variation in site intercepts
    real<lower=0> sigma;        // process error
    vector[N] mu;               // linear predictor
}

model {
    //priors
    gamma ~ normal(0,5);
    phi ~ beta(1,1);
    tau ~ cauchy(0, 2.5);
    sigma ~ normal(0,1);

    for(s in 1:S){
        beta[s] ~ normal(gamma[1], tau);
    }

    mu[1] ~ normal(log_P[1], sigma);

    for(n in 2:N){
        if(new_ts[n] == 1){
            mu[n] ~ normal(log_P[n], sigma);
        }
        else{
            mu[n] ~ normal(beta[ss[n]] + gamma[2]*K600[n] + gamma[3]* light[n] + biomass[n]*gamma[4:(3+K)]  + phi*mu[n-1], sigma);
        }
    }

    P ~ normal(exp(mu), P_sd);
}

generated quantities {
    vector[N] mu_tilde; //linear predictor
    vector[N] y_tilde;

    mu_tilde[1] = normal_rng(log_P[1], sigma);
    y_tilde[1] = normal_rng(exp(mu_tilde[1]), P_sd[1]);

    for(n in 2:N){
        if(new_ts[n] == 1){
            mu_tilde[n] = normal_rng(log_P[n], sigma);
        }
        else{
            mu_tilde[n] = normal_rng(beta[ss[n]] + gamma[2]*K600[n] + gamma[3]*light[n] + biomass[n]*gamma[4:(3+K)] + phi*mu_tilde[n-1], sigma);
        }

        y_tilde[n] = normal_rng(exp(mu_tilde[n]), P_sd[n]);
    }

}
