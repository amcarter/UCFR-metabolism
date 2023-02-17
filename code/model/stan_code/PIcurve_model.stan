data {
    int<lower=0> N;            // total number of observations
    vector[N] biomass;         // algal biomass at each timestep
    vector[N] light;           // light at each timestep
    // vector[N] K600;
    vector[N] P;               // productivity (response)
    int<lower=1> S;            // number of groups (sites)
    int<lower=1,upper=S> ss[N];// site for each observation
    // vector[N] P_sd;        // productivity standard deviation
}

parameters {
    real<lower=0> alpha;       // initial slope
    real<lower=0> Pmax;        // max photosynthetic rate per unit biomass
    real<lower=0> beta;        // intercept
    // vector<lower=0>[S] alpha_s;
    // vector<lower=0>[S] Pmax_s;
    vector<lower=0>[S] beta_s;
    real<lower=0> sigma;       // sd of individual observations
    // real<lower=0> tau_a;
    // real<lower=0> tau_P;
    real<lower=0> tau_b;
}

model {
    // vector[N] mu;

    //priors
    alpha ~ normal(0,1);
    // alpha_s ~ normal(0,1);
    Pmax ~ normal(0,10);
    // Pmax_s ~ normal(0,10);
    beta ~ normal(0,10);
    beta_s ~ normal(0,10);
    sigma ~ normal(0,5);
    // tau_a ~ normal(0,0.1);
    // tau_P ~ normal(0,0.1);
    tau_b ~ normal(0,1);

    for(s in 1:S){
        // alpha_s[s] ~ normal(alpha, tau_a);
        // Pmax_s[s] ~ normal(Pmax, tau_P);
        beta_s[s] ~ normal(Pmax, tau_b);
    }

    //likelihood
    for(i in 1:N){
        P[i] ~ normal(beta_s[ss[i]] + Pmax * biomass[i] * tanh((alpha * light[i])/Pmax), sigma);
    }

    // P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] y_tilde;
    vector[N] mu;

    for(i in 1:N){
        mu[i] = beta_s[ss[i]] + Pmax * biomass[i] * tanh((alpha * light[i])/Pmax) ;
        y_tilde[i] = normal_rng(mu[i], sigma);
    }


}

