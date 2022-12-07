data {
    int<lower=0> N;        // total number of observations
    vector[N] biomass;     // algal biomass at each timestep
    vector[N] light;       // light at each timestep
    vector[N] K600;
    vector[N] P;           // productivity (response)
    int<lower=1> S;             // number of groups (sites)
    int<lower=1,upper=S> ss[N]; // site for each observation
    // vector[N] P_sd;        // productivity standard deviation
}

parameters {
    real<lower=0> alpha;// initial slope
    real<lower=0> Pmax;        // max photosynthetic rate per unit biomass
    real<lower=0> beta;        // parameter for K600
    vector<lower=0>[S] alpha_s;
    vector<lower=0>[S] Pmax_s;
    real<lower=0> sigma;       // sd of individual observations
    real<lower=0> tau_a;
    real<lower=0> tau_Pmax;
}

model {
    // vector[N] mu;

    //priors
    alpha ~ normal(0,1);
    alpha_s ~ normal(0,1);
    Pmax ~ normal(0,10);
    Pmax_s ~ normal(0,10);
    sigma ~ normal(0,5);
    tau_a ~ normal(0,0.5);
    tau_Pmax ~ normal(0,0.05);

    for(s in 1:S){
        alpha_s[s] ~ normal(alpha, tau_a);
        Pmax_s[s] ~ normal(Pmax, tau_Pmax);
    }

    //likelihood
    for(i in 1:N){
        P[i] ~ normal(Pmax_s[ss[i]] * biomass[i] * tanh((alpha_s[ss[i]] * light[i])/Pmax_s[ss[i]]) + beta * K600[i], sigma);
    }

    // P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] y_tilde;
    vector[N] mu;

    for(i in 1:N){
        mu[i] = Pmax_s[ss[i]] * biomass[i] * tanh((alpha_s[ss[i]] * light[i])/Pmax_s[ss[i]]) + beta * K600[i];
        y_tilde[i] = normal_rng(mu[i], sigma);
    }


}

