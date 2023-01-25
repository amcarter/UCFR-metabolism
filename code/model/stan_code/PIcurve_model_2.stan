data {
    int<lower=0> N;             // total number of observations
    vector[N] fila;             // covariate
    vector[N] epil;             // covariate
    vector[N] light;            // light at each timestep
    vector[N] P;                // productivity (response)
    vector[N] P_sd;             // productivity standard deviation
    int<lower=1> S;             // number of groups (sites)
    int<lower=1,upper=S> ss[N]; // site for each observation
    int new_ts[N];              // vector of 0/1 indicating new site years
}

parameters {
    real<lower=0,upper=1> phi;      // autoregressive coefficient
    real<lower=0> alpha_fila;       // initial slope
    real<lower=0> alpha_epil;       // initial slope
    real<lower=0> Pmax_fila;        // max photosynthetic rate per unit biomass
    real<lower=0> Pmax_epil;        // max photosynthetic rate per unit biomass
    vector<lower=0>[S] alpha_fila_s;// site specific alpha
    vector<lower=0>[S] alpha_epil_s;// site specific alpha
    vector<lower=0>[S] Pmax_fila_s; // site specific Pmax
    vector<lower=0>[S] Pmax_epil_s; // site specific Pmax
    real<lower=0> sigma;            // sd of individual observations
    real<lower=0> tau_a;            // sd of site level alphas
    real<lower=0> tau_Pmax;         // sd of site level P max
    vector[N] mu;                   // underlying mean of the process
}

model {
    // vector[N] mu;

    //priors
    phi ~ beta(1,1);
    alpha_fila ~ normal(0,1);
    alpha_epil ~ normal(0,1);
    alpha_fila_s ~ normal(0,1);
    alpha_epil_s ~ normal(0,1);
    Pmax_fila ~ normal(0,10);
    Pmax_epil ~ normal(0,10);
    Pmax_fila_s ~ normal(0,10);
    Pmax_epil_s ~ normal(0,10);
    sigma ~ normal(0,5);
    tau_a ~ normal(0,0.5);
    tau_Pmax ~ normal(0,0.05);

    for(s in 1:S){
        alpha_fila_s[s] ~ normal(alpha_fila, tau_a);
        alpha_epil_s[s] ~ normal(alpha_epil, tau_a);
        Pmax_fila_s[s] ~ normal(Pmax_fila, tau_Pmax);
        Pmax_epil_s[s] ~ normal(Pmax_epil, tau_Pmax);
    }

    //likelihood
    for(i in 1:N){
        if(new_ts[i] == 1){
            mu[i] ~ normal(P[i], sigma);
        }
        else{
            mu[i] ~ normal(phi*mu[i-1] + Pmax_fila_s[ss[i]] * fila[i] * tanh((alpha_fila_s[ss[i]] * light[i])/Pmax_fila_s[ss[i]]) + Pmax_epil_s[ss[i]] * epil[i] * tanh((alpha_epil_s[ss[i]] * light[i])/Pmax_epil_s[ss[i]]), sigma);
        }
    }

    P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] mu_tilde;
    vector[N] y_tilde;

    for(i in 1:N){
        if(new_ts[i] == 1){
            mu_tilde[i] = normal_rng(P[i], sigma);
        }
        else{
            mu_tilde[i] = normal_rng(phi*mu_tilde[i-1] + Pmax_fila_s[ss[i]] * fila[i] * tanh((alpha_fila_s[ss[i]] * light[i])/Pmax_fila_s[ss[i]]) + Pmax_epil_s[ss[i]] * epil[i] * tanh((alpha_epil_s[ss[i]] * light[i])/Pmax_epil_s[ss[i]]), sigma);
        }
        y_tilde[i] = normal_rng(mu_tilde[i], P_sd[i]);
    }


}

