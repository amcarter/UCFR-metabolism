/////////////////////////////////////////////////////////
// Stan script to run a model of GPP as a function of
// Light and biomass based on a Jasby-Platt PI curve
/////////////////////////////////////////////////////////


data {
    int<lower=0> N;             // total number of observations
    int<lower=0> K;             // number of biomass covariates
    vector[N] light;            // light at each timestep
    matrix[N,K] biomass;        // biomass data
    vector[N] P;                // productivity (response)
    vector[N] P_sd;             // productivity standard deviation
    int<lower=1> S;             // number of groups (sites)
    int<lower=1,upper=S> ss[N]; // site for each observation
    int new_ts[N];              // vector of 0/1 indicating new site years
}

parameters {
    real<lower=0,upper=1> phi;   // autoregressive coefficient
    vector<lower=0>[K] alpha;    // initial slope
    vector<lower=0>[K] Pmax;     // max photosynthetic rate per unit biomass
    matrix<lower=0>[S,K] alpha_s;// site specific alpha
    matrix<lower=0>[S,K] Pmax_s; // site specific Pmax
    real<lower=0> sigma;         // sd of individual observations
    real<lower=0> tau_a;         // sd of site level alphas
    real<lower=0> tau_Pmax;      // sd of site level P max
    vector[N] mu;                // underlying mean of the process
}

transformed parameters{
    matrix[N,K] rates;

    for(k in 1:K){
        for(i in 1:N){
            rates[i,k] = Pmax_s[ss[i],k] * tanh((alpha_s[ss[i],k] * light[i])/Pmax_s[ss[i],k]);
        }
    }
}

model {

    //priors
    phi ~ beta(1,1);
    alpha ~ normal(0,1);
    Pmax ~ normal(0,10);
    sigma ~ normal(0,5);
    tau_a ~ normal(0,0.5);
    tau_Pmax ~ normal(0,0.05);

    for(k in 1:K){
        for(s in 1:S){
            alpha_s[s,k] ~ normal(alpha[k], tau_a);
            Pmax_s[s,k] ~ normal(Pmax[k], tau_Pmax);
        }
    }

    //likelihood
    for(i in 1:N){
        if(new_ts[i] == 1){
            mu[i] ~ normal(P[i], sigma);
        }
        else{
            mu[i] ~ normal(phi*mu[i-1] + rates[i,] * biomass[i,]', sigma);
        }
    }

    P ~ normal(mu, P_sd);
}

// generated quantities {
//     vector[N] mu_tilde;
//     vector[N] y_tilde;
//
//     for(i in 1:N){
//         if(new_ts[i] == 1){
//             mu_tilde[i] = normal_rng(P[i], sigma);
//         }
//         else{
//             mu_tilde[i] = normal_rng(phi*mu_tilde[i-1] + Pmax_fila_s[ss[i]] * fila[i] * tanh((alpha_fila_s[ss[i]] * light[i])/Pmax_fila_s[ss[i]]) + Pmax_epil_s[ss[i]] * epil[i] * tanh((alpha_epil_s[ss[i]] * light[i])/Pmax_epil_s[ss[i]]), sigma);
//         }
//         y_tilde[i] = normal_rng(mu_tilde[i], P_sd[i]);
//     }
//
//
// }
//
