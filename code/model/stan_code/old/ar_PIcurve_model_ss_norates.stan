/////////////////////////////////////////////////////////
// Stan script to run a model of GPP as a function of
// Light and biomass based on a Jasby-Platt PI curve
/////////////////////////////////////////////////////////


data {
    int<lower=0> N;             // total number of observations
    int<lower=0> K;             // number of biomass covariates
    vector[N] light;            // light at each timestep
    matrix[N,K] X;              // biomass data
    int<lower=1> S;             // number of groups (sites)
    int<lower=1,upper=S> ss[N]; // site for each observation
    int new_ts[N];              // vector of 0/1 indicating new site years
    vector[N] P;                // productivity (response)
    // vector[N] P_sd;             // productivity standard deviation
}

parameters {
    vector<lower=0>[K] alpha;   // initial slope
    vector<lower=0>[K] Pmax;    // max photosynthetic rate per unit biomass
    real phi;                   // AR term
    real beta;                  // intercept
    vector[S] beta_s;           // site specific intercept
    real<lower=0> tau_b;        // sd of site level intercepts
    real<lower=0> sigma;        // sd of individual observations
    real<lower=0> P_sd;             // productivity standard deviation
    vector[N] mu;                // underlying mean of the process
}

transformed parameters{


}

model {

    //priors
    phi ~ beta(2,2);
    alpha ~ normal(0,1);
    Pmax ~ normal(0,10);
    sigma ~ cauchy(0,5);
    P_sd ~ cauchy(0,5);
    beta ~ normal(0,10);
    tau_b ~ normal(0,1);

    for(s in 1:S){
        beta_s[s] ~ normal(beta, tau_b);
    }

    matrix[N,K] rates;

    for(k in 1:K){
        for(i in 1:N){
            rates[i,k] = Pmax[k] * tanh((alpha[k] * light[i])/Pmax[k]);
        }
    }

    for(i in 1:N){
        if(new_ts[i] == 1){
            mu[i] ~ normal(P[i], sigma);
        }
        else{
            mu[i] ~ normal(beta_s[ss[i]] + mu[i-1] * phi + rates[i,] * X[i,]', sigma);
        }
    }
    //likelihood

    P ~ normal(mu, P_sd);
}

generated quantities {
    vector[N] y_tilde;
    vector[N] mu_tilde;
    vector[N] log_lik;
    matrix[N,K] rates_tilde;

    for(k in 1:K){
        for(i in 1:N){
            rates_tilde[i,k] = Pmax[k] * tanh((alpha[k] * light[i])/Pmax[k]);
        }
    }

    for(i in 1:N){
        if(new_ts[i] == 1){
            mu_tilde[i] = normal_rng(P[i], sigma);
        }
        else{
            mu_tilde[i] = normal_rng(beta_s[ss[i]] + mu[i-1] * phi + rates_tilde[i,] * X[i,]', sigma);
        }
        y_tilde[i] = normal_rng(mu_tilde[i], P_sd);
        log_lik[i] = normal_lpdf(P[i] | mu[i], P_sd);
    }


}

