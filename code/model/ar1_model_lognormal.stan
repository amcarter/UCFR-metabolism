
   data {
    int <lower = 0> N; // Sample size
    vector[N] P_obs;
    vector[N] light;
    real mu_obs;
  }

  parameters {
    real b0; // intercept
    real <lower = 0, upper = 1> phi; //phi
    real <lower = 0, upper = 1> b1; //light coefficient
    real <lower = 0> sigma_proc;
    real <lower = 0> sigma_obs;
    vector <lower = 0> [N] P;

  }

  model {
    P[1] ~ lognormal(log(P_obs[1]), 0.05);
    for(t in 2:N){
      P[t] ~ lognormal(log(b0 + phi * P[t-1] + b1 * light[t]), sigma_proc);
      P_obs[t] ~ lognormal(log(P[t]), sigma_obs);
    }

    b0 ~ normal(0,1);
    phi ~ uniform(0,1);
    b1 ~ uniform(0,1);
    sigma_proc ~ normal(0,0.2) T[0,];
    sigma_obs ~ normal(mu_obs, mu_obs/2) T[0,];
  }

