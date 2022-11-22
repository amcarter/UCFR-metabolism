# simulate AR1 data for hierarcical model
library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)

# prep data ####
dd <- read_csv('data/biomass_metab_model_data.csv')
dd <- dd %>%
    mutate(across(starts_with('epil_chla'),
                  function(x) x/max(dd$epil_chla_mgm2_fit, na.rm = T)),
           across(starts_with('epil_gm'),
                  function(x) x/max(dd$epil_gm2_fit, na.rm = T)),
           across(starts_with('fila_chla'),
                  function(x) x/max(dd$fila_chla_mgm2_fit, na.rm = T)),
           across(starts_with('fila_gm'),
                  function(x) x/max(dd$fila_gm2_fit, na.rm = T)))

sites <- unique(dd$site)

# fill in observation data to use as covariates in simulation
dd <- dd %>%
    mutate(year = lubridate::year(date)) %>%
    group_by(site, year) %>%
    mutate(across(-date, ~zoo::na.approx(., na.rm = F))) %>%
    filter(!is.na(GPP), !is.na(light), !is.na(fila_gm2_fit))

sim_dat <- dd %>% filter(site == 'BN')

# set params
n <- nrow(sim_dat)
gamma <- c(1,6,3)       # population level coefficients
phi <- 0.5              # ar1 parameter
sigma_proc <- .1
sigma_obs <- 0.5
# tau <- 1                # standard deviation of intercepts
# sites <- seq(1:6)       # unique sites
K <- length(gamma)      # number of parameters
# S <- length(sites)      # number of sites
# ss <- rep(sites, n/S)   # site for each observation

light <- sim_dat$light
biomass <- sim_dat$fila_chla_mgm2_fit

# beta <- vector(mode = 'double', length = S)
# for(i in 1:S){
#     beta[i] <- rnorm(1, gamma[1], tau)
# }

# simulate response vector
mu <- vector(mode = 'double', length = n)
mu[1] = (gamma[1] + gamma[2]*light[1] + gamma[3]*biomass[1])/(1-phi)

for(i in 2:n){
    mu[i] <- gamma[1] + gamma[2]*light[i] + gamma[3]*biomass[i] + phi*mu[i-1] +
        rnorm(1, 0, sigma_proc)
}

P <- vector(mode = "double", length = n)
P <- rnorm(n, mean = mu, sd = sigma_obs)
P_sd <- rep(sigma_obs, n)
# sim <- data.frame(site = ss, light = light, biomass = biomass, GPP = P)
# ggplot(sim, aes(light, GPP, col = factor(site))) +
#     geom_point()
# compile data for stan

datlist <- list(
    N = n, K = K,
    # S = S, ss = ss,
    light = light, biomass = biomass,
    P = P, P_sd = P_sd
)

# compile stan model
P_mod <- stan_model("code/model/GPP_biomass_model_ar1_hierarchical.stan")

# fit the model
mfit <- sampling(
    P_mod,
    data = datlist#, chains = 1
)

print(mfit, pars = c('gamma', 'phi', 'sigma'))
plot(mfit, pars = c('gamma', 'phi', 'sigma'))
shinystan::launch_shinystan(mfit)

# run on data:
sim_dat <- dd %>% filter(site == 'PL')
GPP <- sim_dat$GPP
GPP_sd <- (sim_dat$GPP.upper - sim_dat$GPP.lower)/3.92

datlist <- list(
    N = nrow(sim_dat), K = 3,
    # S = S, ss = ss,
    light = sim_dat$light, biomass = sim_dat$fila_chla_mgm2_fit,
    P = GPP, P_sd = GPP_sd
)

dfit <- sampling(
    P_mod,
    data = datlist
)

plot(dfit, pars = c('phi', 'gamma', 'sigma'))
print(dfit, pars = c('phi', 'gamma', 'sigma'))
shinystan::launch_shinystan(dfit)


# look at model residuals
preds <- rstan::summary(dfit, pars = 'mu_tilde')$summary %>%
    data.frame() %>%
    select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
preds_poly <- data.frame(date = c(sim_dat$date, rev(sim_dat$date)),
                         y = c(preds$P_mod_lower, rev(preds$P_mod_upper))) %>%
    left_join(sim_dat, by = 'date')

preds <- bind_cols(sim_dat, preds)


ggplot(preds, aes(date, GPP)) +
    geom_point() +
    geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper)) +
    geom_line(aes(y = P_mod), col = 'steelblue', size = 2) +
    geom_polygon(data = preds_poly, aes(date, y),
                 fill = alpha('steelblue', 0.3),
                 col = alpha('steelblue', 0.3))+
    facet_wrap(.~year, scales = 'free_x') +
    theme_bw()
