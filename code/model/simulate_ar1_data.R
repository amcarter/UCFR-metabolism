# simulate AR1 data for hierarcical model
library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
rle2 <- function(x){ # function for breaking site years to restart AR model

    r <- rle(x)
    ends <- cumsum(r$lengths)

    r <- tibble(values = r$values,
                starts = c(1, ends[-length(ends)] + 1),
                stops = ends,
                lengths = r$lengths)

    return(r)
}

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

# evaluate GPP for normality: ####
AID::boxcoxnc(dd$GPP)

# fill in observation data to use as covariates in simulation
dd <- dd %>%
    mutate(year = lubridate::year(date)) %>%
    group_by(site, year) %>%
    mutate(across(-date, ~zoo::na.approx(., na.rm = F))) %>%
    filter(!is.na(GPP), !is.na(light), !is.na(fila_gm2_fit)) %>%
    mutate(site = factor(site, levels = c('PL','DL','GR','GC','BM','BN')))

# individual site model ########################################################
# Simulate AR1 data ####
sim_dat <- dd %>% filter(site == 'BN')
# set params
n <- nrow(sim_dat)
gamma <- c(1,6,3,0.5)       # population level coefficients
phi <- 0.5              # ar1 parameter
sigma_proc <- 1
sigma_obs <- 0.5
K <- length(gamma)      # number of parameters

light <- sim_dat$light
biomass <- sim_dat$fila_chla_mgm2_fit
K600 <- sim_dat$K600

# simulate response vector
mu <- vector(mode = 'double', length = n)
mu[1] = (gamma[1] + gamma[2]*light[1] + gamma[3]*biomass[1] + gamma[4] * K600[1])/(1-phi)

for(i in 2:n){
    mu[i] <- gamma[1] + gamma[2]*light[i] + gamma[3]*biomass[i] +
        gamma[4] * K600[i] + phi*mu[i-1] +
        rnorm(1, 0, sigma_proc)
}

P <- vector(mode = "double", length = n)
P <- rnorm(n, mean = mu, sd = sigma_obs)
P_sd <- rep(sigma_obs, n)
plot(sim_dat$date, sim_dat$GPP)
lines(sim_dat$date, P)
new_ts <- rep(0, n)
starts <- rle2(paste0(sim_dat$site, sim_dat$year))$starts
new_ts[starts] <- 1

datlist <- list(
    N = n, K = K,
    light = light, biomass = biomass, K600 = K600,
    P = P, P_sd = P_sd,
    new_ts = new_ts
)

# compile stan model
P_mod <- stan_model("code/model/stan_code/GPP_biomass_model_ar1_K600.stan")

# fit the model
mfit <- sampling(
    P_mod,
    data = datlist#, chains = 1
)

print(mfit, pars = c('gamma', 'phi', 'sigma'))
plot(mfit, pars = c('gamma', 'phi', 'sigma'))
pairs(mfit, pars = c('gamma', 'phi', 'sigma'))
# shinystan::launch_shinystan(mfit)


# run on real data: ####
dat <- dd %>% filter(site == 'BN')
GPP <- dat$GPP
GPP_sd <- (dat$GPP.upper - dat$GPP.lower)/3.92
new_ts <- rep(0, nrow(dat))
new_ts[rle2(dat$year)$starts] <- 1

datlist <- list(
    N = nrow(dat), K = 4,
    light = dat$light,
    biomass = dat$fila_chla_mgm2_fit,
    K600 = dat$K600,
    P = GPP, P_sd = GPP_sd,
    new_ts = new_ts
)

dfit <- sampling(
    P_mod,
    data = datlist
)

# plot(dfit, pars = c('phi', 'gamma', 'sigma'))
print(dfit, pars = c('phi', 'gamma', 'sigma'))
# shinystan::launch_shinystan(dfit)


# look at model residuals
preds <- rstan::summary(dfit, pars = 'y_tilde')$summary %>%
    data.frame() %>%
    select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
preds_poly <- data.frame(date = c(dat$date, rev(dat$date)),
                         y = c(preds$P_mod_lower, rev(preds$P_mod_upper))) %>%
    left_join(dat, by = 'date')

preds <- bind_cols(dat, preds)


ggplot(preds, aes(date, GPP)) +
    geom_point() +
    geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper)) +
    geom_line(aes(y = P_mod), col = 'steelblue', linewidth = 1.2) +
    geom_polygon(data = preds_poly, aes(date, y),
                 fill = alpha('steelblue', 0.3),
                 col = alpha('steelblue', 0.3))+
    facet_wrap(.~year, scales = 'free_x') +
    theme_bw()

# ggpubr::ggarrange(PL, DL, GR, GC, BM, BN)

# hierarchical model ###########################################################
# Simulate AR1 data ####

sim_dat <- dd
# set params
n <- nrow(sim_dat)
gamma <- c(-10,6,3,0.2)         # population level coefficients
phi <- 0.5                      # ar1 parameter
sigma_proc <- 1
sigma_obs <- 0.5
tau <- 1                        # standard deviation of intercepts
K <- length(gamma)              # number of parameters
S <- length(unique(dd$site))    # number of sites
ss <- as.numeric(dd$site)       # site for each observation

light <- sim_dat$light
biomass <- sim_dat$fila_chla_mgm2_fit
K600 <- sim_dat$K600

# select random intercept for each site
beta <- vector(mode = 'double', length = S)
for(i in 1:S){
    beta[i] <- rnorm(1, gamma[1], tau)
}

# simulate response vector
mu <- vector(mode = 'double', length = n)
mu[1] = (beta[ss[1]] + gamma[2]*light[1] + gamma[3]*biomass[1] +
             gamma[4] * K600[1])/(1-phi)

for(i in 2:n){
    mu[i] <- beta[ss[i]] + gamma[2]*light[i] + gamma[3]*biomass[i] +
        gamma[4] * K600[i] + phi*mu[i-1] +
        rnorm(1, 0, sigma_proc)
}

P <- vector(mode = "double", length = n)
P <- rnorm(n, mean = mu, sd = sigma_obs)
P_sd <- rep(sigma_obs, n)
sim <- data.frame(site = ss, light = light, biomass = biomass, GPP = P)
ggplot(sim, aes(light, GPP, col = factor(site))) +
    geom_point()
new_ts <- rep(0, n)
new_ts[rle2(paste0(sim_dat$year, sim_dat$site))$starts] <- 1
# compile data for stan

datlist <- list(
    N = n, K = K,
    S = S, ss = ss,
    light = light, biomass = biomass, K600 = K600,
    P = P, P_sd = P_sd,
    new_ts = new_ts
)

# compile stan model
P_mod <- stan_model("code/model/stan_code/GPP_biomass_model_ar1_hierarchical.stan")

# fit the model
mfit <- sampling(
    P_mod,
    data = datlist#, chains = 1
)

print(mfit, pars = c('gamma', 'beta', 'phi', 'sigma', 'tau'))
plot(mfit, pars = c('gamma', 'phi', 'sigma', 'tau'))
pairs(mfit, pars = c('gamma', 'phi', 'sigma', 'tau'))
# shinystan::launch_shinystan(mfit)



# run on real data: ####
GPP <- dd$GPP
GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1

datlist <- list(
    N = nrow(dd), K = 4,
    S = S, ss = as.numeric(dd$site),
    light = dd$light,
    biomass = dd$fila_chla_mgm2_fit,
    K600 = dd$K600,
    P = GPP, P_sd = GPP_sd,
    new_ts = new_ts
)

dfit <- sampling(
    P_mod,
    data = datlist
)

plot(dfit, pars = c('phi', 'gamma[2]', 'gamma[3]', 'gamma[4]', 'sigma', 'tau'),
     show_density = TRUE)+
    xlim(0, 6)
b <- plot(dfit, pars = c('gamma[1]', 'beta'), show_density = TRUE)
ggpubr::ggarrange(a, b)
print(dfit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma'))
# shinystan::launch_shinystan(dfit)


# look at model residuals
preds <- rstan::summary(dfit, pars = 'y_tilde')$summary %>%
    data.frame() %>%
    select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
preds_poly <- data.frame(date = c(dd$date, rev(dd$date)),
                         site = c(dd$site, rev(dd$site)),
                         y = c(preds$P_mod_lower, rev(preds$P_mod_upper))) %>%
    left_join(dd, by = c('site', 'date'))

preds <- bind_cols(dd, preds)

png('figures/hierarchical_biomass_model_fits.png',
    width = 6, height = 8, units = 'in', res = 300)
    ggplot(preds, aes(date, GPP)) +
        geom_point(size = 0.8) +
        geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), linewidth = 0.3) +
        geom_line(aes(y = P_mod), col = 'steelblue', linewidth = 0.75) +
        geom_polygon(data = preds_poly, aes(date, y),
                     fill = alpha('steelblue', 0.3),
                     col = alpha('steelblue', 0.3))+
        facet_grid(site~year, scales = 'free') +
        theme_bw()
dev.off()
# ggpubr::ggarrange(PL, DL, GR, GC, BM, BN)

# hierarchical model with biomass uncertainty ##################################
# this model isn't working yet:
# it gives an error when trying to initialize the biomass draw variable. Will troubleshoot with Bob.

# Simulate AR1 data ####

sim_dat <- dd
# set params
n <- nrow(sim_dat)
gamma <- c(-10,6,3,0.2)         # population level coefficients
phi <- 0.5                      # ar1 parameter
sigma_proc <- 1
sigma_obs <- 0.5
tau <- 1                        # standard deviation of intercepts
K <- length(gamma)              # number of parameters
S <- length(unique(dd$site))    # number of sites
ss <- as.numeric(dd$site)       # site for each observation

light <- sim_dat$light
biomass <- sim_dat$fila_chla_mgm2_fit
biomass_se <- sim_dat$fila_chla_mgm2_se
K600 <- sim_dat$K600

# select random intercept for each site
beta <- vector(mode = 'double', length = S)
for(i in 1:S){
    beta[i] <- rnorm(1, gamma[1], tau)
}

# simulate response vector
mu <- vector(mode = 'double', length = n)
mu[1] = (beta[ss[1]] + gamma[2]*light[1] + gamma[3]*biomass[1] +
             gamma[4] * K600[1])/(1-phi)

for(i in 2:n){
    mu[i] <- beta[ss[i]] + gamma[2]*light[i] + gamma[3]*biomass[i] +
        gamma[4] * K600[i] + phi*mu[i-1] +
        rnorm(1, 0, sigma_proc)
}

P <- vector(mode = "double", length = n)
P <- rnorm(n, mean = mu, sd = sigma_obs)
P_sd <- rep(sigma_obs, n)
sim <- data.frame(site = ss, light = light, biomass = biomass, GPP = P)
ggplot(sim, aes(light, GPP, col = factor(site))) +
    geom_point()
new_ts <- rep(0, n)
new_ts[rle2(paste0(sim_dat$year, sim_dat$site))$starts] <- 1
# compile data for stan

datlist <- list(
    N = n, K = K,
    S = S, ss = ss,
    light = light,
    biomass = biomass,
    biomass_se = biomass_se,
    K600 = K600,
    P = P, P_sd = P_sd,
    new_ts = new_ts
)

# compile stan model
P_mod <- stan_model("code/model/stan_code/GPP_biomass_model_ar1_hierarchical_uncertainty.stan")

# fit the model
mfit_u <- sampling(
    P_mod,
    data = datlist, chains = 1
)

print(mfit_u, pars = c('gamma', 'beta', 'phi', 'sigma', 'tau'))
plot(mfit, pars = c('gamma', 'phi', 'sigma', 'tau'))
pairs(mfit, pars = c('gamma', 'phi', 'sigma', 'tau'))
# shinystan::launch_shinystan(mfit)



# run on real data: ####
GPP <- dd$GPP
GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1

datlist <- list(
    N = nrow(dd), K = 4,
    S = S, ss = as.numeric(dd$site),
    light = dd$light,
    biomass = dd$fila_chla_mgm2_fit,
    K600 = dd$K600,
    P = GPP, P_sd = GPP_sd,
    new_ts = new_ts
)

dfit <- sampling(
    P_mod,
    data = datlist
)

a <- plot(dfit, pars = c('phi', 'gamma', 'sigma'))+
    xlim(0, 7)
b <- plot(dfit, pars = c('beta', 'tau'))
ggpubr::ggarrange(a, b)
print(dfit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma'))
# shinystan::launch_shinystan(dfit)


# look at model residuals
preds <- rstan::summary(dfit, pars = 'y_tilde')$summary %>%
    data.frame() %>%
    select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
preds_poly <- data.frame(date = c(dd$date, rev(dd$date)),
                         site = c(dd$site, rev(dd$site)),
                         y = c(preds$P_mod_lower, rev(preds$P_mod_upper))) %>%
    left_join(dd, by = c('site', 'date'))

preds <- bind_cols(dd, preds)


ggplot(preds, aes(date, GPP)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper)) +
    geom_line(aes(y = P_mod), col = 'steelblue', linewidth = 1) +
    geom_polygon(data = preds_poly, aes(date, y),
                 fill = alpha('steelblue', 0.3),
                 col = alpha('steelblue', 0.3))+
    facet_grid(site~year, scales = 'free') +
    theme_bw()

# ggpubr::ggarrange(PL, DL, GR, GC, BM, BN)
