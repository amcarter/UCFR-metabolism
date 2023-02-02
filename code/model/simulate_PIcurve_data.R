#######################################################################
# Script to simulate data to test different PI curve models of
# GPP as a function of light and biomass
#######################################################################

library(rstan)
library(tidyverse)

# prep data ####
dd <- read_csv('data/model_fits/biomass_metab_model_data.csv')

dd <- dd %>%
    mutate(year = lubridate::year(date),
           # across(starts_with(c('epil', 'fila')), ~ log(.+ 0.1),
           #        .names = 'log_{.col}'),
           # across(starts_with(c('GPP')), ~ log(.), .names = 'log_{.col}'),
           # log_GPP_sd = (log_GPP.upper - log_GPP.lower)/3.92,
           GPP_sd = (GPP.upper - GPP.lower)/3.92) %>%
    # mutate(across(contains(c('epil','fila', 'light')), ~ scale(.)[,1])) %>%
    group_by(site, year) %>%
    mutate(
        across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
        site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    ungroup()

sites <- unique(dd$site)

# Simulate data for a basic PI curve model
N <- nrow(dd)
S <- length(sites)
ss <- as.numeric(dd$site)
biomass <- dd$fila_chla_mgm2_fit
light <- dd$light

#set parameters
Pmax <- 1
alpha <- 1
tau_a <- 0.5
tau_P <- 1
Pmax_s <- rnorm(S, Pmax, tau_P)
alpha_s <- rnorm(S, alpha, tau_a)
sigma <- 1

mu <- Pmax_s[ss] * biomass * tanh((alpha_s[ss] * light)/Pmax_s[ss])
dd$P <- rnorm(N, mu, sigma)

ggplot(dd, aes(light, P, col = fila_chla_mgm2_fit))+
    geom_point()+
    geom_point(aes(y = GPP), col = 'brown')+
    facet_grid(site~year)

PI_curve <- stan_model("code/model/stan_code/PIcurve_model.stan",
                       model_name = 'PI_curve')

datlist <- list(
    N = N, biomass = biomass,
    light = light, S = S, ss = ss, P = dd$P)#, P_sd = rep(P_sd, N)

fit <- sampling(
    PI_curve,
    datlist
)

p1 <- plot(fit, pars = c("Pmax", "alpha", "tau_a", "tau_Pmax", "sigma"))

datlist <- list(
    N = N, biomass = biomass,
    light = light, S = S, ss = ss, P = dd$GPP
)

fit_dat <- sampling(
    PI_curve,
    datlist
)
p2 <- plot(fit_dat, pars = c("Pmax", "alpha", "tau_a", "tau_Pmax", "sigma"))

preds <- get_model_preds(fit_dat, dd)
plot(preds$P_mod, preds$GPP)
plot(preds$GPP, preds$P_mod)
abline(0,1)
1-sum((preds$P_mod - preds$GPP)^2)/sum((preds$GPP - mean(preds$GPP))^2)
1-sum((preds$GPP - preds$P_mod)^2)/sum((preds$GPP - mean(preds$GPP))^2)
calculate_r2_adj(fit_dat, dd, log = F)

