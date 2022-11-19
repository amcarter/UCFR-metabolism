# compare a basic linear model of GPP to one with an AR1 term:

library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)

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


# insert fake data test here

## run on datasets ####
BN <- filter(dd, site == 'BN') %>%
    mutate(across(c(starts_with('GPP'), 'light', starts_with('epil'),
                           starts_with('fila')),
                  zoo::na.approx, na.rm = F)) %>%
    filter(!is.na(GPP) & !is.na(epil_gm2_fit) & !is.na(fila_gm2_fit))

GPP <- BN$GPP
GPP_sd <- (BN$GPP.upper - BN$GPP.lower)/3.92
# compile stan models
l_mod <- stan_model("code/model/GPP_biomass_model_one_site.stan")
ar1_mod <- stan_model("code/model/GPP_biomass_model_ar1.stan")

mdat <- list(N = nrow(BN),
               K = 3,
               light = BN$light,
               biomass = BN$fila_chla_mgm2_fit,
               P = GPP,
               P_sd = GPP_sd)

l_fit <- sampling(l_mod, data = mdat)
plot(l_fit, pars = 'gamma')


ar1_fit <- sampling(ar1_mod, data = mdat)#, chains = 1)
# ar1_fit <- sampling(ar1_mod, data = mdat, chains = 1)
print(ar1_fit, pars = c('phi', 'sigma', 'gamma'))
plot(ar1_fit, pars = c('phi', 'sigma', 'gamma'))
pairs(ar1_fit, pars = c('phi', 'sigma', 'gamma'))
