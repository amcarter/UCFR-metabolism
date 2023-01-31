# Script to run time series models on UCFR data

library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
library(tidyverse)

# read in helper functions
source('code/UCFR_helpers.R')

# prep data ####
dd <- read_csv('data/model_fits/biomass_metab_model_data.csv')

dd<- dd %>%
    mutate(across(starts_with(c('epil', 'fila')), ~ .+ 0.1)) %>%
    mutate(across(starts_with(c('GPP', 'epil', 'fila')), ~ log(.))) %>%
    mutate(across(starts_with(c('epil','fila', 'light')), ~scale(.)[,1])) %>%
    group_by(site, year) %>%
           mutate(
               across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
           year = lubridate::year(date),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light))

sites <- unique(dd$site)

# compile time series models
lmod <- stan_model("co")
ar1_lmod <- stan_model("code/model/stan_code/AR1_linear_model.stan")
ar1_mod_nb <- stan_model("code/model/stan_code/GPP_biomass_model_ar1_noK600_hierarchical.stan")
# ar1_mod_nb <- stan_model("code/model/stan_code/GPP_NObiomass_model_ar1_noK600_hierarchical.stan")
ar1_lmod_int <- stan_model("")
ar1_int_only <- stan_model("")
PIcurve <- stan_model("")
ar1_PIcurve <- stan_model("")

# run on real data: ####
biomass <- matrix(dd$fila_chla_mgm2_fit, ncol = 1)
X <- matrix(c(dd$light, dd$fila_chla_mgm2_fit), ncol = 2)
dd$GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1

datlist <- list(
    N = nrow(dd), K = ncol(X),
    S = length(sites), ss = as.numeric(dd$site),
    X = X, P = dd$GPP, P_sd = dd$GPP_sd,
    new_ts = new_ts
)


fit <- sampling(
    ar1_lmod,
    data = datlist
)

datlist2 <- list(
    N = nrow(dd), K = ncol(biomass),
    S = length(sites), ss = as.numeric(dd$site),
    biomass = biomass, light = dd$light,
    P = dd$GPP, P_sd = dd$GPP_sd,
    new_ts = new_ts
)


fit2 <- sampling(
    ar1_mod_nb,
    data = datlist2
)

p1 <- plot(fit, pars = c('phi', 'gamma', 'tau', 'beta'))
p2 <- plot(fit2, pars = c('phi', 'gamma', 'tau', 'beta'))
p3 <- plot(fit, pars = c('phi', 'gamma', 'tau', 'beta'))

ggpubr::ggarrange(p1,p2, p3)
#
# print(dfit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma'))
# plot(dfit, pars = c('phi', 'gamma[2]', 'gamma[3]', 'gamma[4]', 'sigma', 'tau'),
#      show_density = TRUE)
# plot(dfit, pars = c('gamma[1]', 'beta'), show_density = TRUE)
#
# # shinystan::launch_shinystan(dfit)
# preds <- get_model_preds(dfit, dd)
# rstan::loo(dfit)

# run different model combinations: ####
mod_ests <- data.frame()
chains <- data.frame()
# Chlorophyll:
# fila
biomass <- matrix(dd$fila_chla_mgm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod, K600 = FALSE)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_chla')
ests <- extract_model_params(fit, 'fila', 'chla_mgm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains(fit, 'fila', 'chla_mgm2')
chains <- bind_rows(chains, cc)

# epil
biomass <- matrix(dd$epil_chla_mgm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod, K600 = FALSE)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'epil_chla')
ests <- extract_model_params(fit, 'epil', 'chla_mgm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains(fit, 'epil', 'chla_mgm2')
chains <- bind_rows(chains, cc)

# fila + epil
biomass <- matrix(c(dd$fila_chla_mgm2_fit, dd$epil_chla_mgm2_fit), ncol = 2)
fit <- fit_biomass_model(dd, biomass, P_mod, K600 = FALSE)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_epil_chla')
ests <- extract_model_params(fit, 'fila_epil', 'chla_mgm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains(fit, 'fila_epil', 'chla_mgm2')
chains <- bind_rows(chains, cc)

# biomass:
# fila
biomass <- matrix(dd$fila_gm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod, K600 = FALSE)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_gm2')
ests <- extract_model_params(fit, 'fila', 'gm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains(fit, 'fila', 'gm2')
chains <- bind_rows(chains, cc)

# epil
biomass <- matrix(dd$epil_gm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod, K600 = FALSE)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'epil_gm2')
ests <- extract_model_params(fit, 'epil', 'gm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains(fit, 'epil', 'gm2')
chains <- bind_rows(chains, cc)

# fila + epil
biomass <- matrix(c(dd$fila_gm2_fit, dd$epil_gm2_fit), ncol = 2)
fit <- fit_biomass_model(dd, biomass, P_mod, K600 = FALSE)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_epil_gm2')
ests <- extract_model_params(fit, 'fila_epil', 'gm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains(fit, 'fila_epil', 'gm2')
chains <- bind_rows(chains, cc)


# model with no biomass
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1
dd$GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92
datlist <- list(
    N = nrow(dd),
    S = length(sites), ss = as.numeric(dd$site),
    # K600 = dd$K600,
    light = dd$light,
    P = dd$GPP, P_sd = dd$GPP_sd,
    new_ts = new_ts
)

fit <- sampling(
    P_mod_nb,
    data = datlist
)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'no_biomass')
ests <- extract_model_params(fit, NA, NA, phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains(fit, NA, NA)
chains <- bind_rows(chains, cc)

# compare parameter estimates across models

mod_ests <- mod_ests %>%
    mutate(parameter = sub('\\.{3}[0-9]+$', '', mod_ests$parameter),
           parameter = str_replace_all(parameter, '[\\[\\]]', ''))
mod_ests <- mod_ests %>%
    relocate(biomass_vars, units, parameter) %>%
    mutate(parameter = case_when(parameter == 'gamma1' ~ 'intercept',
                                 parameter == 'gamma2' ~ 'gamma_light',
                                 parameter == 'gamma3' &
                                     biomass_vars %in% c('fila', 'fila_epil') ~ 'gamma_fila',
                                 parameter == 'gamma3' &
                                     biomass_vars == 'epil' ~ 'gamma_epil',
                                 parameter == 'gamma4' ~ 'gamma_epil',
                                 TRUE ~ parameter))
beepr::beep(5)

write_csv(mod_ests, 'data/model_fits/hmodel_GPP_bm_linGAM_pars.csv')
write_csv(chains, 'data/model_fits/hmodel_GPP_bm_linGAM_post_dists_for_plot.csv')


mod_ests <- read_csv( 'data/model_fits/hmodel_logGPP_bm_logGAM_pars.csv')
chains <- read_csv('data/model_fits/hmodel_logGPP_bm_logGAM_post_dists_for_plot.csv')

mod_ests %>%
    filter(!grepl('^beta', parameter) ,
           parameter != 'intercept') %>%
    ggplot( aes(parameter, mean, fill = units, col = units))+
    geom_boxplot(aes(ymin = X2.5., lower = X25., middle = X50.,
                     upper = X75., ymax = X97.5.)) +
    # ylim(0,2) +
    facet_wrap(.~biomass_vars)


mod_ests %>% filter(parameter %in% c('gamma_epil', 'gamma_fila'))

plot(dd$epil_gm2_fit, dd$epil_chla_mgm2_fit)
abline(0,1)
plot(dd$fila_gm2_fit, dd$fila_chla_mgm2_fit)


library(viridis)
chains <- data.frame(parameter = c(rep('gamma_fila',2),
                                   rep('gamma_epil',2)),
                     chains = rep(0,4),
                     biomass_vars = rep(NA_character_,4),
                     units = rep(NA_character_,4)) %>%
    bind_rows(chains)

chains %>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = factor(parameter,
                              levels=c("phi", "gamma_light", "gamma_fila",
                                       "gamma_epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ylab("Parameter Estimate") +
    xlab("") +
    theme_bw()
