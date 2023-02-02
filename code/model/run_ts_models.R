# Script to run time series models on UCFR data

library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
library(tidyverse)

# read in helper functions
source('code/UCFR_helpers.R')

# prep data ####
dd <- read_csv('data/model_fits/biomass_metab_model_data.csv')

dd <- dd %>%
    mutate(year = lubridate::year(date),
           across(starts_with(c('epil', 'fila')), ~ log(.+ 0.1),
                  .names = 'log_{.col}'),
           across(starts_with(c('GPP')), ~ log(.), .names = 'log_{.col}'),
           GPP_sd = (GPP.upper - GPP.lower)/3.92,
           log_GPP_sd = (log_GPP.upper - log_GPP.lower)/3.92) %>%
    mutate(across(contains(c('epil','fila', 'light')), ~ scale(.)[,1])) %>%
    group_by(site, year) %>%
           mutate(
               across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    ungroup()

sites <- unique(dd$site)

# compile time series models
lmod <- stan_model(file = "code/model/stan_code/linear_model.stan",
                   model_name = 'lmod')
ar1_lmod <- stan_model("code/model/stan_code/AR1_linear_model.stan",
                       model_name = 'ar1_lmod')
quab_PIcurve <- stan_model("code/model/stan_code/PIcurve_model.stan")
PIcurve <- stan_model("code/model/stan_code/PIcurve_model_3.stan")
ar1_PIcurve <- stan_model("code/model/stan_code/PIcurve_model_2.stan")

# run on real data: ####
X <- matrix(c(dd$light,
              dd$log_epil_gm2_fit, dd$log_fila_gm2_fit,
              dd$log_epil_chla_mgm2_fit, dd$log_fila_chla_mgm2_fit),
            ncol = 5)
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1

datlist <- list(
    N = nrow(dd), K = 1,
    S = length(sites), ss = as.numeric(dd$site),
    X = X[,1, drop = F], P = dd$log_GPP, P_sd = dd$log_GPP_sd
    # new_ts = new_ts
)

fit <- sampling(
    lmod,
    data = datlist
)

preds <- get_model_preds(fit, dd)
plot_model_fit(preds)

# Compile set of linear model runs:
master_X <- select(dd, light, epil_afdm = log_epil_gm2_fit,
                   fila_afdm = log_fila_gm2_fit,
                   epil_chla = log_epil_chla_mgm2_fit,
                   fila_chla = log_fila_chla_mgm2_fit) %>%
    mutate(across(-light, .fn = ~light * ., .names = 'light_{.col}'))

model_combinations <- list(1, c(1,2), c(1,3), c(1,2,3),
                           c(1,4), c(1,5), c(1,4,5),
                           c(1,2,6), c(1,3,7), c(1,2,3,6,7),
                           c(1,4,8), c(1,5,9), c(1,4,5,8,9),
                           6, 7, 8, 9, c(6,7), c(8,9))
mod_ests <- data.frame()
chains <- data.frame()

for(i in 1:length(model_combinations)){
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(lmod, dd, X = as.matrix(X))
    preds <- get_model_preds(fit, dd)
    plot_model_fit(preds) +
        ggtitle(paste0('lmod: ', paste(colnames(X), collapse = ' + ')))
    ests <- extract_model_params(fit, X, dd, extract.chains = TRUE)

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
}

write_csv(mod_ests, 'data/model_fits/linear_model_parameter_ests.csv')
write_csv(chains, 'data/model_fits/linear_model_chains.csv')

mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)

mod_ests <- data.frame()
chains <- data.frame()
for(i in 1:length(model_combinations)){
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(ar1_lmod, dd, X = as.matrix(X))
    preds <- get_model_preds(fit, dd)
    p <- plot_model_fit(preds) +
        ggtitle(paste0('ar1_lmod: ', paste(colnames(X), collapse = ' + ')))
    print(p)
    ests <- extract_model_params(fit, X, dd, phi_corrected = T,
                                 extract.chains = TRUE)
    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
}

write_csv(mod_ests, 'data/model_fits/ar1_linear_model_parameter_ests.csv')
write_csv(chains, 'data/model_fits/ar1_linear_model_chains.csv')

mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)



model with no biomass
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
