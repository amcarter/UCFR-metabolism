# Script to run time series models on UCFR data

library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
library(tidyverse)
library(viridis)

# read in helper functions
source('code/UCFR_helpers.R')

# prep data ####
dat <- read_csv('data/model_fits/biomass_metab_model_data.csv')

dat <- mutate(dat,
              biomass = epil_chla_mgm2_fit + fila_chla_mgm2_fit,
              log_biomass = log(epil_chla_mgm2_fit + 0.1) +
                  log(fila_chla_mgm2_fit + 0.1))
q1 <- dat %>% select(starts_with(c('epil', 'fila')), -ends_with('_se')) %>%
    mutate(across(.fns = ~case_when(. == 0 ~ NA_real_,
                                    TRUE ~ .))) %>%
    summarize(across(.fns =  ~quantile(., 0.01, na.rm = T)))
dat <- dat %>%
    mutate(year = lubridate::year(date),
           month = lubridate::month(date),
           epil_gm2_fit = epil_gm2_fit + q1$epil_gm2_fit,
           epil_chla_mgm2_fit = epil_chla_mgm2_fit + q1$epil_chla_mgm2_fit,
           fila_gm2_fit = fila_gm2_fit + q1$fila_gm2_fit,
           fila_chla_mgm2_fit = fila_chla_mgm2_fit + q1$fila_chla_mgm2_fit,
           across(starts_with(c('epil', 'fila', 'light', 'GPP')), ~ log(.),
                  .names = 'log_{.col}'),
           GPP_sd = (GPP.upper - GPP.lower)/3.92,
           log_GPP_sd = (log_GPP.upper - log_GPP.lower)/3.92)
           # light = light/max(light),
           # log_light = log_light - min(log_light),
           # log_light = log_light/max(log_light)) %>%
dd <- dat %>%
    mutate(across(contains(c('epil', 'fila', 'light')), ~ scale(.)[,1])) %>%
    group_by(site, year) %>%
           mutate(
               across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_chla_mgm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    ungroup()

sites <- unique(dd$site)

# compile time series model
ar1_lmod_ss <- stan_model("code/model/stan_code/AR1_linear_model_ss.stan",
                       model_name = 'ar1_lmod_ss')

# run on real data: ####
master_X <- select(dd, log_light,
                   epil_afdm = log_epil_gm2_fit,
                   fila_afdm = log_fila_gm2_fit,
                   epil_chla = log_epil_chla_mgm2_fit,
                   fila_chla = log_fila_chla_mgm2_fit) %>%
    mutate(across(-log_light, .fn = ~log_light * .,
                  .names = 'log_light_{.col}'))

colnames(master_X)

model_combinations <- list(1, c(1,2), c(1,3), c(1,4), c(1,5),
                           c(1,2,3), c(1,4,5),
                           c(1,2,6), c(1,3,7),
                           c(1,4,8), c(1,5,9),
                           c(1,2,3,6,7), c(1,4,5,8,9),
                           6,7,8,9, c(6,7), c(8,9))

# fit AR1 linear models ####
mod_ests <- data.frame()
chains <- data.frame()
mod_comp <- data.frame()

for(i in 1:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(ar1_lmod_ss, dd, X = as.matrix(X))
    preds <- get_model_preds(fit, dd)
    p <- plot_model_fit(preds) +
        ggtitle(paste0('ar1_lmod: ', paste(colnames(X), collapse = ' + ')))
    print(p)
    ggsave(paste('figures/biomass_models/unconditioned_preds/ar1_lmod_ss',
              paste(colnames(X),collapse = '_'),'2.png'), p,
        width = 5, height = 7, units = 'in')

    ests <- extract_model_params(fit, X, dd, phi_corrected = T,
                                 extract.chains = TRUE)

    preds <- mutate(preds, P_delta = P_mod - ests$ests$mean[1] *
                   c(preds$P_mod[1], preds$P_mod[1:(nrow(preds)-1)]))

    comp <- calculate_model_metrics(fit, X, dd, 'ar1_lmod')

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
    mod_comp <- bind_rows(mod_comp, comp)

}
beepr::beep(5)

write_csv(mod_ests, 'data/model_fits/ar1_linear_model_parameter_ests_ss.csv')
write_csv(chains, 'data/model_fits/ar1_linear_model_chains_ss.csv')
write_csv(mod_comp, 'data/model_fits/ar1_linear_model_metrics_ss_cond.csv')

mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)

mod_ests %>% filter(biomass_vars == 'log_light + epil_chla + fila_chla') %>%
    select(parameter, mean)
mod_comp %>% arrange(waic)
plot(mod_comp$waic, mod_comp$rmse)


# fit PI curve models ####

master_X <- select(dd, light, epil_afdm = epil_gm2_fit,
                   fila_afdm = fila_gm2_fit,
                   epil_chla = epil_chla_mgm2_fit,
                   fila_chla = fila_chla_mgm2_fit) %>%
    mutate(across(.fns  = ~. - min(.)))
model_combinations <- list(2, 3, 4, 5, c(2,3), c(4,5))

mod_ests <- data.frame()
chains <- data.frame()
mod_comp <- data.frame()
for(i in 1:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(PIcurve, dd, X = as.matrix(X), log_GPP = FALSE)
    preds <- get_model_preds(fit, dd)
    plot_model_fit(preds, log_ests = FALSE) +
        ggtitle(paste0('PI_curve: ', paste(colnames(X), collapse = ' + ')))
    ests <- extract_model_params(fit, X, dd, extract.chains = TRUE)
    comp <- calculate_model_metrics(fit, X, dd, 'PI_curve', log = FALSE)

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
    mod_comp <- bind_rows(mod_comp, comp)

}
 beepr::beep(5)

write_csv(mod_ests, 'data/model_fits/PIcurve_parameter_ests.csv')
write_csv(chains, 'data/model_fits/PIcurve_chains.csv')
write_csv(mod_comp, 'data/model_fits/PIcurve_model_metrics.csv')

# fit AR PI curve models ####

mod_ests <- data.frame()
chains <- data.frame()
mod_comp <- data.frame()

for(i in 2:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(ar1_PIcurve_ss, dd, X = as.matrix(X), log_GPP = FALSE)
    preds <- get_model_preds(fit, dd)
    p <- plot_model_fit(preds, log_ests = FALSE) +
        ggtitle(paste0('AR1_PI_curve:', paste(colnames(X), collapse = ' + ')))
    print(p)
    ggsave(paste('figures/biomass_models/conditioned_preds/ar1_PIcurve_ss',
              paste(colnames(X),collapse = '_'),'.png'), p,
        width = 5, height = 7, units = 'in')

    ests <- extract_model_params(fit, X, dd, phi_corrected = TRUE,
                                 extract.chains = TRUE)
    comp <- calculate_model_metrics(fit, X, dd, 'ar1_PI_curve', log = FALSE)

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
    mod_comp <- bind_rows(mod_comp, comp)

}
beepr::beep(5)


write_csv(mod_ests, 'data/model_fits/ar_PIcurve_parameter_ests.csv')
write_csv(chains, 'data/model_fits/ar_PIcurve_chains.csv')
write_csv(mod_comp, 'data/model_fits/ar_PIcurve_model_metrics_conditioned.csv')

####

