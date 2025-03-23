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
           month = lubridate::month(date),
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

dd <- dd %>%
    select(site, date, GPP, log_GPP, GPP_sd, log_GPP_sd, light,
           epil_gm2_fit, epil_chla_mgm2_fit,
           fila_gm2_fit, fila_chla_mgm2_fit,
           log_epil_gm2_fit, log_epil_chla_mgm2_fit,
           log_fila_gm2_fit, log_fila_chla_mgm2_fit)

log_GPP <- dd %>% select(site, date, log_GPP) %>%
    pivot_wider(names_from = 'site', values_from = 'log_GPP')

log_GPP_sd <- dd %>% select(site, date, log_GPP_sd) %>%
    pivot_wider(names_from = 'site', values_from = 'log_GPP_sd') %>%
    select(-date)


start_idx <- unname(apply(log_GPP[,-1], 2,
                          function(x) Position(function(y) ! is.na(y), x)))
end_idx <- unname(apply(log_GPP[,-1], 2,
                        function(x) Position(function(y) ! is.na(y), x,
                                             right = TRUE)))
new_year <- Position(function(y) y > lubridate::ymd('2021-01-01'), log_GPP$date)

S <- ncol(log_GPP)-1
N <- nrow(log_GPP)
K = 9


master_X <- array(rep(NA, N * K * S), dim = c(S, N,K))

for(i in 1:length(sites)){
    ss = sites[i]
    smax <- dd %>%
        filter(site == ss) %>%
        select(date, light, epil_afdm = log_epil_gm2_fit,
               fila_afdm = log_fila_gm2_fit,
               epil_chla = log_epil_chla_mgm2_fit,
               fila_chla = log_fila_chla_mgm2_fit) %>%
        right_join(log_GPP[,1], by = 'date') %>%
        mutate(across(-date, zoo::na.approx, na.rm = F),
               across(-date, .fn = ~case_when(is.na(.) ~ mean(., na.rm = T),
                                              TRUE ~ .)),
               across(-c(light, date), .fn = ~light * ., .names = 'light_{.col}')) %>%
        select(-date)
    master_X[i,,] <- as.matrix(smax)
}

log_GPP <- log_GPP %>%
    select(-date) %>%
    mutate(across(.fn = ~case_when(is.na(.) ~ mean(., na.rm = T),
                                              TRUE ~ .)))

log_GPP_sd <- log_GPP_sd %>%
    mutate(across(.fn = ~case_when(is.na(.) ~ mean(., na.rm = T),
                                              TRUE ~ .)))


ar1_lmod_ss <- stan_model("code/model/stan_code/AR1_linear_model_ss_matrix.stan",
                       model_name = 'ar1_lmod_ss')

datlist <- list(
    S = S, N = N, K = 1,
    X = master_X[,,1, drop = F],
    P = as.matrix(log_GPP),
    P_sd = as.matrix(log_GPP_sd),
    start_idx = start_idx,
    end_idx = end_idx,
    new_year = new_year
)

fit2 <- sampling(
    ar1_lmod_ss,
    datlist,
    chains = 4,
    cores = 4
)

# compile time series models
lmod <- stan_model(file = "code/model/stan_code/linear_model.stan",
                   model_name = 'lmod')
ar1_lmod <- stan_model("code/model/stan_code/AR1_linear_model.stan",
                       model_name = 'ar1_lmod')
# quab_PIcurve <- stan_model("code/model/stan_code/PIcurve_model.stan")
PIcurve <- stan_model("code/model/stan_code/PIcurve_model_3.stan",
                      model_name = 'PI_curve')
ar1_PIcurve <- stan_model("code/model/stan_code/ar_PIcurve_model.stan",
                      model_name = 'PI_curve')

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

preds <- get_model_preds(fit, dd_train)
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

# fit linear models ####
mod_ests <- data.frame()
chains <- data.frame()
mod_comp <- data.frame()

for(i in 1:length(model_combinations)){
    print(paste('model', i, ' of ', length(model_combinations)))
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(lmod, dd, X = as.matrix(X))
    preds <- get_model_preds(fit, dd)
    plot_model_fit(preds) +
        ggtitle(paste0('lmod: ', paste(colnames(X), collapse = ' + ')))
    ests <- extract_model_params(fit, X, dd, extract.chains = TRUE)
    comp <- calculate_model_metrics(fit, X, dd, 'lmod')

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
    mod_comp <- bind_rows(mod_comp, comp)
}

write_csv(mod_ests, 'data/model_fits/linear_model_parameter_ests.csv')
write_csv(chains, 'data/model_fits/linear_model_chains.csv')
write_csv(mod_comp, 'data/model_fits/linear_model_metrics.csv')

mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)
mod_comp %>% arrange(waic)
plot(mod_comp$waic, mod_comp$rmse)

# fit AR1 linear models ####
mod_ests <- data.frame()
chains <- data.frame()
mod_comp <- data.frame()

for(i in 1:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(ar1_lmod_ss, dd, X = as.matrix(X), ss = T)
    preds <- get_model_preds(fit, dd)
    plot_model_fit(preds) +
        ggtitle(paste0('ar1_lmod: ', paste(colnames(X), collapse = ' + ')))
    ests <- extract_model_params(fit, X, dd, phi_corrected = T,
                                 extract.chains = TRUE)

    comp <- calculate_model_metrics(fit, X, dd, 'ar1_lmod')

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
    mod_comp <- bind_rows(mod_comp, comp)

}
beepr::beep(5)

write_csv(mod_ests, 'data/model_fits/ar1_linear_model_parameter_ests.csv')
write_csv(chains, 'data/model_fits/ar1_linear_model_chains.csv')
write_csv(mod_comp, 'data/model_fits/ar1_linear_model_metrics.csv')

mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)
mod_comp %>% arrange(waic)
plot(mod_comp$waic, mod_comp$rmse)


# fit PI curve models ####

master_X <- select(dd, light, epil_afdm = epil_gm2_fit,
                   fila_afdm = fila_gm2_fit,
                   epil_chla = epil_chla_mgm2_fit,
                   fila_chla = fila_chla_mgm2_fit)
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
 # beepr::beep(5)

write_csv(mod_ests, 'data/model_fits/PIcurve_parameter_ests.csv')
write_csv(chains, 'data/model_fits/PIcurve_chains.csv')
write_csv(mod_comp, 'data/model_fits/PIcurve_model_metrics.csv')

# fit AR PI curve models ####

mod_ests <- data.frame()
chains <- data.frame()
mod_comp <- data.frame()

for(i in 1:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(ar1_PIcurve, dd, X = as.matrix(X), log_GPP = FALSE)
    preds <- get_model_preds(fit, dd)
    plot_model_fit(preds, log_ests = FALSE) +
        ggtitle(paste0('AR1_PI_curve:', paste(colnames(X), collapse = ' + ')))
    ests <- extract_model_params(fit, X, dd, phi_corrected = TRUE,
                                 extract.chains = TRUE)
    comp <- calculate_model_metrics(fit, X, dd, 'ar1_PI_curve', log = FALSE)

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
    mod_comp <- bind_rows(mod_comp, comp)

}
# beepr::beep(5)


write_csv(mod_ests, 'data/model_fits/ar_PIcurve_parameter_ests2.csv')
write_csv(chains, 'data/model_fits/ar_PIcurve_chains2.csv')
write_csv(mod_comp, 'data/model_fits/ar_PIcurve_model_metrics2.csv')

####

