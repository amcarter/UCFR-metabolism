library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
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

# functions to interact with model output: ####
fit_biomass_model <- function(dd, biomass, model, K600 = TRUE){

    new_ts <- rep(0, nrow(dd))
    new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1
    dd$GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92

    datlist <- list(
        N = nrow(dd), K = ncol(biomass),
        S = length(sites), ss = as.numeric(dd$site),
        K600 = dd$K600,
        light = dd$light,
        biomass = biomass,
        P = dd$GPP,
        P_sd = dd$GPP_sd,
        new_ts = new_ts
    )
    if(!K600){
    datlist <- list(
        N = nrow(dd), K = ncol(biomass),
        S = length(sites), ss = as.numeric(dd$site),
        light = dd$light,
        biomass = biomass,
        P = dd$GPP,
        P_sd = dd$GPP_sd,
        new_ts = new_ts
    )

    }

    fit <- sampling(
        model,
        data = datlist
    )

    return(fit)
}

get_model_preds <- function(fit, dat){
    preds <- rstan::summary(fit, pars = 'y_tilde')$summary %>%
        data.frame() %>%
        select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
    preds <- bind_cols(dat, preds)

    return(preds)
}

plot_model_fit <- function(preds, dat, mod = 'fila_epil_chla'){
    preds_poly <- data.frame(date = c(dat$date, rev(dat$date)),
                             site = c(dat$site, rev(dat$site)),
                             y = c(preds$P_mod_lower, rev(preds$P_mod_upper))) %>%
        left_join(dat, by = c('site', 'date'))

     p <- ggplot(preds, aes(date, GPP)) +
            geom_point(size = 0.8) +
            geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), linewidth = 0.3) +
            geom_line(aes(y = P_mod), col = 'steelblue', linewidth = 0.75) +
            geom_polygon(data = preds_poly, aes(date, y),
                         fill = alpha('steelblue', 0.3),
                         col = alpha('steelblue', 0.3))+
            facet_grid(site~year, scales = 'free') +
            ylab(expression(paste('GPP (g ', O[2], ' ', m^-2, d^-1, ')')))+
            xlab('')+
            theme_bw()

    ggsave(paste0('figures/biomass_models/hierarchical_biomass_model_fit_',mod,'.png'),
           plot = p, width = 6, height = 7, units = 'in', dpi = 300)
}

extract_model_params <- function(fit, bmv = 'fila + epil', units = 'gm2',
                                 phi_corrected = FALSE){
    ss <- summary(fit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma'))$summary %>%
        data.frame()
    ss$biomass_vars = bmv
    ss$units = units
    ss$parameter = row.names(ss)
    row.names(ss) <- NULL

    w <- which(ss$parameter != 'phi')

    if(phi_corrected){
    ss <- mutate(ss,
                 across(c('mean', 'se_mean', 'sd', starts_with('X')),
                        ~case_when(parameter != 'phi' ~ ./(1-ss$mean[1]),
                                   TRUE ~ .)))
    }

    return(ss)
}

extract_model_chains <- function(fit, bmv = 'fila + epil', units = 'gm2'){
    ss <- rstan::extract(fit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma')) %>%
        data.frame()
    ss <- ss %>%
        select(phi, starts_with('gamma')) %>%
        select(-gamma.1) %>%
        mutate(across(starts_with('gamma'),
                      ~ ./(1-phi))) %>%
        pivot_longer(cols = everything(), names_to = 'parameter',
                     values_to = 'chains') %>%
        mutate(biomass_vars = bmv,
               units = units,
               parameter = case_when(parameter == 'gamma.2' ~ 'gamma_light',
                                     parameter == 'gamma.3' &
                                        biomass_vars %in% c('fila', 'fila_epil') ~ 'gamma_fila',
                                     parameter == 'gamma.3' &
                                        biomass_vars == 'epil' ~ 'gamma_epil',
                                     parameter == 'gamma.4' ~ 'gamma_epil',
                                     TRUE ~ parameter))


    return(ss)
}
extract_model_chains_K600 <- function(fit, bmv = 'fila + epil', units = 'gm2'){
    ss <- rstan::extract(fit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma')) %>%
        data.frame()
    ss <- ss %>%
        select(phi, starts_with('gamma')) %>%
        select(-gamma.1) %>%
        mutate(across(starts_with('gamma'),
                      ~ ./(1-phi))) %>%
        pivot_longer(cols = everything(), names_to = 'parameter',
                     values_to = 'chains') %>%
        mutate(biomass_vars = bmv,
               units = units,
               parameter = case_when(parameter == 'gamma.2' ~ 'gamma_k600',
                                     parameter == 'gamma.3' ~ 'gamma_light',
                                     parameter == 'gamma.4' &
                                        biomass_vars %in% c('fila', 'fila_epil') ~ 'gamma_fila',
                                     parameter == 'gamma.4' &
                                        biomass_vars == 'epil' ~ 'gamma_epil',
                                     parameter == 'gamma.5' ~ 'gamma_epil',
                                     TRUE ~ parameter))


    return(ss)
}
calculate_r2_adj <- function(preds, npar = 13){

    N = nrow(preds)
    r2 = 1 - sum((preds$GPP - preds$P_mod)^2)/
        sum((preds$GPP - mean(preds$GPP))^2)

    r2_adj = 1 - ((1-r2)*(N-1))/(N-npar-1)

    return(r2_adj)
}
calculate_rmse <- function(preds){

    N = nrow(preds)
    rmse = sqrt((1/N) * sum((preds$GPP - preds$P_mod)^2))

    return(rmse)
}


# prep data ####
dd <- read_csv('data/model_fits/biomass_metab_model_data.csv')
par(mfrow = c(3,3))
hist(dd$GPP)
hist(dd$ER)
hist(dd$fila_chla_mgm2_fit)
hist(dd$fila_gm2_fit)
hist(dd$epil_chla_mgm2_fit)
hist(dd$epil_gm2_fit)
hist(dd$light)
hist(dd$K600)
hist(log(dd$q.cms))

dev.off()

# manually scale so that the se's are correct
dd<- dd %>%
    mutate(across(starts_with(c('GPP', 'q')), ~log(.))) %>%
    mutate(across(starts_with(c('epil','fila', 'light', 'K600', 'q')), ~scale(.)[,1]),
           across(starts_with(c('epil','fila', 'light', 'K600', 'q')), ~ . - min(., na.rm = T)),
    # mutate(across(starts_with(c('epil','fila')), ~scale(.)[,1]),
           across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
           year = lubridate::year(date),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),#!is.na(q.cms),
           !is.na(light))

sites <- unique(dd$site)
P_mod <- stan_model("code/model/stan_code/GPP_biomass_model_ar1_hierarchical.stan")
P_mod_nb <- stan_model("code/model/stan_code/GPP_NObiomass_model_ar1_hierarchical.stan")

# rstan::loo(dfit)
# dd <- select(dd, -K600) %>% rename(K600 = q.cms)
# run different model combinations: ####
mod_ests <- data.frame()
chains <- data.frame()
# Chlorophyll:
# fila
biomass <- matrix(dd$fila_chla_mgm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_chla')
ests <- extract_model_params(fit, 'fila', 'chla_mgm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains_K600(fit, 'fila', 'chla_mgm2')
chains <- bind_rows(chains, cc)

# epil
biomass <- matrix(dd$epil_chla_mgm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'epil_chla')
ests <- extract_model_params(fit, 'epil', 'chla_mgm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains_K600(fit, 'epil', 'chla_mgm2')
chains <- bind_rows(chains, cc)

# fila + epil
biomass <- matrix(c(dd$fila_chla_mgm2_fit, dd$epil_chla_mgm2_fit), ncol = 2)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_epil_chla')
ests <- extract_model_params(fit, 'fila_epil', 'chla_mgm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains_K600(fit, 'fila_epil', 'chla_mgm2')
chains <- bind_rows(chains, cc)

# biomass:
# fila
biomass <- matrix(dd$fila_gm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_gm2')
ests <- extract_model_params(fit, 'fila', 'gm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains_K600(fit, 'fila', 'gm2')
chains <- bind_rows(chains, cc)

# epil
biomass <- matrix(dd$epil_gm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'epil_gm2')
ests <- extract_model_params(fit, 'epil', 'gm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains_K600(fit, 'epil', 'gm2')
chains <- bind_rows(chains, cc)

# fila + epil
biomass <- matrix(c(dd$fila_gm2_fit, dd$epil_gm2_fit), ncol = 2)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_epil_gm2')
ests <- extract_model_params(fit, 'fila_epil', 'gm2', phi_corrected = TRUE)
ests$rmse = calculate_rmse(preds)
ests$r2_adj = calculate_r2_adj(preds, npar = nrow(ests))
mod_ests <- bind_rows(mod_ests, ests)
cc <- extract_model_chains_K600(fit, 'fila_epil', 'gm2')
chains <- bind_rows(chains, cc)


# model with no biomass
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1
dd$GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92
datlist <- list(
    N = nrow(dd),
    S = length(sites), ss = as.numeric(dd$site),
    K600 = dd$K600,
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
cc <- extract_model_chains_K600(fit, NA, NA)
chains <- bind_rows(chains, cc)

# compare parameter estimates across models

mod_ests <- mod_ests %>%
    mutate(parameter = sub('\\.{3}[0-9]+$', '', mod_ests$parameter),
           parameter = str_replace_all(parameter, '[\\[\\]]', ''))
mod_ests <- mod_ests %>%
    relocate(biomass_vars, units, parameter) %>%
    mutate(parameter = case_when(parameter == 'gamma1' ~ 'intercept',
                                 parameter == 'gamma2' ~ 'gamma_K600',
                                 parameter == 'gamma3' ~ 'gamma_light',
                                 parameter == 'gamma4' &
                                     biomass_vars %in% c('fila', 'fila_epil') ~ 'gamma_fila',
                                 parameter == 'gamma4' &
                                     biomass_vars == 'epil' ~ 'gamma_epil',
                                 parameter == 'gamma5' ~ 'gamma_epil',
                                 TRUE ~ parameter))
beepr::beep()

write_csv(mod_ests, 'data/model_fits/hierarchical_model_parameters_loggpp_logbm.csv')
write_csv(chains, 'data/model_fits/hierarchical_model_posterior_distributions_for_plot_loggpp_logbm.csv')
# write_csv(chains, 'data/model_fits/hierarchical_model_posterior_distributions_for_plot_fixedK.csv')
# write_csv(mod_ests, 'data/model_fits/hierarchical_model_parameters_logGPP.csv')

mod_ests <- read_csv( 'data/model_fits/hierarchical_model_parameters_logbm.csv')
chains <- read_csv('data/model_fits/hierarchical_model_posterior_distributions_for_plot_logbm.csv')

mod_ests %>%
    filter(!grepl('^beta', parameter) ,
           parameter != 'intercept') %>%
ggplot( aes(parameter, mean, fill = units, col = units))+
    geom_boxplot(aes(ymin = X2.5., lower = X25., middle = X50.,
                     upper = X75., ymax = X97.5.)) +
    ylim(0,2) +
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
