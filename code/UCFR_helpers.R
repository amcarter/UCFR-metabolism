# UCFR project helper functions

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


