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
fit_biomass_model <- function(model, dat, X, log_GPP = TRUE,
                              iter = 2000, ss = FALSE){

    if(!log_GPP){
        dat <- mutate(dat,
                      log_GPP = GPP,
                      log_GPP_sd = GPP_sd)
    }

    datlist <- list(
        N = nrow(X), K = ncol(X),
        S = length(unique(dat$site)),
        ss = as.numeric(dat$site),
        X = X, P = dat$log_GPP
    )


    if(model@model_name %in% c('ar1_lmod','PI_curve', 'ar1_lmod_ss')){
        new_ts <- rep(0, nrow(dat))
        new_ts[rle2(paste0(dat$year, dat$site))$starts] <- 1
        datlist <- c(datlist, list(new_ts = new_ts))
        if(model@model_name == 'PI_curve'){
            datlist <- c(datlist, list(light = dat$light - min(dat$light)))
        }
    }

    if(ss){ datlist <- c(datlist, list(P_sd = dat$log_GPP_sd))}

    fit <- sampling(
        model,
        datlist,
        chains = 4,
        iter = iter
    )

    return(fit)
}
# fit_biomass_model <- function(dd, biomass, model, K600 = TRUE){
#
#     new_ts <- rep(0, nrow(dd))
#     new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1
#     dd$GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92
#
#     datlist <- list(
#         N = nrow(dd), K = ncol(biomass),
#         S = length(sites), ss = as.numeric(dd$site),
#         K600 = dd$K600,
#         light = dd$light,
#         biomass = biomass,
#         P = dd$GPP,
#         P_sd = dd$GPP_sd,
#         new_ts = new_ts
#     )
#     if(!K600){
#         datlist <- list(
#             N = nrow(dd), K = ncol(biomass),
#             S = length(sites), ss = as.numeric(dd$site),
#             light = dd$light,
#             biomass = biomass,
#             P = dd$GPP,
#             P_sd = dd$GPP_sd,
#             new_ts = new_ts
#         )
#
#     }
#
#     fit <- sampling(
#         model,
#         data = datlist
#     )
#
#     return(fit)
# }

get_model_preds <- function(fit, dat){
    preds <- rstan::summary(fit, pars = 'y_tilde')$summary %>%
        data.frame() %>%
        select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
    preds <- bind_cols(dat, preds)

    return(preds)
}

plot_model_fit <- function(preds, log_ests = TRUE, mod = NA,
                           scales = 'free_x', write = FALSE){
    if(log_ests){
        preds <- mutate(preds,
                        across(starts_with('P_mod'), ~exp(.)))
    }

    p <- ggplot(preds, aes(date, GPP)) +
        geom_line(aes(y = P_mod), col = 'steelblue', linewidth = 0.5) +
        geom_ribbon(aes(ymin = P_mod_lower, ymax = P_mod_upper),
                     fill = alpha('steelblue', 0.3), col = NA)+
        geom_point(size = 0.3) +
        geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), linewidth = 0.3) +
        facet_grid(site~year, scales = scales) +
        ylab(expression(paste('GPP (g ', O[2], ' ', m^-2, d^-1, ')')))+
        xlab('')+
        theme_bw()

    if(write){
    ggsave(paste0('figures/biomass_models/hierarchical_biomass_model_fit_',mod,'.png'),
           plot = p, width = 6, height = 7, units = 'in', dpi = 300)
    }else{
        return(p)
    }
}

extract_model_params <- function(fit, X, dd,
                                 pars = c('phi', 'gamma','beta', 'tau', 'sigma'),
                                 extract.chains = FALSE,
                                 phi_corrected = FALSE){
    avail_pars <- fit@model_pars

    if(fit@model_name == 'PI_curve'){
        pars = c('alpha', 'Pmax', 'phi', 'beta', 'beta_s', 'tau_b', 'sigma')
    }

    if(!("phi" %in% avail_pars)){
        pars = grep('phi', pars, value = T, invert = T)
    }

    ss <- summary(fit, pars = pars)$summary %>%
        data.frame() %>%
        select(mean, se_mean, sd, X_2.5 = X2.5., X_50 = X50., X_97.5 = X97.5.)
    ss$parameter = row.names(ss)
    row.names(ss) <- NULL

    if(phi_corrected){
        ss <- ss %>%
            mutate(across(-parameter,
                          ~case_when(parameter != 'phi' ~
                                     ./(1-ss$mean[which(ss$parameter == 'phi')[1]]),
                                     TRUE ~ .)))
    }

    # rename parameters
    ss$parameter <- sub('^beta\\[([0-9]+)\\]$', 'beta_\\1', ss$parameter)
    ss$parameter[grepl('gamma', ss$parameter)] <- paste0('gamma_',
                                                         c('intercept', colnames(X)))

    ss$biomass_vars = paste(colnames(X), collapse = ' + ')
    ss$model <- fit@model_name
    ss$r2_adj = calculate_r2_adj(fit, dd)
    ss$rmse = calculate_rmse(fit, dd)

    if(fit@model_name == 'PI_curve'){
        ss$parameter <- sub('^beta$', 'intercept', ss$parameter)
        ss$parameter <- sub('^beta_s\\[([0-9]+)\\]$', 'beta_\\1', ss$parameter)
        ss$parameter[grepl('alpha', ss$parameter)] <- paste0('alpha_', colnames(X))
        ss$parameter[grepl('Pmax', ss$parameter)] <- paste0('Pmax_', colnames(X))
    }

    if(extract.chains){
        cc <- rstan::extract(fit, pars = pars) %>% data.frame()
        colnames(cc) <- sub('^beta\\.([0-9]+)$', 'beta_\\1', colnames(cc))
        colnames(cc)[grepl('gamma', colnames(cc))] <- paste0('gamma_', c('intercept', colnames(X)))
        if(fit@model_name == 'PI_curve'){
            colnames(cc) <- sub('^beta$', 'intercept', colnames(cc))
            colnames(cc) <- sub('^beta_s\\.([0-9]+)$', 'beta_\\1', colnames(cc))
            colnames(cc)[grepl('alpha', colnames(cc))] <- paste0('alpha_', colnames(X))
            colnames(cc)[grepl('Pmax', colnames(cc))] <- paste0('Pmax_', colnames(X))
        }
        if(phi_corrected){
            cc <- mutate(cc,across(-phi,
                                   ~./(1-ss$mean[which(ss$parameter == 'phi')[1]])))
        }

        cc <- pivot_longer(cc, cols = everything(),
                           names_to = 'parameter',
                           values_to = 'chains') %>%
            mutate(biomass_vars = paste(colnames(X), collapse = ' + '),
                   model = fit@model_name)

        return(list(ests = ss,
                    chains = cc))
    }

    return(ss)
}


calculate_r2_adj <- function(fit, dd, log = TRUE){
    preds <- get_model_preds(fit, dd)
    npar <- length(grep('mu|y_tilde|lp|log_lik|rates',
                        row.names(summary(fit)$summary), invert = T))
    N = nrow(preds)
    if(log) preds$P_mod <- exp(preds$P_mod)
    r2 = 1 - sum((preds$GPP - preds$P_mod)^2)/
        sum((preds$GPP - mean(preds$GPP))^2)

    r2_adj = 1 - ((1-r2)*(N-1))/(N-npar-1)

    return(r2_adj)
}

calculate_rmse <- function(fit = NULL, dd = NULL, preds = NULL, log = TRUE){
    if(is.null(preds)) preds <- get_model_preds(fit, dd)
    if(log) preds$P_mod <- exp(preds$P_mod)

    N = nrow(preds)
    rmse = sqrt((1/N) * sum((preds$GPP - preds$P_mod)^2))

    return(rmse)
}

calculate_model_metrics <- function(fit, X, dd, model, log = TRUE){

    log_lik <- rstan::extract(fit, pars = 'log_lik')$log_lik
    waic <- loo::waic(log_lik)$estimates
    r_eff <- loo::relative_eff(exp(log_lik), chain_id = 1:nrow(log_lik))

    loo_cv <- loo::loo(log_lik, r_eff = r_eff)$estimates
    comp <- data.frame(
        model = model,
        covariates = paste(colnames(X), collapse = ' + '),
        waic = waic[3,1],
        waic_se = waic[3,2],
        loo = loo_cv[3,1],
        loo_se = loo_cv[3,2],
        r2adj = calculate_r2_adj(fit, dd, log = log),
        rmse = calculate_rmse(fit, dd, log = log)
    )

    return(comp)

}
