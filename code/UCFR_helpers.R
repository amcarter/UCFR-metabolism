# UCFR project helper functions

library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
library(tidyverse)



calculate_ts_mean_se <- function(ts){
    ts <- zoo::na.approx(ts)
    rho <- acf(ts, plot = FALSE)[[1]][2]
    n <- length(ts)
    neff <- n*(1-rho)/(1+rho)
    sem <- sd(ts/sqrt(neff))
    return(sem)

}

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
        # P_sd = dat$log_GPP_sd
    )


    if(model@model_name %in% c('ar1_lmod','ar_mod', 'ar1_lmod_ss','arma_mod', 'ma_mod')){
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
    preds <- rstan::summary(fit, pars = 'y_hat')$summary %>%
        data.frame() %>%
        select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
    preds <- bind_cols(dat, preds)

    return(preds)
}

predict_holdout <- function(fit, dat_full, X_full, holdout_rows){

    X <- as.matrix(X_full)
    site <- dat_full[holdout_rows,]$site[1]
    phi_post <- rstan::extract(fit, pars = 'phi')$phi
    draws <- nrow(phi_post)
    beta_post <- matrix(rstan::extract(fit, pars = 'gamma')$gamma[,-1],
                        nrow = draws)
    beta0_post <- rstan::extract(fit, pars = 'beta')$beta[,as.numeric(site)]
    sigma_post <- rstan::extract(fit, pars = 'sigma')$sigma


    # forecast the held out observations
    post_preds <- matrix(nrow = draws, ncol = length(holdout_rows))

    # fill in first observation:
    post_preds[,1] <- matrix(rep(log(dat_full$GPP[holdout_rows[1]]), draws),
                             nrow = draws, ncol = 1)
    for(i in 1:draws){
        for(t in 2:length(holdout_rows)){
            y_past <- as.double(post_preds[i, t-1])
            post_preds[i, t] <- X[holdout_rows[t],] %*% beta_post[i,] + beta0_post[i] +
                phi_post[i] * (y_past - X[holdout_rows[t-1],] %*% beta_post[i,] - beta0_post[i]) +
                rnorm(1, sd = sigma_post[i])
        }
    }

    preds <- dat_full[holdout_rows,] %>%
        select(site, date, GPP) %>%
        mutate(log_GPP = log(GPP),
               estim = apply(post_preds, 2, mean),
               low = apply(post_preds, 2, quantile, probs = 0.025),
               high = apply(post_preds, 2, quantile, probs = 0.975)
               )
    post_preds <- post_preds[sample(1:8000, 1000),]
    plot(preds$date, preds$log_GPP, type = 'l',
         ylim = c(min(preds$low), max(preds$high)))
    lines(preds$date, preds$estim, lty = 2)
    # for(i in sample(1:1000, 800)){
    #     lines(preds$date, exp(post_preds[i,]), col = alpha('brown', 0.02))
    # }
    # lines(preds$date, exp(preds$estim), lty = 2)
    # lines(preds$date, preds$GPP, lty = 2)
    polygon(c(preds$date, rev(preds$date)),
            c(preds$low, rev(preds$high)), col = alpha('brown', 0.3),
            border = NA)

    RMSE_forecast <- sqrt(mean((preds$log_GPP - preds$estim)^2))
    post_preds <- bind_cols(select(preds, site, date), data.frame(t(post_preds)))

    return(list(preds = preds,
                post_preds = post_preds)
    )

}
predict_holdout_ma <- function(fit, dat_full, X_full, holdout_rows){

    X <- as.matrix(X_full)
    site <- dat_full[holdout_rows,]$site[1]
    theta_post <- rstan::extract(fit, pars = 'theta')$theta
    draws <- nrow(theta_post)
    beta_post <- matrix(rstan::extract(fit, pars = 'gamma')$gamma[,-1],
                        nrow = draws)
    beta0_post <- rstan::extract(fit, pars = 'beta')$beta[,as.numeric(site)]
    sigma_post <- rstan::extract(fit, pars = 'sigma')$sigma


    # forecast the held out observations
    post_preds <- matrix(nrow = draws, ncol = length(holdout_rows))

    # fill in first observation:
    post_preds[,1] <- matrix(rep(log(dat_full$GPP[holdout_rows[1]]), draws),
                             nrow = draws, ncol = 1)
    for(i in 1:draws){
        eps_past <- as.double(post_preds[i, 1] -
                                  X[holdout_rows[1],] %*% beta_post[i,] -
                                  beta0_post[i])
        for(t in 2:length(holdout_rows)){
            eps <- rnorm(1, sd = sigma_post[i])
            post_preds[i, t] <- X[holdout_rows[t],] %*% beta_post[i,] +
                beta0_post[i] + theta_post[i] * eps_past + eps
            eps_past <- eps
        }
    }

    preds <- dat_full[holdout_rows,] %>%
        select(site, date, GPP) %>%
        mutate(log_GPP = log(GPP),
               estim = apply(post_preds, 2, mean),
               low = apply(post_preds, 2, quantile, probs = 0.025),
               high = apply(post_preds, 2, quantile, probs = 0.975)
               )
    post_preds <- post_preds[sample(1:8000, 1000),]
    plot(preds$date, preds$GPP, type = 'l',
         ylim = exp(c(min(preds$low), max(preds$high))))
    for(i in sample(1:1000, 800)){
        lines(preds$date, exp(post_preds[i,]), col = alpha('brown', 0.02))
    }
    lines(preds$date, exp(preds$estim), lty = 2)
    lines(preds$date, preds$GPP, lty = 2)
    # polygon(c(preds$date, rev(preds$date)),
    #         c(exp(preds$low), rev(exp(preds$high))), col = alpha('brown', 0.3),
    #         border = NA)

    RMSE_forecast <- sqrt(mean((preds$GPP - exp(preds$estim))^2))
    post_preds <- bind_cols(select(preds, site, date), data.frame(t(post_preds)))

    return(list(preds = preds,
                post_preds = post_preds))

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

    return(r2)
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
        r2 = calculate_r2_adj(fit, dd, log = log),
        rmse = calculate_rmse(fit, dd, log = log)
    )

    return(comp)

}


stan_psum <- function(fit){

    iter <- fit@stan_args[[1]]$iter
    s_init <- as.data.frame(rstan::summary(fit)$summary)
    s_init$pars <- row.names(s_init)
    s_init$n_eff_pct <- s_init$n_eff/iter  ## effective samples are the number of independent samples with the same estimation power as the N autocorrelated samples
    s_init$n_eff_less10pct <- ifelse(s_init$n_eff_pct < 0.10, yes = "true", no = "false") # 10% is often used as a threshold, below which the chains for a parameter did not properly converge
    s <- s_init[,c("pars","Rhat","n_eff_less10pct")]
    s <- s[!grepl('^mu', row.names(s_init)),]
    s <- s[s$Rhat > 1.1 | s$n_eff_less10pct == 'true',]
    spars <- get_sampler_params(fit, inc_warmup = FALSE)
    divtrans <- sum(sapply(spars, function(x) sum(x[,'divergent__'])))

    return(list(par_conv = s,
                divergent = divtrans)
    )
}
