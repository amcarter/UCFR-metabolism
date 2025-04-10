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
ar1_lmod_ss <- stan_model("code/model/stan_code/AR1_error_model.stan",
                       model_name = 'ar1_lmod_ss')

# simulate data: ####

N = nrow(dd)
K = 2
S = 6
ss <- as.numeric(as.factor(dd$site))
X <- as.matrix(select(dd,log_light, log_fila_chla_mgm2_fit), ncol = K)
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1

# parameters
gamma = rnorm(K+1,0,2)
phi = rbeta(1,1,1)
tau = abs(rnorm(1,0, 2))
sigma = abs(rnorm(1,0,2))
P_sd = abs(rnorm(N,0,0.2))
beta = rnorm(S, gamma[1], tau);

P = rep(NA_real_, N)
P_state = rep(NA_real_, N)
P[new_ts == 1] <- beta[ss[new_ts == 1]] + X[new_ts == 1,]%*% gamma[2:(K+1)]
mu = beta[ss] + X %*% gamma[2:(K+1)]

for(n in 1:N){
    if(new_ts[n] == 1){
        # restart the AR process on each new time series
        P_state[n] = rnorm(1, P[n], sigma)
    }
    else{
        P_state[n] = rnorm(1, mu[n] + phi * (P[n-1] - mu[n-1]), sigma)
    }

    P[n] = rnorm(1, P_state[n], P_sd[n]);
}



datlist <- list(
    N = N, K = K, S = S,
    ss = ss, X = X, new_ts = new_ts,
    P = P, P_sd = P_sd
)

simfit <- sampling(
    ar1_lmod_ss,
    datlist,
    # iter = 4000,
    chains = 1
)

shinystan::launch_shinystan(simfit)
print(simfit, pars = c('gamma', 'sigma', 'tau', 'phi', 'beta'))
gamma
sigma
tau
phi
beta

ggplot(dd, aes(date, log_fila_chla_mgm2_fit))+
    geom_line() + facet_grid(site~year, scales = 'free')
dd %>% mutate(log_GPP = P) %>%
ggplot( aes(log_fila_chla_mgm2_fit, log_GPP, col = site))+
    geom_point()

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
    fit <- fit_biomass_model(ar1_lmod_ss, dd,
                             X = as.matrix(X),
                             log_GPP = TRUE,
                             iter = 4000)

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

write_csv(mod_ests, 'data/model_fits/ar1_log_model_parameter_ests_ss2.csv')
write_csv(chains, 'data/model_fits/ar1_log_model_chains_ss2.csv')
write_csv(mod_comp, 'data/model_fits/ar1_log_model_metrics_ss_cond2.csv')

# mod_ests <- read_csv('data/model_fits/ar1_log_model_parameter_ests_ss.csv')
mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)

mod_ests %>% filter(biomass_vars == 'log_light + epil_chla + fila_chla') %>%
    select(parameter, mean)
mod_comp %>% arrange(waic)
plot(mod_comp$waic, mod_comp$rmse)


# Fit models with holdout of 2021 years
mod_holdout_preds <- data.frame()
mod_holdout_RMSE <- data.frame()

for(i in 1:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    full_holdout_preds <- data.frame()
    RMSE_holdout <- data.frame()
    for(s in sites){
        print(paste('Holdout site', s, sep = ' '))
        # remove data from 2021 at each site iteratively and fit model:
        holdout_rows <- which(dd$site == s & dd$year == 2021)

        X <- master_X[,model_combinations[[i]]]
        fit <- fit_biomass_model(ar1_lmod_ss, dd[-holdout_rows,],
                             X = as.matrix(X[-holdout_rows,]),
                             log_GPP = TRUE,
                             iter = 4000)
        holdout_preds <- predict_holdout(fit, dd, X, holdout_rows)
        div <- stan_psum(fit)
        RMSE_holdout <-
            data.frame(RMSE = sqrt(mean((holdout_preds$log_GPP - holdout_preds$estim)^2)),
                       divergent_trans = div$divergent) %>%
            bind_rows(RMSE_holdout)

        full_holdout_preds <- bind_rows(full_holdout_preds, holdout_preds)
    }


    mod_holdout_preds <- full_holdout_preds %>%
        mutate(biomass_vars = paste(colnames(X), collapse = ' + '),
               model = fit@model_name) %>%
        bind_rows(mod_holdout_preds)
    mod_holdout_RMSE <-
        RMSE_holdout %>%
        mutate(site = sites,
               biomass_vars = rep(paste(colnames(X), collapse = ' + '), 6),
               model = rep(fit@model_name, 6)) %>%
        bind_rows(mod_holdout_RMSE)

}

beepr::beep(5)

write_csv(mod_holdout_preds, 'data/model_fits/ar1_log_model_holdout_predictions2.csv')
write_csv(mod_holdout_RMSE, 'data/model_fits/ar1_log_model_holdout_RMSEs2.csv')

mod_holdout_preds <- read_csv('data/model_fits/ar1_log_model_holdout_predictions.csv')

mod_holdout_preds  %>% tibble() %>%
    group_by(biomass_vars, model, site) %>%
    mutate(se = (log_GPP - (estim))^2) %>%
    summarize(mse = mean(se)) %>%
    mutate(rmse = sqrt(mse)) %>%
    group_by(biomass_vars) %>%
    summarize(RMSE = mean(rmse)) %>% arrange(RMSE)
    ggplot(aes(biomass_vars, rmse, col = site)) +
    geom_point()

    mod_holdout_preds %>%
        mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))%>%
        mutate(model = case_when(biomass_vars == 'log_light' ~ '0. Baseline (Light)',
                                 biomass_vars == 'log_light + fila_chla' ~
                                     '1. Light + Biomass',
                                 biomass_vars == 'log_light + fila_chla + log_light_fila_chla' ~
                                     '3. Light + Biomass + \nLight \U00D7 Biomass',
                                 biomass_vars == 'log_light_fila_chla' ~
                                     '2. Light \U00D7 Biomass',
                                TRUE ~ NA_character_)) %>%
        filter(!is.na(model)) %>%
        ggplot(aes(date, GPP)) +
        geom_line() +
        geom_line(aes(y = exp(estim)), lty = 2)+
        geom_ribbon(aes(ymin = exp(low), ymax = exp(high)),
                    fill = alpha('brown', 0.3))+
        facet_grid(site~model) +
        theme_classic() +
        # ylim(0.88, 3.5)+
        # ylab(expression("log(GPP)"~(g~O[2]~m^-2~d^-1)))+
        ylab("log(GPP)")+
        xlab('Date')+
        theme(panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, units = 'pt'))

png('figures/SI/holdout_site_predictions.png', width = 7.5, height = 6.5,
    units = 'in', res = 300)

    # p <-
    mod_holdout_preds %>%
        mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))%>%
        mutate(model = case_when(biomass_vars == 'log_light' ~ '0. Baseline (Light)',
                                 biomass_vars == 'log_light + fila_chla' ~
                                     '1. Light + Biomass',
                                 biomass_vars == 'log_light + fila_chla + log_light_fila_chla' ~
                                     '3. Light + Biomass + \nLight \U00D7 Biomass',
                                 biomass_vars == 'log_light_fila_chla' ~
                                     '2. Light \U00D7 Biomass',
                                TRUE ~ NA_character_)) %>%
        filter(!is.na(model)) %>%
        ggplot(aes(date, log_GPP)) +
        geom_line() +
        geom_line(aes(y = estim), lty = 2)+
        geom_ribbon(aes(ymin = low, ymax = high),
                    fill = alpha('brown', 0.3))+
        facet_grid(site~model) +
        theme_classic() +
        ylim(0.88, 3.5)+
        # ylab(expression("log(GPP)"~(g~O[2]~m^-2~d^-1)))+
        ylab("log(GPP)")+
        xlab('Date')+
        theme(panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, units = 'pt'))



    # dummy_df <- data.frame(GPP = rep(c('Actual', 'Predicted'), each = 2),
    #                        x = c(1,2,1,2), y = c(1,1,1,2))
    # dummy_plot <- ggplot(dummy_df, aes(x,y, lty = GPP)) +
    #     geom_line() + theme_classic()+
    #     # geom_ribbon(aes(ymin = y/2, ymax = 2*y, fill = GPP))+
    #     # scale_fill_manual(values = c('transparent', alpha('brown', 0.3)),
    #     #                   name = "") +
    #     theme(legend.direction = 'horizontal')
    #
    # lgnd <- ggpubr::get_legend(dummy_plot)
    # pp <- ggpubr::ggarrange(p, legend.grob = lgnd)

dev.off()
