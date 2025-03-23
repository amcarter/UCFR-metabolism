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
              biomass_chla = epil_chla_mgm2_fit + fila_chla_mgm2_fit,
              biomass_afdm = epil_gm2_fit + fila_gm2_fit,
              log_biomass = log(epil_chla_mgm2_fit + 0.1) +
                  log(fila_chla_mgm2_fit + 0.1))
q1 <- dat %>% select(starts_with(c('epil', 'fila', 'biomass')), -ends_with('_se')) %>%
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
           biomass_chla = biomass_chla + q1$biomass_chla,
           biomass_afdm = biomass_afdm + q1$biomass_afdm,
           across(starts_with(c('epil', 'fila', 'biomass_', 'light', 'GPP')), ~ log(.),
                  .names = 'log_{.col}'),
           GPP_sd = (GPP.upper - GPP.lower)/3.92,
           log_GPP_sd = (log_GPP.upper - log_GPP.lower)/3.92)
# light = light/max(light),
# log_light = log_light - min(log_light),
# log_light = log_light/max(log_light)) %>%
dd <- dat %>%
    mutate(across(contains(c('epil', 'fila', 'light', 'biomass_')), ~ scale(.)[,1])) %>%
    group_by(site, year) %>%
    mutate(
        across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
        site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_chla_mgm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    ungroup()

sites <- unique(dd$site)

# compile time series model
# ar1_lmod_ss <- stan_model("code/model/stan_code/AR1_linear_model_ss.stan",
#                           model_name = 'ar1_lmod_ss')
# ar1_lmod_ss <- stan_model("code/model/stan_code/AR1_error_model.stan",
#                           model_name = 'ar1_lmod_ss')

# simulate data: ####

N = nrow(dd)
K = 2
S = 6
ss = as.numeric(as.factor(dd$site))
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1
y_sd = abs(rnorm(N,0,0.2))

beta = rnorm(3)
tau = 0.5
beta0 = rnorm(S, beta[1], tau)
X <- as.matrix(select(dd,log_light, log_fila_chla_mgm2_fit), ncol = K)
X = matrix(rnorm(N*K), ncol = K, nrow = N)
theta = 0.6
sigma = 1

mu = c(X %*% beta[2:(K+1)])
y = vector(length = N)
for(i in 1:length(which(new_ts == 1))){
    start = which(new_ts == 1)[i]
    if(i == length(which(new_ts == 1))){
        end = N
    } else{
        end = which(new_ts == 1)[i + 1] - 1
    }
    y[start:end] = arima.sim(model = list(ma = theta), n = end-start+1) +
        mu[start:end] + beta0[ss[start]]
}

y = rnorm(N, y, y_sd)

ma1_mod <- stan_model('code/model/stan_code/MA1.stan',
                      model_name = 'ma_mod')
arma_mod <- stan_model('code/model/stan_code/ARMA1.stan',
                       model_name = 'arma_mod')
datlist <- list(
    N = N, K = K, S = S,
    ss = ss, X = X, new_ts = new_ts,
    y = y#, y_sd = y_sd
)

ma_fit <- sampling(
    ma1_mod,
    datlist,
    iter = 4000,
    chains = 4
)

print(ma_fit, pars = c('sigma', 'theta', 'gamma', 'beta', 'tau'))
pairs(ma_fit, pars = c('sigma', 'theta', 'gamma', 'tau'))

shinystan::launch_shinystan(ma_fit)

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
                   fila_chla = log_fila_chla_mgm2_fit,
                   biomass_chla = log_biomass_chla,
                   biomass_afdm = log_biomass_afdm) %>%
    mutate(across(-log_light, .fn = ~log_light * .,
                  .names = 'log_light_{.col}'))

colnames(master_X)

model_combinations <- list(1, c(1,2), c(1,3), c(1,4), c(1,5),
                           c(1,2,3), c(1,4,5),
                           c(1,2,8), c(1,3,9),
                           c(1,4,10), c(1,5,11),
                           c(1,2,3,8,9), c(1,4,5,10,11),
                           8,9,10,11, c(8,9), c(10,11),
                           c(1,6), c(1,7), c(1,6,12), c(1,7,13), 12,13)

# fit AR1 linear models ####
mod_ests <- data.frame()
chains <- data.frame()
mod_comp <- data.frame()

for(i in 1:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    X <- master_X[,model_combinations[[i]]]
    fit <- fit_biomass_model(ma1_mod, dd,
                             X = as.matrix(X),
                             log_GPP = TRUE,
                             iter = 4000)

    preds <- get_model_preds(fit, dd)
    p <- plot_model_fit(preds) +
        ggtitle(paste0('ma1_lmod: ', paste(colnames(X), collapse = ' + ')))
    print(p)
    ggsave(paste('figures/biomass_models/unconditioned_preds/ma1_lmod',
                 paste(colnames(X),collapse = '_'),'2.png'), p,
           width = 5, height = 7, units = 'in')

    ests <- extract_model_params(fit, pars = c('theta','gamma', 'beta', 'sigma', 'tau'),
                                 X, dd, phi_corrected = FALSE,
                                 extract.chains = TRUE)

    # preds <- mutate(preds, P_delta = P_mod - ests$ests$mean[1] *
    #                     c(preds$P_mod[1], preds$P_mod[1:(nrow(preds)-1)]))

    comp <- calculate_model_metrics(fit, X, dd, 'ma1_lmod')

    mod_ests <- bind_rows(mod_ests, ests$ests)
    chains <- bind_rows(chains, ests$chains)
    mod_comp <- bind_rows(mod_comp, comp)

}
beepr::beep(5)

write_csv(mod_ests, 'data/model_fits/ma1_log_model_parameter_ests_ss.csv')
write_csv(chains, 'data/model_fits/ma1_log_model_chains_ss.csv')
write_csv(mod_comp, 'data/model_fits/ma1_log_model_metrics_ss_cond.csv')

mod_ests <- read_csv('data/model_fits/ma1_log_model_parameter_ests_ss.csv')
mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)

mod_ests %>% filter(biomass_vars == 'log_light + epil_chla + fila_chla') %>%
    select(parameter, mean, X_2.5)
mod_comp %>% arrange(waic)
plot(mod_comp$waic, mod_comp$rmse)


# Fit models with holdout of 2021 years
mod_holdout_preds <- data.frame()
mod_post_preds <- data.frame()
mod_holdout_RMSE <- data.frame()

for(i in 1:length(model_combinations)){
    print(paste('model', i, 'of', length(model_combinations), sep = ' '))
    full_holdout_preds <- data.frame()
    full_post_preds <- data.frame()
    RMSE_holdout <- data.frame()
    for(s in sites){
        print(paste('Holdout site', s, sep = ' '))
        # remove data from 2021 at each site iteratively and fit model:
        holdout_rows <- which(dd$site == s & dd$year == 2021)

        X <- master_X[,model_combinations[[i]]]
        fit <- fit_biomass_model(ma1_mod, dd[-holdout_rows,],
                                 X = as.matrix(X[-holdout_rows,]),
                                 log_GPP = TRUE,
                                 iter = 4000)
        holdout_preds <- predict_holdout_ma(fit, dd, X, holdout_rows)
        # div <- stan_psum(fit)
        # RMSE_holdout <-
        #     data.frame(RMSE = sqrt(mean((holdout_preds$log_GPP - holdout_preds$estim)^2)),
        #                divergent_trans = div$divergent) %>%
        #     bind_rows(RMSE_holdout)

        full_holdout_preds <- bind_rows(full_holdout_preds, holdout_preds$preds)
        full_post_preds <- bind_rows(full_post_preds, holdout_preds$post_preds)
    }


    mod_holdout_preds <- full_holdout_preds %>%
        mutate(biomass_vars = paste(colnames(X), collapse = ' + '),
               model = fit@model_name) %>%
        bind_rows(mod_holdout_preds)

    mod_post_preds <- full_post_preds %>%
        mutate(biomass_vars = paste(colnames(X), collapse = ' + ')) %>%
        bind_rows(mod_post_preds)
    # mod_holdout_RMSE <-
    #     RMSE_holdout %>%
    #     mutate(site = sites,
    #            biomass_vars = rep(paste(colnames(X), collapse = ' + '), 6),
    #            model = rep(fit@model_name, 6)) %>%
    #     bind_rows(mod_holdout_RMSE)

}

beepr::beep(5)

write_csv(mod_holdout_preds, 'data/model_fits/ma1_log_model_holdout_predictions.csv')
write_csv(mod_post_preds, 'data/model_fits/ma1_log_model_holdout_post_predictions.csv')
write_csv(mod_holdout_RMSE, 'data/model_fits/ma1_log_model_holdout_RMSEs.csv')

mod_holdout_preds2 <- read_csv('data/model_fits/ar1_log_model_holdout_predictions.csv')

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
    group_by(biomass_vars) %>%
    mutate(estim = exp(estim),
           sqe = (GPP - estim)^2) %>%
    filter(grepl('biomass', biomass_vars)) %>%
    summarize(rmse = sqrt(mean(sqe))) %>% arrange(rmse)

N = 300
mod_post <- mod_post_preds %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))%>%
    mutate(model = case_when(biomass_vars == 'log_light' ~ '0. Baseline (Light)',
                             biomass_vars == 'log_light + epil_afdm + fila_afdm' ~
                                 '1. Light + Biomass',
                             biomass_vars == 'log_light + epil_afdm + fila_afdm + log_light_epil_afdm + log_light_fila_afdm' ~
                                 '3. Light + Biomass + \nLight \U00D7 Biomass',
                             biomass_vars == 'log_light_epil_chla' ~
                                 '2. Light \U00D7 Biomass',
                             TRUE ~ NA_character_),
           across(starts_with('X'), function(x) exp(x))) %>%
    filter(!is.na(model)) %>%
    select(site, date, model, paste0('X', sample(1:1000, N)))

colnames(mod_post) <- c('site', 'date', 'model', paste0('X', 1:N))
p <- mod_holdout_preds %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))%>%
    mutate(model = case_when(biomass_vars == 'log_light' ~ '0. Baseline (Light)',
                             biomass_vars == 'log_light + epil_afdm + fila_afdm' ~
                                 '1. Light + Biomass',
                             biomass_vars == 'log_light + epil_afdm + fila_afdm + log_light_epil_afdm + log_light_fila_afdm' ~
                                 '3. Light + Biomass + \nLight \U00D7 Biomass',
                             biomass_vars == 'log_light_epil_chla' ~
                                 '2. Light \U00D7 Biomass',
                             TRUE ~ NA_character_),
           across(starts_with('X'), function(x) exp(x))) %>%
    filter(!is.na(model)) %>%
    ggplot(aes(date, GPP)) +
    # geom_line() +
    # geom_ribbon(aes(ymin = exp(low), ymax = exp(high)),
    #             fill = alpha('#9E3A14', 0.03))+
    facet_grid(site~model) +
    theme_classic() +
    ylim(min(exp(mod_holdout_preds$low)),
         max(exp(mod_holdout_preds$high)))+
    ylab(expression("GPP"~(g~O[2]~m^-2~d^-1)))+
    # ylab("log(GPP)")+
    xlab('Date')+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'pt'))
err_col <- rgb(158/255, 57/255, 21/255)
for(i in 1:150){
    p <- p + geom_line(data = mod_post, aes_string(y = paste0('X', i)),
              col = alpha(err_col, 0.025))
}

p <- p + geom_line(aes(date, GPP), col = 'black', linewidth = 0.4)+
    geom_line(aes(y = exp(estim)), lty = 2, linewidth = 0.4)
png('figures/SI/holdout_site_predictions_MA.png', width = 7.5, height = 6.5,
    units = 'in', res = 300)
    p
dev.off()

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
