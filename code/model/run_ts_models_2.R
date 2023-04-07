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
P[new_ts == 1] <- beta[ss[new_ts == 1]] + X[new_ts == 1,]%*% gamma[2:(K+1)]
mu = rep(NA_real_, N)

for(n in 1:N){
    if(new_ts[n] == 1){
        # restart the AR process on each new time series
        mu[n] = rnorm(1, P[n], sigma)
    }
    else{
        mu[n] = rnorm(1, beta[ss[n]] + X[n,] %*% gamma[2:(K+1)] + phi * P[n-1], sigma)
    }

    P[n] = rnorm(1, mu[n], P_sd[n]);
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
    chains = 4
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
    # write_lines(paste('model', i, 'of', length(model_combinations), ' vars = ',
    #                   colnames(X), sep = ' '),
    #             'model_warnings.txt', append = TRUE)
    # fit <- tryCatch({
    fit <- fit_biomass_model(ar1_lmod_ss, dd, X = as.matrix(X), iter = 4000)
    # }, warning = function(w) write_lines(w, 'model_warnings.txt', append = TRUE))
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

write_csv(mod_ests, 'data/model_fits/ar1_log_model_parameter_ests_ss.csv')
write_csv(chains, 'data/model_fits/ar1_log_model_chains_ss.csv')
write_csv(mod_comp, 'data/model_fits/ar1_log_model_metrics_ss_cond.csv')

mod_ests %>% distinct(biomass_vars, .keep_all = TRUE) %>%
    arrange(r2_adj)

mod_ests %>% filter(biomass_vars == 'log_light + epil_chla + fila_chla') %>%
    select(parameter, mean)
mod_comp %>% arrange(waic)
plot(mod_comp$waic, mod_comp$rmse)


