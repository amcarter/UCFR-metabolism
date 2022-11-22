# compare a basic linear model of GPP to one with an AR1 term:

library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)

dd <- read_csv('data/biomass_metab_model_data.csv')
dd <- dd %>%
    mutate(across(starts_with('epil_chla'),
                  function(x) x/max(dd$epil_chla_mgm2_fit, na.rm = T)),
           across(starts_with('epil_gm'),
                  function(x) x/max(dd$epil_gm2_fit, na.rm = T)),
           across(starts_with('fila_chla'),
                  function(x) x/max(dd$fila_chla_mgm2_fit, na.rm = T)),
           across(starts_with('fila_gm'),
                  function(x) x/max(dd$fila_gm2_fit, na.rm = T)))


# insert fake data test here

## run on datasets ####
# compile stan models
l_mod <- stan_model("code/model/GPP_biomass_model_one_site.stan")
ar1_mod <- stan_model("code/model/GPP_biomass_model_ar1.stan")


compare_fits <- function(dd, s){
    BN <- filter(dd, site == s) %>%
        mutate(across(c(starts_with('GPP'), 'light', starts_with('epil'),
                               starts_with('fila')),
                      zoo::na.approx, na.rm = F)) %>%
        filter(!is.na(GPP) & !is.na(epil_gm2_fit) & !is.na(fila_gm2_fit) &
                   !is.na(light))

    GPP <- BN$GPP
    GPP_sd <- (BN$GPP.upper - BN$GPP.lower)/3.92
    mdat <- list(N = nrow(BN),
                   K = 3,
                   light = BN$light,
                   biomass = BN$fila_chla_mgm2_fit,
                   P = GPP,
                   P_sd = GPP_sd)
    l_fit <- sampling(l_mod, data = mdat)
    # plot(l_fit, pars = 'gamma')
    ests <- rstan::summary(l_fit, pars = 'gamma')$summary %>%
        data.frame() %>%
        select(mean, se_mean, sd, lower = 'X2.5.', upper = 'X97.5.') %>%
        mutate(model = 'lm',
               site = s)

    ar1_fit <- sampling(ar1_mod, data = mdat)
    ests <- rstan::summary(ar1_fit, pars = c('gamma', 'sigma', 'phi'))$summary %>%
        data.frame() %>%
        select(mean, se_mean, sd, lower = 'X2.5.', upper = 'X97.5.') %>%
        mutate(model = 'ar1',
               site = s) %>%
        bind_rows(ests)

    return(ests)
}

ests <- data.frame()
e <- compare_fits(dd, 'PL')
ests <- bind_rows(ests, e)
e <- compare_fits(dd, 'DL')
ests <- bind_rows(ests, e)
e <- compare_fits(dd, 'GR')
ests <- bind_rows(ests, e)
e <- compare_fits(dd, 'GC')
ests <- bind_rows(ests, e)
e <- compare_fits(dd, 'BM')
ests <- bind_rows(ests, e)
e <- compare_fits(dd, 'BN')
ests <- bind_rows(ests, e)

ests <- ests %>%
    mutate(parameter = sub('\\.{3}[0-9]+$', '', row.names(ests)),
           parameter = str_replace_all(parameter, '[\\[\\]]', ''))

row.names(ests) = NULL

ests %>% filter(parameter == 'phi')

png('figures/parameter_comparison_lm_ar1.png')
ests %>%
    select(-se_mean, -sd) %>%
    pivot_wider(id_cols = c(model, site), names_from = parameter,
                values_from = c(mean, lower, upper)) %>%
    mutate(across(ends_with(c('gamma1', 'gamma2', 'gamma3')),
                  ~case_when(model == 'ar1' ~ ./(1-mean_phi),
                             model == 'lm' ~ .))) %>%
    select(-ends_with(c('sigma', 'phi'))) %>%
    pivot_longer(cols = c(-model, -site), names_to = c('stat', 'parameter'),
                 values_to = 'estimate', names_sep = '_') %>%
    ggplot(aes(parameter, estimate, fill = model)) +
    geom_boxplot() +
    facet_wrap(.~site)
dev.off()

print(ar1_fit, pars = c('phi', 'sigma', 'gamma'))
plot(ar1_fit, pars = c('phi', 'sigma', 'gamma'))
pairs(ar1_fit, pars = c('phi', 'sigma', 'gamma'))
