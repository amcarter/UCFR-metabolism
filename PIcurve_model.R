library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)


# functions ####
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
get_model_preds <- function(fit, dat){
    preds <- rstan::summary(fit, pars = 'y_tilde')$summary %>%
        data.frame() %>%
        select(P_mod = mean, P_mod_upper = 'X97.5.', P_mod_lower = 'X2.5.')
    preds <- bind_cols(dat, preds)

    return(preds)
}

# data ####
dd <- read_csv('data/biomass_metab_model_data.csv')
dd <- dd %>%
    mutate(#across(starts_with(c('epil','fila', 'light', 'K600')), ~scale(.)[,1]),
           # across(starts_with(c('epil','fila', 'light', 'K600')), ~ . - min(., na.rm = T)),
           # mutate(across(starts_with(c('epil','fila')), ~scale(.)[,1]),
           across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
           year = lubridate::year(date),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    mutate(GPP_sd = (GPP.upper - GPP.lower)/3.92)
# dd <- filter(dd,site == "GC")
GPP = dd$GPP
sites <- unique(dd$site)

PImod <- stan_model('code/model/stan_code/PIcurve_model.stan')


datlist <- list(
    N = nrow(dd),
    S = length(sites),
    ss = as.numeric(dd$site),
    light = dd$light,
    K600 = dd$K600,
    biomass = dd$epil_chla_mgm2_fit + dd$fila_chla_mgm2_fit,
    P = dd$GPP#, P_sd = dd$GPP_sd
)

fit <- sampling(
    PImod,
    data = datlist
    # chains = 1
)

plot(fit, pars = c('Pmax', 'alpha','beta', 'sigma','tau_Pmax', 'tau_a'))
print(fit, pars = c('Pmax', 'tau_Pmax', 'Pmax_s', 'tau_a', 'alpha', 'alpha_s'))
shinystan::launch_shinystan(fit)

preds <- get_model_preds(fit, dd)
print(fit)
calculate_r2_adj(preds, 18)
calculate_rmse(preds)

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

plot_model_fit(preds, dd, 'PI_curve_fit')


lmmod <- stan_model('code/model/stan_code/GPP_biomass_model_hierarchical.stan')
datlist <- list(
    N = nrow(dd),
    S = length(sites),
    ss = as.numeric(dd$site),
    light = dd$light,
    K600 = dd$K600,
    biomass = dd$epil_chla_mgm2_fit + dd$fila_chla_mgm2_fit,
    P = dd$GPP, P_sd = dd$GPP_sd
)

fit <- sampling(
    lmmod,
    data = datlist
)

print(fit)
preds <- get_model_preds(fit, dd)
calculate_r2_adj(preds, 12)
calculate_rmse(preds)

lmimod <- stan_model('code/model/stan_code/linear_model_with_interactions.stan')
datlist <- list(
    N = nrow(dd),
    S = length(sites),
    ss = as.numeric(dd$site),
    light = dd$light,
    K600 = dd$K600,
    biomass = dd$epil_chla_mgm2_fit + dd$fila_chla_mgm2_fit,
    P = dd$GPP, P_sd = dd$GPP_sd
)

fit <- sampling(
    lmimod,
    data = datlist
)

print(fit)
preds <- get_model_preds(fit, dd)
calculate_r2_adj(preds, 13)
calculate_rmse(preds)


lmiomod <- stan_model('code/model/stan_code/linear_model_only_interactions.stan')
datlist <- list(
    N = nrow(dd),
    S = length(sites),
    ss = as.numeric(dd$site),
    light = dd$light,
    K600 = dd$K600,
    biomass = dd$epil_chla_mgm2_fit + dd$fila_chla_mgm2_fit,
    P = dd$GPP, P_sd = dd$GPP_sd
)

fit <- sampling(
    lmiomod,
    data = datlist
)

print(fit)
preds <- get_model_preds(fit, dd)
calculate_r2_adj(preds, 11)
calculate_rmse(preds)

