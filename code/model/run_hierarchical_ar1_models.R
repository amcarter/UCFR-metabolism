library(rstan)
options(mc.cores = parallel::detectCores())
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
fit_biomass_model <- function(dd, biomass, model){

    new_ts <- rep(0, nrow(dd))
    new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1
    dd$GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92

    datlist <- list(
        N = nrow(dd), K = ncol(biomass),
        S = length(sites), ss = as.numeric(dd$site),
        K600 = dd$K600,
        light = dd$light,
        biomass = biomass,
        P = dd$GPP, P_sd = dd$GPP_sd,
        new_ts = new_ts
    )

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

extract_model_params <- function(fit, bmv = 'fila + epil', units = 'gm2'){
    ss <- summary(fit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma'))$summary %>%
        data.frame()
    ss$biomass_vars = bmv
    ss$units = units
    ss$parameter = row.names(ss)
    row.names(ss) <- NULL

    return(ss)
}

# prep data ####
dd <- read_csv('data/biomass_metab_model_data.csv')
dd <- dd %>%
    # mutate(across(starts_with(c('epil','fila', 'light', 'K600')), ~scale(.)[,1]),
    mutate(across(starts_with(c('epil','fila')), ~scale(.)[,1]),
           across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
           year = lubridate::year(date),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light))

sites <- unique(dd$site)
P_mod <- stan_model("code/model/GPP_biomass_model_ar1_hierarchical.stan")
P_mod_nb <- stan_model("code/model/GPP_NObiomass_model_ar1_hierarchical.stan")

# run on real data: ####

biomass = matrix(dd$fila_chla_mgm2_fit, ncol = 1)
dfit <- fit_biomass_model(dd, biomass, P_mod)

print(dfit, pars = c('phi', 'gamma', 'beta', 'tau', 'sigma'))
plot(dfit, pars = c('phi', 'gamma[2]', 'gamma[3]', 'gamma[4]', 'sigma', 'tau'),
     show_density = TRUE)
plot(dfit, pars = c('gamma[1]', 'beta'), show_density = TRUE)

# shinystan::launch_shinystan(dfit)
preds <- get_model_preds(dfit, dd)

# run different model combinations: ####
mod_ests <- data.frame()
# Chlorophyll:
# fila
biomass <- matrix(dd$fila_chla_mgm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_chla')
ests <- extract_model_params(fit, 'fila', 'chla_mgm2')
mod_ests <- bind_rows(mod_ests, ests)
# epil
biomass <- matrix(dd$epil_chla_mgm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'epil_chla')
ests <- extract_model_params(fit, 'epil', 'chla_mgm2')
mod_ests <- bind_rows(mod_ests, ests)
# fila + epil
biomass <- matrix(c(dd$fila_chla_mgm2_fit, dd$epil_chla_mgm2_fit), ncol = 2)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_epil_chla')
ests <- extract_model_params(fit, 'fila_epil', 'chla_mgm2')
mod_ests <- bind_rows(mod_ests, ests)

# biomass:
# fila
biomass <- matrix(dd$fila_gm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_gm2')
ests <- extract_model_params(fit, 'fila', 'gm2')
mod_ests <- bind_rows(mod_ests, ests)
# epil
biomass <- matrix(dd$epil_gm2_fit, ncol = 1)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'epil_gm2')
ests <- extract_model_params(fit, 'epil', 'gm2')
mod_ests <- bind_rows(mod_ests, ests)
# fila + epil
biomass <- matrix(c(dd$fila_gm2_fit, dd$epil_gm2_fit), ncol = 2)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_epil_gm2')
ests <- extract_model_params(fit, 'fila_epil', 'gm2')
mod_ests <- bind_rows(mod_ests, ests)


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
ests <- extract_model_params(fit, NA, NA)
mod_ests <- bind_rows(mod_ests, ests)

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
write_csv(mod_ests, 'data/model_fits/hierarchical_model_parameters.csv')

mod_ests %>%
    filter(!grepl('^beta', parameter) ,
           parameter != 'intercept') %>%
ggplot( aes(parameter, mean, fill = units))+
    geom_boxplot() +
    ylim(0,2)


mod_ests %>% filter(parameter %in% c('gamma_epil', 'gamma_fila'))

par(mfrow = c(1,3),
    mar = c(4,0.3,5,0.3),
    oma = c(1,1,1,1))

plot(mod_ests$mean[mod_ests$parameter == 'gamma_light'], seq(1:7),
     xlim = c(-6.5, 6.5), pch = 19,
     xlab = 'Light', yaxt = 'n', bty = 'n')
segments(x0 = mod_ests$X2.5.[mod_ests$parameter == 'gamma_light'], y0 = seq(1:7),
         x1 = mod_ests$X97.5.[mod_ests$parameter == 'gamma_light'], y1 = seq(1:7))
abline(v = 0)
mtext('Coefficient', line = 3.1, adj = 0)
mtext('Estimates', line = 1.9, adj = 0)
plot(mod_ests$mean[mod_ests$parameter == 'gamma_K600'], seq(1:7),
     pch = 19, xlab = 'K600', xlim = c(-0.5, 0.5),
     yaxt = 'n', bty = 'n')
segments(x0 = mod_ests$X2.5.[mod_ests$parameter == 'gamma_K600'], y0 = seq(1:7),
         x1 = mod_ests$X97.5.[mod_ests$parameter == 'gamma_K600'], y1 = seq(1:7))
abline(v = 0)
plot(mod_ests$mean[mod_ests$parameter == 'gamma_fila'], seq(1:4),
     pch = as.numeric(factor(mod_ests$units[mod_ests$parameter == 'gamma_epil'])),
     xlab = 'Connectivity',
     yaxt = 'n', bty = 'n')
segments(x0 = mod_table$C_mean - mod_table$C_se, y0 = seq(1:nrow(mod_table)),
         x1 = mod_table$C_mean + mod_table$C_se, y1 = seq(1:nrow(mod_table)),
         col = mod_table$C_col)
abline(v = 0)
legend(x = -.2, y = nrow(mod_table) + 2.5,
       legend = c('width to area', 'precipitation',
                  'terrestrial NPP', 'con. drainage dens'),
       pch = c(18, 17, 16, 15), col = c('black', rep('brown3', 3)),
       bty = 'n', xpd = TRUE)
