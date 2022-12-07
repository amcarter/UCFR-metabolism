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

    dd$GPP_sd <- (dd$GPP.upper - dd$GPP.lower)/3.92

    datlist <- list(
        N = nrow(dd), K = ncol(biomass),
        S = length(sites), ss = as.numeric(dd$site),
        light = dd$light,
        biomass = biomass,
        P = dd$GPP, P_sd = dd$GPP_sd
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
    ss <- summary(fit, pars = c('gamma', 'beta', 'tau', 'sigma'))$summary %>%
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
    mutate(across(starts_with(c('epil','fila', 'light', 'K600')), ~scale(.)[,1]),
           across(starts_with(c('epil','fila', 'light', 'K600')), ~ . - min(., na.rm = T)),
    # mutate(across(starts_with(c('epil','fila')), ~scale(.)[,1]),
           across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
           year = lubridate::year(date),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light))

sites <- unique(dd$site)
P_mod <- stan_model("code/model/stan_code/GPP_biomass_model_hierarchical.stan")

# run on real data: ####

# fila + epil
biomass <- matrix(c(dd$fila_chla_mgm2_fit, dd$epil_chla_mgm2_fit), ncol = 2)
fit <- fit_biomass_model(dd, biomass, P_mod)
preds <- get_model_preds(fit, dd)
plot_model_fit(preds, dd, mod = 'fila_epil_chla_noar1')
ests <- extract_model_params(fit, 'fila_epil', 'chla_mgm2')
