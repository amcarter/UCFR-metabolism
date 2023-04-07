# plot parameter estimates from hierarchical models

library(tidyverse)
library(viridis)

mod_ests <- read_csv('data/model_fits/ar_PIcurve_parameter_ests.csv')
mod_ests <- read_csv('data/model_fits/ar1_linear_model_parameter_ests_ss.csv')
chains <- read_csv('data/model_fits/ar_PIcurve_chains.csv')
chains <- read_csv('data/model_fits/ar1_linear_model_chains_ss.csv') %>%
    bind_rows(chains)
mod_metrics <- read_csv('data/model_fits/ar_PIcurve_model_metrics_unconditioned.csv')
mod_metrics <- read_csv('data/model_fits/ar1_linear_model_metrics_ss_uncond.csv') %>%
    bind_rows(mod_metrics)

best_mods_by_cat <- mod_metrics %>%
    group_by(model)%>%
    filter(r2adj > max(r2adj) - 0.012) %>%
    mutate(across(c('r2adj', 'rmse'), round, digits = 2),
           across(starts_with(c('loo', 'waic')), round, digits = 0))


mod_ests %>%
    filter(!grepl('^beta', parameter) ,
           parameter != 'gamma_intercept') %>%
    ggplot( aes(parameter, mean, fill = units, col = units))+
    geom_boxplot(aes(ymin = X_2.5,  middle = X_50,
                      ymax = X_97.5)) +
    # ylim(0,2) +
    facet_wrap(.~biomass_vars)


mod_ests %>% filter(parameter %in% c('gamma_epil_chla', 'gamma_fila_chla'))


ch <- data.frame(parameter = c(rep('gamma_fila_chla',2),
                                   rep('gamma_epil_chla',2)),
                     chains = rep(0,4),
                     biomass_vars = rep("light",4)) %>%
    bind_rows(chains)

nl <- chains %>%
    filter(grepl('^gamma', parameter)|parameter == 'phi',
           parameter != 'gamma_intercept')%>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = case_when(parameter == "gamma_light" ~ 'light',
                                 parameter == "gamma_epil" ~"epil",
                                 parameter == "gamma_fila" ~ "fila",
                                 TRUE ~ parameter),
           parameter = factor(parameter,
                              levels=c("phi", "light",
                                       "fila", "epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ggtitle('GPP=f(bm)')+
    ylab("Parameter Estimate") +
    xlab("") +
    # ylim(-0.1,0.8)+
    # ylim(-0.5,5)+
    theme_bw()
chains <- read_csv('data/model_fits/hmodel_GPP_bm_linGAM_post_dists_for_plot.csv')
chains <- data.frame(parameter = c(rep('gamma_fila',2),
                                   rep('gamma_epil',2)),
                     chains = rep(0,4),
                     biomass_vars = rep(NA_character_,4),
                     units = rep(NA_character_,4)) %>%
    bind_rows(chains)

nllin <- chains %>%
    # filter(parameter != 'gamma_k600')%>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = case_when(parameter == "gamma_light" ~ 'light',
                                 parameter == "gamma_epil" ~"epil",
                                 parameter == "gamma_fila" ~ "fila",
                                 TRUE ~ parameter),
           parameter = factor(parameter,
                              levels=c("phi", "light",
                                       "fila", "epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ggtitle('GPP=f(bm)')+
    ylab("Parameter Estimate") +
    xlab("") +
    # ylim(-0.1,0.8)+
    # ylim(-0.5,5)+
    theme_bw()

chains <- read_csv('data/model_fits/hmodel_logGPP_bm_linGAM_post_dists_for_plot.csv')
chains <- data.frame(parameter = c(rep('gamma_fila',2),
                                   rep('gamma_epil',2)),
                     chains = rep(0,4),
                     biomass_vars = rep(NA_character_,4),
                     units = rep(NA_character_,4)) %>%
    bind_rows(chains)

log_plin <- chains %>%
    # filter(parameter != 'gamma_k600')%>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = case_when(parameter == "gamma_light" ~ 'light',
                                 parameter == "gamma_epil" ~"epil",
                                 parameter == "gamma_fila" ~ "fila",
                                 TRUE ~ parameter),
           parameter = factor(parameter,
                              levels=c("phi", "light",
                                       "fila", "epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ylab("") +
    ggtitle('log(GPP) = f(bm)')+
    xlab("") +
    ylim(-0.1,0.8)+
    # ylim(-0.5,5)+
    theme_bw()

chains <- read_csv('data/model_fits/hmodel_logGPP_logbm_linGAM_post_dists_for_plot.csv')
chains <- data.frame(parameter = c(rep('gamma_fila',2),
                                   rep('gamma_epil',2)),
                     chains = rep(0,4),
                     biomass_vars = rep(NA_character_,4),
                     units = rep(NA_character_,4)) %>%
    bind_rows(chains)

loglog_plin <- chains %>%
    # filter(parameter != 'gamma_k600')%>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = case_when(parameter == "gamma_light" ~ 'light',
                                 parameter == "gamma_epil" ~"epil",
                                 parameter == "gamma_fila" ~ "fila",
                                 TRUE ~ parameter),
           parameter = factor(parameter,
                              levels=c("phi", "light",
                                       "fila", "epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ylab("") +
    ggtitle("log(GPP) = fn(log(bm))")+
    ylim(-0.1,0.8)+
    xlab("") +
    theme_bw()
chains <- read_csv('data/model_fits/hmodel_logGPP_bm_logGAM_post_dists_for_plot.csv')
chains <- data.frame(parameter = c(rep('gamma_fila',2),
                                   rep('gamma_epil',2)),
                     chains = rep(0,4),
                     biomass_vars = rep(NA_character_,4),
                     units = rep(NA_character_,4)) %>%
    bind_rows(chains)

log_p <- chains %>%
    # filter(parameter != 'gamma_k600')%>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = case_when(parameter == "gamma_light" ~ 'light',
                                 parameter == "gamma_epil" ~"epil",
                                 parameter == "gamma_fila" ~ "fila",
                                 TRUE ~ parameter),
           parameter = factor(parameter,
                              levels=c("phi", "light",
                                       "fila", "epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ylab("") +
    ggtitle('log(GPP) = f(bm)')+
    xlab("") +
    ylim(-0.1,0.8)+
    # ylim(-0.5,5)+
    theme_bw()

chains <- read_csv('data/model_fits/hmodel_logGPP_logbm_logGAM_post_dists_for_plot.csv')
chains <- data.frame(parameter = c(rep('gamma_fila',2),
                                   rep('gamma_epil',2)),
                     chains = rep(0,4),
                     biomass_vars = rep(NA_character_,4),
                     units = rep(NA_character_,4)) %>%
    bind_rows(chains)

loglog_p <- chains %>%
    # filter(parameter != 'gamma_k600')%>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = case_when(parameter == "gamma_light" ~ 'light',
                                 parameter == "gamma_epil" ~"epil",
                                 parameter == "gamma_fila" ~ "fila",
                                 TRUE ~ parameter),
           parameter = factor(parameter,
                              levels=c("phi", "light",
                                       "fila", "epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ylab("") +
    ggtitle("log(GPP) = fn(log(bm))")+
    ylim(-0.1,0.8)+
    xlab("") +
    theme_bw()

ggpubr::ggarrange(nllin, log_plin, loglog_plin, nl, log_p, loglog_p, ncol = 3, nrow = 2,
                  common.legend = TRUE)
