# plot parameter estimates from hierarchical models

library(tidyverse)
library(viridis)

mod_ests <- read_csv( 'data/model_fits/hierarchical_model_parameters_fixedK.csv')
mod_ests <- read_csv( 'data/model_fits/hierarchical_model_parameters.csv')
chains <- read_csv('data/model_fits/hierarchical_model_posterior_distributions_for_plot_logbm_fixedK.csv')

mod_ests %>%
    filter(!grepl('^beta', parameter) ,
           parameter != 'intercept') %>%
    ggplot( aes(parameter, mean, fill = units, col = units))+
    geom_boxplot(aes(ymin = X2.5., lower = X25., middle = X50.,
                     upper = X75., ymax = X97.5.)) +
    # ylim(0,2) +
    facet_wrap(.~biomass_vars)


mod_ests %>% filter(parameter %in% c('gamma_epil', 'gamma_fila'))

plot(dd$epil_gm2_fit, dd$epil_chla_mgm2_fit)
abline(0,1)
plot(dd$fila_gm2_fit, dd$fila_chla_mgm2_fit)


chains <- read_csv('data/model_fits/hierarchical_model_posterior_distributions_for_plot_logbm_fixedK.csv')
chains <- data.frame(parameter = c(rep('gamma_fila',2),
                                   rep('gamma_epil',2)),
                     chains = rep(0,4),
                     biomass_vars = rep(NA_character_,4),
                     units = rep(NA_character_,4)) %>%
    bind_rows(chains)

lg_K <- chains %>%
    filter(parameter != 'gamma_k600')%>%
    filter(biomass_vars == 'fila_epil' | is.na(biomass_vars)) %>%
    mutate(model = case_when(units == 'chla_mgm2' ~ 'c Algal chl a',
                             units == 'gm2' ~ 'b Algal mass',
                             TRUE ~ 'a No biomass')) %>%
    mutate(parameter = factor(parameter,
                              levels=c("phi", "gamma_k600", "gamma_light",
                                       "gamma_fila", "gamma_epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ylab("Parameter Estimate") +
    xlab("") +
    ylim(-0.5,5)+
    theme_bw()

ggpubr::ggarrange(lin, lg, lin_K, lg_K, nrow = 2, ncol = 2,common.legend = TRUE)
row = c(2,1))
lin
