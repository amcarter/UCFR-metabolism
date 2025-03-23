#############################################################
# Script to compare across different GAM fits to biomass data
#############################################################

library(tidyverse)

fit_linear <- read_csv('data/biomass_data/linear_gam_fits_biomass.csv') %>%
    mutate(fit = 'linear')
k_linear <- read_csv('data/biomass_data/linear_gam_smoothness_parameter_checks.csv') %>%
    mutate(fit = 'linear')
fit_log <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv') %>%
    mutate(fit = 'log')
k_log <- read_csv('data/biomass_data/log_gamma_gam_smoothness_parameter_checks.csv') %>%
    mutate(fit = 'log')
fit_inv <- read_csv('data/biomass_data/invgamma_gam_fits_biomass.csv') %>%
    mutate(fit = 'inverse')
k_inv <- read_csv('data/biomass_data/invgamma_gam_smoothness_parameter_checks.csv') %>%
    mutate(fit = 'inverse')
fit_idgam <- read_csv('data/biomass_data/idgamma_gam_fits_biomass.csv') %>%
    mutate(fit = 'identity')
k_idgam <- read_csv('data/biomass_data/idgamma_gam_smoothness_parameter_checks.csv') %>%
    mutate(fit = 'identity')

fit <- bind_rows(fit_linear, fit_log)
colnames(fit)
k <- bind_rows(k_linear, k_log, k_inv, k_idgam)

k %>% filter(k. < 2*edf | p.value < 0.05)
# these all look reasonable, some of the p values are significant, but
# none of the edfs are close to their limits (k.)

fit <- fit %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 names_to = c('biomass_type', 'units', 'stat'),
                 names_pattern = '(^[a-z]+)_([a-z0-9_]+)_([a-z]+$)',
                 values_to = 'value') %>%
    rename(model = fit) %>%
    pivot_wider(names_from = 'stat', values_from = 'value')

ep_afdm <- fit %>%
    filter(units == 'gm2', biomass_type == 'epil') %>%
    ggplot(aes(date, fit, col = model)) +
    geom_line() +
    geom_ribbon(aes(ymax= fit + se, ymin = fit - se, fill = model),
                col = NA, alpha = 0.3) +
    facet_grid(site~year, scales = 'free_x')+
    scale_y_log10()+
    ggtitle('Epilitheon, AFDM') +
    theme_bw()
ep_chla <- fit %>%
    filter(units == 'chla_mgm2', biomass_type == 'epil') %>%
    ggplot(aes(date, fit, col = model)) +
    geom_line() +
    geom_ribbon(aes(ymax= fit + se, ymin = fit - se, fill = model),
                col = NA, alpha = 0.3) +
    facet_grid(site~year, scales = 'free_x')+
    scale_y_log10()+
    ggtitle('Epilitheon, chl') +
    theme_bw()
fi_afdm <- fit %>%
    filter(units == 'gm2', biomass_type == 'fila') %>%
    ggplot(aes(date, fit, col = model)) +
    geom_line() +
    geom_ribbon(aes(ymax= fit + se, ymin = fit - se, fill = model),
                col = NA, alpha = 0.3) +
    facet_grid(site~year, scales = 'free_x')+
    scale_y_log10()+
    ggtitle('Filamentous, AFDM') +
    theme_bw()
fi_chla <- fit %>%
    filter(units == 'chla_mgm2', biomass_type == 'fila') %>%
    ggplot(aes(date, fit, col = model)) +
    geom_line() +
    geom_ribbon(aes(ymax= fit + se, ymin = fit - se, fill = model),
                col = NA, alpha = 0.3) +
    facet_grid(site~year, scales = 'free_x')+
    scale_y_log10(limits = c(0.3, 1000))+
    ggtitle('Filamentous, chl') +
    theme_bw()

ggpubr::ggarrange(ep_afdm, ep_chla, fi_afdm, fi_chla, nrow = 2, ncol = 2,
                  common.legend = TRUE)

ep_afdm
ep_chla
fi_afdm
fi_chla

# Conclusion: they're all pretty similar...
# the two with the best fits to the data are the inverse and log
# Let's compare those two

fit %>%
    select(-se) %>%
    filter(model %in% c('log', 'inverse'), biomass_type == 'fila') %>%
    pivot_wider(values_from = 'fit', names_from = 'model') %>%
    mutate(delta = inverse - log) %>%
    ggplot(aes(date, delta, col = units)) +
    geom_line() +
    facet_grid(site~year, scales = 'free_x')+
    ylim(-20,90)+
    geom_hline(yintercept = 0)

# okay, so they are similar, except the inverse model captures more
# of the dynamics in filamentous chlorophyll. We will go with that one







