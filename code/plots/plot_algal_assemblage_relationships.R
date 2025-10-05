# Partition productivity into algal groups
# plot the role of algal assemblage in fluxes and turnovers

library(tidyverse)
library(lme4)


met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK.csv') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')))

q90 <- read_csv('data/quantile_PR_fits_summary_brms.csv') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site))

# met <- left_join(met, select(q90, site, year, ARf), by = c('site', 'year')) %>%
met <- left_join(met, q90, by = c('site', 'year')) %>%
    rename(ARf = slope, ARf.lower = slope.lower, ARf.upper = slope.upper) %>%
    mutate(year = factor(year),
           NPP = GPP * (1-ARf),
           # NPP_globalARf = GPP * (1 + beta[2,1]),
           AR = -GPP*(ARf),
           HR = ER - AR) %>%
    select(-msgs.fit, -warnings, -errors, -K600, -DO_fit)

biogams <- read_csv('data/biomass_data/log_gamma_brms_gam_fits_biomass.csv') %>%
    mutate(frac_fila = if_else(is.na(frac_fila),
                               filamentous_gm2/(filamentous_gm2 + epilithon_gm2),
                               frac_fila))
bio <- read_csv('data/biomass_data/biomass_working_data_summary.csv') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site))
min_epil_gm2 <- min(bio$epilithon_gm2)
min_fila_gm2 <- filter(bio, filamentous_gm2>0) %>% select(filamentous_gm2) %>% min()


# Calculate average site year metrics:
NPP <- met %>%
    group_by(site, year) %>%
    summarize(across(starts_with(c('GPP', 'ER', 'ARf', 'NPP',
                              'AR', 'HR')),
                     mean, na.rm = T))


# compare NPP to biomass Growth ####
NPP <- bio %>%
    mutate(year = factor(lubridate::year(date))) %>%
    group_by(site, year) %>%
    summarize(fila_chla_meas = median(filamentous_chla_mgm2),
              epil_chla_meas = median(epilithon_chla_mgm2),
              fila_meas = median(filamentous_gm2),
              epil_meas = median(epilithon_gm2)) %>%
    left_join(NPP, by = c('site', 'year'))
NPP <- mutate(biogams, year = factor(lubridate::year(date))) %>%
    group_by(site, year) %>%
    summarize(fila_chla = median(filamentous_chla_mgm2),
              epil_chla = median(epilithon_chla_mgm2),
              upper_fila_chla = median(filamentous_chla_upper),
              lower_fila_chla = median(filamentous_chla_lower),
              upper_epil_chla = median(epilithon_chla_upper),
              lower_epil_chla = median(epilithon_chla_lower),
              fila_gm = median(filamentous_gm2),
              epil_gm = median(epilithon_gm2),
              upper_fila_gm = median(filamentous_gm2_upper),
              lower_fila_gm = median(filamentous_gm2_lower),
              upper_epil_gm = median(epilithon_gm2_upper),
              lower_epil_gm = median(epilithon_gm2_lower),
              frac_fila = median(frac_fila),
              frac_fila_upper = median(frac_fila_upper, na.rm = T),
              frac_fila_lower = median(frac_fila_lower, na.rm = T)) %>%
    left_join(NPP, by = c('site', 'year'))

# ggplot(NPP, aes(epil_meas, epil_gm, col = site)) +
#     geom_point(size = 2)
# ggplot(NPP, aes(fila_meas, fila_gm, col = site)) +
#     geom_point(size = 2)
# ggplot(NPP, aes(epil_gm, ARf, col = site)) +
#     geom_point(size = 2)
# ggplot(NPP, aes(fila_chla, ARf, col = year)) +
#     geom_point(size = 2)
# ggplot(NPP, aes((fila_chla + epil_chla)/(fila_gm + epil_gm), ARf, col = year)) +
#     geom_point(size = 2)
# ggplot(NPP, aes((fila_chla)/(fila_gm), ARf, col = site)) +
#     geom_point(size = 2)
# ggplot(NPP, aes(fila_chla/(fila_chla + epil_chla), ARf, col = year)) +
#     geom_point(size = 2)
# ggplot(NPP, aes(epil_chla/(fila_chla + epil_chla), ARf, col = year)) +
#     geom_point(size = 2)

npp <- NPP %>%
    select(-ends_with('_meas')) %>%
    mutate(biomass_C = (epil_gm + fila_gm)/2,
           biomass_C_upper = (upper_epil_gm + upper_fila_gm)/2,
           biomass_C_lower = (lower_epil_gm + lower_fila_gm)/2,
           biomass_chla = epil_chla + fila_chla,
           GPP_C = GPP * 12/32,
           GPP_C_upper = GPP.upper * 12/32,
           GPP_C_lower = GPP.lower * 12/32,
           ER_C = ER * 12/32,
           ER_C_upper = ER.upper * 12/32,
           ER_C_lower = ER.lower * 12/32,
           NPP_C = NPP * 12/32,
           # NPP_C_globalARf = NPP_globalARf * 12/32,
           days_bio = biomass_C/NPP_C,
           days_bio_upper = biomass_C_upper/NPP_C,
           days_bio_lower = biomass_C_lower/NPP_C
           ) %>%
    ungroup() %>%
    mutate(fila_bloom = c(0,1,0,1,0,0,1,1,1,0,0,0),
           fila_bloom = case_when(fila_bloom == 0 ~ 'No Bloom',
                                  fila_bloom == 1 ~ 'Bloom'))

ggplot(npp, aes(NPP_C, fila_gm, col = factor(year)))+
    geom_point(size = 2)
summary(npp)

mean(npp$days_bio)
plot(density(npp$days_bio))

ggplot(npp, aes(frac_fila, GPP))+
    geom_point(aes(col = fila_bloom), size = 1.5)
npp %>%
    pivot_longer(cols = c('GPP', 'ER', 'ARf', 'NPP_C',
                          'biomass_C', 'biomass_chla', 'days_bio'),
                 names_to = 'response',
                 values_to = 'value') %>%
    pivot_longer(cols = starts_with('frac'),
                 names_to = 'fila_units',
                 names_prefix = 'frac_fila_',
                 values_to = 'frac_fila') %>%
    select(site, year, fila_bloom, frac_fila, fila_units, response, value) %>%
    ggplot(aes(frac_fila, value, col = fila_bloom)) +
    geom_point() +
    facet_grid(response~fila_units, scales = 'free')

ff_cor <- npp %>%
    select(frac_fila, GPP_C, ER_C,
           ARf, NPP_C, biomass_C,
           days_bio) %>% cor()
cor.test(npp$frac_fila, npp$biomass_C, method = 'pearson')
cor.test(npp$frac_fila, npp$days_bio, method = 'pearson')

corrplot::corrplot(ff_cor)
labs <- data.frame(response = c('A_GPP_C', 'B_ER_C', 'C_ARf', 'D_NPP_C',
                     'E_biomass_C', 'F_days_bio'),
                   corr = c(#'', '', '', '',
                            paste0('r = ', round(ff_cor[1,2], 2)),# '0'),
                            paste0('r = ', round(ff_cor[1,3], 2)),
                            paste0('r = ', round(ff_cor[1,4], 2)),
                            paste0('r = ', round(ff_cor[1,5], 2)),
                            paste0('r = ', round(ff_cor[1,6], 2)),
                            paste0('r = ', round(ff_cor[1,7], 2))))



png('figures/algal_assemblage_relationships_meds.png',
    width = 4.5, height = 5, units = 'in', res = 300)

extra_row <- data.frame(fila_bloom = 'Bloom', frac_fila = -1,
                        A_GPP_C.mean = 0, B_ER_C.mean = -0.8, C_ARf.mean = 50,
                        D_NPP_C.mean = 0, E_biomass_C.mean = 0, F_days_bio.mean = 0)

    npp %>%
        select(site, year, fila_bloom, fila_gm,
               frac_fila, frac_fila_upper, frac_fila_lower,
               A_GPP_C.mean = GPP_C, A_GPP_C.upper = GPP_C_upper,
               A_GPP_C.lower = GPP_C_lower, B_ER_C.mean = ER_C,
               B_ER_C.upper = ER_C_upper, B_ER_C.lower = ER_C_lower,
               C_ARf.mean = ARf, C_ARf.lower = ARf.lower, C_ARf.upper = ARf.upper,
               D_NPP_C.mean = NPP_C, E_biomass_C.mean = biomass_C,
               E_biomass_C.upper = biomass_C_upper, E_biomass_C.lower = biomass_C_lower,
               F_days_bio.mean = days_bio, F_days_bio.upper = days_bio_upper,
               F_days_bio.lower = days_bio_lower) %>%
        mutate(C_ARf.mean = 100*C_ARf.mean, C_ARf.upper = C_ARf.upper * 100,
               C_ARf.lower = 100*C_ARf.lower,
               fila_bloom = factor(fila_bloom, levels = c('No Bloom', 'Bloom'))) %>%
        bind_rows(extra_row) %>%
        pivot_longer(cols = starts_with(c('A_GPP_C', 'B_ER_C', 'C_ARf', 'D_NPP_C',
                                     'E_biomass_C', 'F_days_bio')),
                     names_to = 'response',
                     values_to = 'value') %>%
        separate_wider_delim(response, ".", names = c("response", "stat")) %>%
        pivot_wider(names_from = stat, values_from = value) %>%
        ggplot(aes(frac_fila, mean)) +
        # geom_point(aes(col = fila_gm), size = 1.6) +
        # scale_color_continuous(expression('Filamentous Algae (g m'^{-2}*')'))+
        # geom_errorbar(aes(xmin = frac_fila_lower, xmax = frac_fila_upper,
        #                   col = fila_bloom), alpha = 0.4) +
        # geom_errorbar(aes(ymin = lower, ymax = upper,
        #                   col = fila_bloom), alpha = 0.4) +
        geom_point(aes(shape = fila_bloom,
                       col = fila_bloom), size = 1.4) +
        scale_color_manual("", values = c( '#7FB43A','#111011'))+
        scale_shape_manual("", values = c(17,19))+
        facet_wrap(response~., scales = 'free_y',
                   strip.position = 'left',
                   nrow = 3,
                   labeller = as_labeller(c(A_GPP_C = 'GPP~(g~C~m^{-2}~d^{-1})',
                                            B_ER_C = 'ER~(g~C~m^{-2}~d^{-1})',
                                            C_ARf = 'AR~Fraction~("%")',
                                            D_NPP_C = 'P[N]~(g~C~m^{-2}~d^{-1})',
                                            E_biomass_C = 'Biomass~(g~C~m^{-2})',
                                            F_days_bio = 'Turnover~Time~(d)'),
                                          default = label_parsed)) +
        labs(y = NULL, x = 'Filamentous fraction of total biomass')+
        # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, col = 'grey50')+
        # xlim(0,0.91) +
        xlim(0,0.82) +
        # scale_x_continuous(expand = c(0.04, 0.04)) +
        scale_y_continuous(expand =  expand_scale(mult = c(.08, 0.2),
                                                  add = c(0, 0.2))) +
        geom_text(data = labs, aes(label = corr, x = -Inf, y = Inf),
                  size = 3, hjust = -0.25, vjust = 1.7)+
        theme_classic()+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(color = "grey50", fill = NA),
              legend.title = element_text(size = 9, vjust = 0.8),
              axis.line = element_line(color = "grey50"),     # Color of axis lines
              axis.text = element_text(color = "grey50"),     # Color of axis labels
              axis.ticks = element_line(color = "grey50"),    # Color of axis ticks
              strip.text = element_text(color = "grey10"),
              # legend.title = element_blank(),
              legend.position = 'top')

dev.off()

tiff('figures/algal_assemblage_relationships_meds.tiff',
    width = 8.5, height = 9.45, units = 'cm', res = 600)

extra_row <- data.frame(fila_bloom = 'Bloom', frac_fila = -1,
                        A_GPP_C.mean = 0, B_ER_C.mean = -0.8, C_ARf.mean = 50,
                        D_NPP_C.mean = 0, E_biomass_C.mean = 0, F_days_bio.mean = 0)

    npp %>%
        select(site, year, fila_bloom, fila_gm,
               frac_fila, frac_fila_upper, frac_fila_lower,
               A_GPP_C.mean = GPP_C, A_GPP_C.upper = GPP_C_upper,
               A_GPP_C.lower = GPP_C_lower, B_ER_C.mean = ER_C,
               B_ER_C.upper = ER_C_upper, B_ER_C.lower = ER_C_lower,
               C_ARf.mean = ARf, C_ARf.lower = ARf.lower, C_ARf.upper = ARf.upper,
               D_NPP_C.mean = NPP_C, E_biomass_C.mean = biomass_C,
               E_biomass_C.upper = biomass_C_upper, E_biomass_C.lower = biomass_C_lower,
               F_days_bio.mean = days_bio, F_days_bio.upper = days_bio_upper,
               F_days_bio.lower = days_bio_lower) %>%
        mutate(C_ARf.mean = 100*C_ARf.mean, C_ARf.upper = C_ARf.upper * 100,
               C_ARf.lower = 100*C_ARf.lower,
               fila_bloom = factor(fila_bloom, levels = c('No Bloom', 'Bloom'))) %>%
        bind_rows(extra_row) %>%
        pivot_longer(cols = starts_with(c('A_GPP_C', 'B_ER_C', 'C_ARf', 'D_NPP_C',
                                     'E_biomass_C', 'F_days_bio')),
                     names_to = 'response',
                     values_to = 'value') %>%
        separate_wider_delim(response, ".", names = c("response", "stat")) %>%
        pivot_wider(names_from = stat, values_from = value) %>%
        ggplot(aes(frac_fila, mean)) +
        # geom_point(aes(col = fila_gm), size = 1.6) +
        # scale_color_continuous(expression('Filamentous Algae (g m'^{-2}*')'))+
        # geom_errorbar(aes(xmin = frac_fila_lower, xmax = frac_fila_upper,
        #                   col = fila_bloom), alpha = 0.4) +
        # geom_errorbar(aes(ymin = lower, ymax = upper,
        #                   col = fila_bloom), alpha = 0.4) +
        geom_point(aes(shape = fila_bloom,
                       col = fila_bloom), size = 1.4) +
        scale_color_manual("", values = c( '#7FB43A','#111011'))+
        scale_shape_manual("", values = c(17,19))+
        facet_wrap(response~., scales = 'free_y',
                   strip.position = 'left',
                   nrow = 3,
                   labeller = as_labeller(c(A_GPP_C = 'GPP~(g~C~m^{-2}~d^{-1})',
                                            B_ER_C = 'ER~(g~C~m^{-2}~d^{-1})',
                                            C_ARf = 'AR~Fraction~("%")',
                                            D_NPP_C = 'P[N]~(g~C~m^{-2}~d^{-1})',
                                            E_biomass_C = 'Biomass~(g~C~m^{-2})',
                                            F_days_bio = 'Turnover~Time~(d)'),
                                          default = label_parsed)) +
        labs(y = NULL, x = 'Filamentous fraction of total biomass')+
        # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, col = 'grey50')+
        # xlim(0,0.91) +
        xlim(0,0.82) +
        # scale_x_continuous(expand = c(0.04, 0.04)) +
        scale_y_continuous(expand =  expand_scale(mult = c(.08, 0.2),
                                                  add = c(0, 0.2))) +
        geom_text(data = labs, aes(label = corr, x = -Inf, y = Inf),
                  size = 3, hjust = -0.25, vjust = 1.7)+
        theme_classic()+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(color = "grey50", fill = NA),
              legend.title = element_text(size = 9, vjust = 0.8),
              axis.line = element_line(color = "grey50"),     # Color of axis lines
              axis.text = element_text(color = "grey50"),     # Color of axis labels
              axis.ticks = element_line(color = "grey50"),    # Color of axis ticks
              strip.text = element_text(color = "grey10"),
              # legend.title = element_blank(),
              legend.position = 'top')

dev.off()


png('figures/algal_assemblage_relationships_meds_error.png',
    width = 4.5, height = 5, units = 'in', res = 300)

extra_row <- data.frame(fila_bloom = 'Bloom', frac_fila = -1,
                        A_GPP_C.mean = 0, B_ER_C.mean = -0.8, C_ARf.mean = 50,
                        D_NPP_C.mean = 0, E_biomass_C.mean = 0, F_days_bio.mean = 0)

    npp %>%
        select(site, year, fila_bloom, fila_gm,
               frac_fila, frac_fila_upper, frac_fila_lower,
               A_GPP_C.mean = GPP_C, A_GPP_C.upper = GPP_C_upper,
               A_GPP_C.lower = GPP_C_lower, B_ER_C.mean = ER_C,
               B_ER_C.upper = ER_C_upper, B_ER_C.lower = ER_C_lower,
               C_ARf.mean = ARf, C_ARf.lower = ARf.lower, C_ARf.upper = ARf.upper,
               D_NPP_C.mean = NPP_C, E_biomass_C.mean = biomass_C,
               E_biomass_C.upper = biomass_C_upper, E_biomass_C.lower = biomass_C_lower,
               F_days_bio.mean = days_bio, F_days_bio.upper = days_bio_upper,
               F_days_bio.lower = days_bio_lower) %>%
        mutate(C_ARf.mean = 100*C_ARf.mean, C_ARf.upper = C_ARf.upper * 100,
               C_ARf.lower = 100*C_ARf.lower,
               fila_bloom = factor(fila_bloom, levels = c('No Bloom', 'Bloom')),
               D_NPP_C.upper = D_NPP_C.mean + 0.48, D_NPP_C.lower = D_NPP_C.mean - 0.48) %>%
        bind_rows(extra_row) %>%
        pivot_longer(cols = starts_with(c('A_GPP_C', 'B_ER_C', 'C_ARf', 'D_NPP_C',
                                     'E_biomass_C', 'F_days_bio')),
                     names_to = 'response',
                     values_to = 'value') %>%
        separate_wider_delim(response, ".", names = c("response", "stat")) %>%
        pivot_wider(names_from = stat, values_from = value) %>%
        ggplot(aes(frac_fila, mean)) +
        geom_errorbar(aes(xmin = frac_fila_lower, xmax = frac_fila_upper,
                          col = fila_bloom), alpha = 0.34) +
        geom_errorbar(aes(ymin = lower, ymax = upper,
                          col = fila_bloom), alpha = 0.34) +
        geom_point(aes(shape = fila_bloom,
                       col = fila_bloom), size = 1.4) +
        scale_color_manual("", values = c( '#7FB43A','#111011'))+
        scale_shape_manual("", values = c(17,19))+
        facet_wrap(response~., scales = 'free_y',
                   strip.position = 'left',
                   nrow = 3,
                   labeller = as_labeller(c(A_GPP_C = 'GPP~(g~C~m^{-2}~d^{-1})',
                                            B_ER_C = 'ER~(g~C~m^{-2}~d^{-1})',
                                            C_ARf = 'AR~Fraction~("%")',
                                            D_NPP_C = 'P[N]~(g~C~m^{-2}~d^{-1})',
                                            E_biomass_C = 'Biomass~(g~C~m^{-2})',
                                            F_days_bio = 'Turnover~Time~(d)'),
                                          default = label_parsed)) +
        labs(y = NULL, x = 'Filamentous fraction of total biomass')+
        xlim(0,0.91) +
        scale_y_continuous(expand =  expand_scale(#mult = c(.08),
                                                  add = 0.9)) +
        geom_text(data = labs, aes(label = corr, x = -Inf, y = Inf),
                  size = 3, hjust = -0.25, vjust = 1.7)+
        theme_classic()+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(color = "grey50", fill = NA),
              legend.title = element_text(size = 9, vjust = 0.8),
              axis.line = element_line(color = "grey50"),     # Color of axis lines
              axis.text = element_text(color = "grey50"),     # Color of axis labels
              axis.ticks = element_line(color = "grey50"),    # Color of axis ticks
              strip.text = element_text(color = "grey10"),
              legend.position = 'top')

dev.off()
tiff('figures/algal_assemblage_relationships_meds_error.tiff',
    width = 8.5, height = 9.4, units = 'cm', res = 600)

extra_row <- data.frame(fila_bloom = 'Bloom', frac_fila = -1,
                        A_GPP_C.mean = 0, B_ER_C.mean = -0.8, C_ARf.mean = 50,
                        D_NPP_C.mean = 0, E_biomass_C.mean = 0, F_days_bio.mean = 0)

    npp %>%
        select(site, year, fila_bloom, fila_gm,
               frac_fila, frac_fila_upper, frac_fila_lower,
               A_GPP_C.mean = GPP_C, A_GPP_C.upper = GPP_C_upper,
               A_GPP_C.lower = GPP_C_lower, B_ER_C.mean = ER_C,
               B_ER_C.upper = ER_C_upper, B_ER_C.lower = ER_C_lower,
               C_ARf.mean = ARf, C_ARf.lower = ARf.lower, C_ARf.upper = ARf.upper,
               D_NPP_C.mean = NPP_C, E_biomass_C.mean = biomass_C,
               E_biomass_C.upper = biomass_C_upper, E_biomass_C.lower = biomass_C_lower,
               F_days_bio.mean = days_bio, F_days_bio.upper = days_bio_upper,
               F_days_bio.lower = days_bio_lower) %>%
        mutate(C_ARf.mean = 100*C_ARf.mean, C_ARf.upper = C_ARf.upper * 100,
               C_ARf.lower = 100*C_ARf.lower,
               fila_bloom = factor(fila_bloom, levels = c('No Bloom', 'Bloom')),
               D_NPP_C.upper = D_NPP_C.mean + 0.48, D_NPP_C.lower = D_NPP_C.mean - 0.48) %>%
        bind_rows(extra_row) %>%
        pivot_longer(cols = starts_with(c('A_GPP_C', 'B_ER_C', 'C_ARf', 'D_NPP_C',
                                     'E_biomass_C', 'F_days_bio')),
                     names_to = 'response',
                     values_to = 'value') %>%
        separate_wider_delim(response, ".", names = c("response", "stat")) %>%
        pivot_wider(names_from = stat, values_from = value) %>%
        ggplot(aes(frac_fila, mean)) +
        geom_errorbar(aes(xmin = frac_fila_lower, xmax = frac_fila_upper,
                          col = fila_bloom), alpha = 0.34) +
        geom_errorbar(aes(ymin = lower, ymax = upper,
                          col = fila_bloom), alpha = 0.34) +
        geom_point(aes(shape = fila_bloom,
                       col = fila_bloom), size = 1.4) +
        scale_color_manual("", values = c( '#7FB43A','#111011'))+
        scale_shape_manual("", values = c(17,19))+
        facet_wrap(response~., scales = 'free_y',
                   strip.position = 'left',
                   nrow = 3,
                   labeller = as_labeller(c(A_GPP_C = 'GPP~(g~C~m^{-2}~d^{-1})',
                                            B_ER_C = 'ER~(g~C~m^{-2}~d^{-1})',
                                            C_ARf = 'AR~Fraction~("%")',
                                            D_NPP_C = 'P[N]~(g~C~m^{-2}~d^{-1})',
                                            E_biomass_C = 'Biomass~(g~C~m^{-2})',
                                            F_days_bio = 'Turnover~Time~(d)'),
                                          default = label_parsed)) +
        labs(y = NULL, x = 'Filamentous fraction of total biomass')+
        xlim(0,0.91) +
        scale_y_continuous(expand =  expand_scale(#mult = c(.08),
                                                  add = 0.9)) +
        geom_text(data = labs, aes(label = corr, x = -Inf, y = Inf),
                  size = 2.5, hjust = -0.25, vjust = 1.7)+
        theme_classic()+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(color = "grey50", fill = NA),
              legend.title = element_text(size = 7.5, vjust = 1),
              legend.text = element_text(size = 7.5),
              axis.line = element_line(color = "grey50"),     # Color of axis lines
              axis.text = element_text(color = "grey50", size = 8),     # Color of axis labels
              axis.title = element_text( size = 9),     # Color of axis labels
              axis.ticks = element_line(color = "grey50"),    # Color of axis ticks
              strip.text = element_text(color = "grey10", size = 7.5),
              legend.position = 'top', vjust = 1)

dev.off()
