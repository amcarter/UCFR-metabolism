# Partition productivity into algal groups
# plot the role of algal assemblage in fluxes and turnovers

library(tidyverse)
library(lme4)


met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK.csv') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

q90 <- read_csv('data/quantile_PR_fits_summary.csv')

met <- left_join(met, select(q90, site, year, ARf), by = c('site', 'year')) %>%
    mutate(year = factor(year),
           NPP = GPP * (1-ARf),
           # NPP_globalARf = GPP * (1 + beta[2,1]),
           AR = -GPP*(ARf),
           HR = ER - AR) %>%
    select(-msgs.fit, -warnings, -errors, -K600, -DO_fit)

biogams <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass_epip.csv')
bio <- read_csv('data/biomass_data/biomass_working_data_summary_epip.csv') %>%
    mutate(epi_gm2 = epilithon_gm2 + epiphyton_gm2,
           epi_chla_mgm2 = epilithon_chla_mgm2 + epiphyton_chla_mgm2,
           filaepip_gm2 = epiphyton_gm2 + filamentous_gm2,
           filaepip_chla_mgm2 = epiphyton_chla_mgm2 + filamentous_chla_mgm2)

min_vals <- bio %>%
    select(epilithon_gm2, epilithon_chla_mgm2,
           epiphyton_gm2, epiphyton_chla_mgm2,
           epi_gm2, epi_chla_mgm2,
           filamentous_gm2, filamentous_chla_mgm2,
           filaepip_gm2, filaepip_chla_mgm2) %>%
    mutate(across(.fns = ~case_when(. == 0 ~ NA_real_,
                                    TRUE ~ .))) %>%
    pivot_longer(cols = everything(), names_to = 'biomass', values_to = 'value') %>%
    group_by(biomass) %>%
    summarize(min = quantile(value, 0.01, na.rm = T),
              q2 = quantile(value, 0.02, na.rm = T))


# Calculate average site year metrics:
NPP <- met %>%
    group_by(site, year) %>%
    summarize(across(any_of(c('GPP', 'ER', 'ARf', 'NPP',
                              'AR', 'HR')),
                     mean, na.rm = T))

# compare NPP to biomass Growth ####
NPP <- bio %>%
    mutate(year = factor(lubridate::year(date))) %>%
    group_by(site, year) %>%
    summarize(fila_chla_meas = mean(filamentous_chla_mgm2),
              filaepip_meas = mean(filaepip_chla_mgm2),
              epip_chla_meas = mean(epiphyton_chla_mgm2),
              epil_chla_meas = mean(epilithon_chla_mgm2),
              fila_meas = mean(filamentous_gm2),
              filaepip_meas = mean(filaepip_gm2),
              epip_meas = mean(epiphyton_gm2),
              epil_meas = mean(epilithon_gm2)) %>%
    left_join(NPP, by = c('site', 'year'))
NPP <- mutate(biogams, year = factor(year)) %>%
    group_by(site, year) %>%
    summarize(fila_chla = mean(fila_chla_mgm2_fit),
              filaepip_chla = mean(filaepip_chla_mgm2_fit),
              epip_chla = mean(epip_chla_mgm2_fit),
              epil_chla = mean(epil_chla_mgm2_fit),
              fila_gm = mean(fila_gm2_fit),
              filaepip_gm = mean(filaepip_gm2_fit),
              epip_gm = mean(epip_gm2_fit),
              epil_gm = mean(epil_gm2_fit),
              frac_fila_chla = fila_chla/(epil_chla + epip_chla + fila_chla),
              frac_fila_gmB = filaepip_gm/(epil_gm + filaepip_gm),
              frac_fila_gm = fila_gm/(epil_gm + epip_chla + fila_gm)) %>%
    left_join(NPP, by = c('site', 'year'))


npp <- NPP %>%
    select(-ends_with('_meas')) %>%
    mutate(biomass_C = (epil_gm + epip_gm + fila_gm)/2,
           biomass_chla = epil_chla + epip_chla + fila_chla,
           GPP_C = GPP * 12/32,
           ER_C = ER * 12/32,
           NPP_C = NPP * 12/32,
           # NPP_C_globalARf = NPP_globalARf * 12/32,
           days_bio = biomass_C/NPP_C) %>%
    ungroup() %>%
    mutate(fila_bloom = c(0,1,0,1,0,0,1,1,1,0,0,0),
           fila_bloom = case_when(fila_bloom == 0 ~ 'No Bloom',
                                  fila_bloom == 1 ~ 'Bloom'))


plot(density(npp$days_bio))

ggplot(npp, aes(frac_fila_gm, GPP))+
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
    select(frac_fila_gmB, GPP_C, ER_C,
           ARf, NPP_C, biomass_C,
           days_bio) %>% cor()
cor.test(npp$frac_fila_gm, npp$biomass_C, method = 'pearson')
cor.test(npp$frac_fila_gm, npp$days_bio, method = 'pearson')

corrplot::corrplot(ff_cor)
labs <- data.frame(response = c('A_GPP_C', 'B_ER_C', 'C_ARf', 'D_NPP_C',
                     'E_biomass_C', 'F_days_bio'),
                   corr = c(#'', '', '', '',
                            paste0('r = ', round(ff_cor[1,2], 2)),
                            paste0('r = ', round(ff_cor[1,3], 2)),
                            paste0('r = ', round(ff_cor[1,4], 2)),
                            paste0('r = ', round(ff_cor[1,5], 2)),
                            paste0('r = ', round(ff_cor[1,6], 2), '0'),
                            paste0('r = ', round(ff_cor[1,7], 2))))


png('figures/algal_assemblage_relationships.png',
    width = 4.5, height = 5, units = 'in', res = 300)
    npp %>%
        select(site, year, fila_bloom, fila_gm, frac_fila_gmB, A_GPP_C = GPP_C, B_ER_C = ER_C,
               C_ARf = ARf, D_NPP_C = NPP_C, E_biomass_C = biomass_C,
               F_days_bio = days_bio) %>%
        mutate(C_ARf = 100*C_ARf,
               fila_bloom = factor(fila_bloom, levels = c('No Bloom', 'Bloom'))) %>%
        pivot_longer(cols = any_of(c('A_GPP_C', 'B_ER_C', 'C_ARf', 'D_NPP_C',
                                     'E_biomass_C', 'F_days_bio')),
                     names_to = 'response',
                     values_to = 'value') %>%
        ggplot(aes(frac_fila_gmB, value)) +
        # geom_point(aes(col = fila_gm), size = 1.6) +
        # scale_color_continuous(expression('Filamentous Algae (g m'^{-2}*')'))+
        geom_point(aes(shape = fila_bloom,
                   col = fila_bloom), size = 1.4) +
        scale_color_manual("", values = c('#111011', '#7FB43A'))+
        scale_shape_manual("", values = c(19,17))+
        facet_wrap(response~., scales = 'free_y',
                   strip.position = 'left',
                   nrow = 3,
                   labeller = as_labeller(c(A_GPP_C = 'GPP~(g~C~m^{-2}~d^{-1})',
                                            B_ER_C = 'ER~(g~C~m^{-2}~d^{-1})',
                                            C_ARf = 'AR~Fraction~("%")',
                                            D_NPP_C = 'NPP~(g~C~m^{-2}~d^{-1})',
                                            E_biomass_C = 'Biomass~(g~C~m^{-2})',
                                            F_days_bio = 'Residence~Time~(d)'),
                                          default = label_parsed)) +
        labs(y = NULL, x = 'Filamentous fraction of total biomass (%)')+
        # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, col = 'grey50')+
        scale_x_continuous(expand = c(0.04, 0.04)) +
        scale_y_continuous(expand = c(0.1, 0.5)) +
        geom_text(data = labs, aes(label = corr, x = -Inf, y = Inf),
                  size = 3, hjust = -0.25, vjust = 1.8)+
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


# Growth rate calculations ####
light <- read_csv('data/site_data/daily_modeled_light_all_sites.csv')
bm_met <- select(biogams, site, date, epil_gm2_fit, epip_gm2_fit, epi_gm2_fit,
                 fila_gm2_fit, filaepip_gm2_fit, epil_chla_mgm2_fit,
                 epip_chla_mgm2_fit, epi_chla_mgm2_fit,
                 fila_chla_mgm2_fit, filaepip_chla_mgm2_fit) %>%
    rename_with(~gsub('_fit', '', .x)) %>%
    left_join(select(ungroup(met), site, date, year, GPP, ER, ARf, NPP),
              by = c('site', 'date')) %>%
    left_join(select(light, site, date, PAR_surface)) %>%
    mutate(light = PAR_surface/max(PAR_surface)) %>%
    mutate(#light = PAR_surface,
           NPP_light = NPP/light) %>%
    filter(!is.na(GPP))

ggplot(bm_met, aes(log(epil_chla_mgm2+1), NPP_light, col = fila_chla_mgm2))+
    geom_point()
ggplot(bm_met, aes(log(fila_chla_mgm2+1),  NPP_light, col = site))+
    geom_point()
ggplot(bm_met, aes(log(epip_chla_mgm2+1), NPP_light, col = epi_chla_mgm2))+
    geom_point()
ggplot(bm_met, aes(log(epi_chla_mgm2)+log(fila_chla_mgm2), log(NPP_light), col = site))+
    geom_point()

mod1 <- lm(NPP_light ~ 0 + epip_chla_mgm2 + epil_chla_mgm2 +
               fila_chla_mgm2, data = bm_met)
modmix1 <- lmer(NPP_light ~ 0 + (0 + epil_chla_mgm2|site) +
                    (0 + epip_chla_mgm2|site) + (0 + fila_chla_mgm2|site), data = bm_met)
ranef(modmix1)
mod2 <- lm(NPP_light ~ 0  + filaepip_chla_mgm2 + epil_chla_mgm2, data = bm_met)
mod2b <- lm(NPP_light ~ 0  + log(filaepip_chla_mgm2+1) + log(epil_chla_mgm2+1), data = bm_met)
modmix2 <- lmer(NPP_light ~ 0  + filaepip_chla_mgm2 + epil_chla_mgm2 +
                    (0 + filaepip_chla_mgm2|site) + (0 + epil_chla_mgm2|site),
                data = bm_met)
mod3 <- lm(NPP_light ~ 0  + fila_chla_mgm2 + epi_chla_mgm2, data = bm_met)
modmix3 <- lmer(NPP_light ~ 0  + fila_chla_mgm2 + epi_chla_mgm2 +
                    (0 + fila_chla_mgm2|site) + (0 + epi_chla_mgm2|site),
                data = bm_met)

summary(mod1)
summary(mod2)
summary(mod3)

AIC(mod1)
AIC(mod2)
AIC(mod2b)
AIC(mod3)
plot(mod3)

# modmix <- lme4::lmer(NPP/light ~ 0 + epil_chla_mgm2 + epip_chla_mgm2 + fila_chla_mgm2 +
#                          (0+epil_chla_mgm2|site) + (0+epip_chla_mgm2|site) +
#                          (0+fila_chla_mgm2|site),
#                      data = bm_met,
#                      REML = FALSE,
#                      control = lmerControl(optimizer ="Nelder_Mead"))

# library(brms)
# brm.f1 <- brmsformula(NPP_light ~ 0 + epi_chla_mgm2 + fila_chla_mgm2)
#
# get_prior(brm.f1, data = bm_met)
# prior.f1 <- c(prior(gamma(0.5, 1), class = b, coef = epi_chla_mgm2),# lb = 0),
#               prior(gamma(0.5, 1), class = b, coef = fila_chla_mgm2),# lb = 0),
#               prior(normal(0, 0.1), class = sigma))
#
# modmix1 <- brms::brm(NPP_light ~ 0 + epi_chla_mgm2 + fila_chla_mgm2,
#                      data = bm_met)
#
#
# modmix <- lme4::lmer(NPP/light ~ 0 + epil_chla_mgm2 + epip_chla_mgm2 + fila_chla_mgm2 +
#                          (0+epil_chla_mgm2|site) + (0+epip_chla_mgm2|site) +
#                          (0+fila_chla_mgm2|site),
#                      data = bm_met,
#                      REML = FALSE,
#                      control = lmerControl(optimizer ="Nelder_Mead"))
# summary(mod2a)
# summary(mod3)
# summary(modmix)
# ranef(modmix)

bm_met %>%
    mutate(preds = predict(mod2, newdata = bm_met)) %>%
    ggplot(aes(NPP_light, preds, col = site)) +
    geom_point() + geom_abline(intercept = 0, slope = 1)

coefs <- summary(mod1)$coefficients[,1]
coefs2 <- summary(mod2)$coefficients[,1]
coefs3 <- summary(mod3)$coefficients[,1]
# coefs <- fixef(modmix)

bm_met_epil_fila <- bm_met %>%
    mutate(fila_prod_gCd = (coefs2[1] * filaepip_chla_mgm2*light) * 12/32 ,
           fila_turnover = case_when(filaepip_gm2 > min_vals$min[min_vals$biomass == 'filaepip_gm2'] ~
                                         fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epil_prod_gCd = (coefs2[2] * epil_chla_mgm2*light) * 12/32 ,
           epil_turnover = epil_gm2/2/epil_prod_gCd)

bm_met_epi_fila <- bm_met %>%
    mutate(fila_prod_gCd = (coefs3[1] * fila_chla_mgm2*light) * 12/32 ,
           fila_turnover = case_when(fila_gm2 > min_vals$min[min_vals$biomass == 'filamentous_gm2'] ~
                                         fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epi_prod_gCd = (coefs3[2] * epi_chla_mgm2*light) * 12/32 ,
           epi_turnover = case_when(epi_gm2 > min_vals$min[min_vals$biomass == 'epi_gm2'] ~
                                        epi_gm2/2/epi_prod_gCd,
                                     TRUE ~ NA_real_))
bm_met_epil_epip_fila <- bm_met %>%
    mutate(fila_prod_gCd = (coefs[3] * fila_chla_mgm2*light) * 12/32 ,
           fila_turnover = case_when(fila_gm2 > min_vals$min[min_vals$biomass == 'filamentous_gm2'] ~
                                         fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epip_prod_gCd = (coefs[1] * epip_chla_mgm2*light) * 12/32 ,
           epip_turnover = case_when(epip_gm2 > min_vals$min[min_vals$biomass == 'epiphyton_gm2'] ~
                                        epip_gm2/2/epip_prod_gCd,
                                     TRUE ~ NA_real_),
           epil_prod_gCd = (coefs[2] * epil_chla_mgm2*light) * 12/32 ,
           epil_turnover = epil_gm2/2/epil_prod_gCd)


# compare across different model fits:
plot(density(bm_met_epil_fila$epil_turnover))
lines(density(bm_met_epi_fila$epi_turnover))
lines(density(bm_met_epil_epip_fila$epil_turnover))

plot(density(bm_met_epil_fila$fila_turnover, na.rm = T))
lines(density(bm_met_epi_fila$fila_turnover, na.rm = T))
lines(density(bm_met_epil_epip_fila$fila_turnover, na.rm = T))


bm_met_epil_fila %>%
    filter(year == 2021) %>%
    select(-year, -ARf, -PAR_surface) %>%
    write_csv('data/biomass_data/2021_turnovers_for_Rafa.csv')

bm_met_epi_fila %>%
    filter(year == 2021) %>%
    select(-year, -ARf, -PAR_surface) %>%
    write_csv('data/biomass_data/2021_turnovers_for_Rafa_epi_grouped.csv')

bm_met_epil_epip_fila %>%
    filter(year == 2021) %>%
    select(-year, -ARf, -PAR_surface) %>%
    write_csv('data/biomass_data/2021_turnovers_for_Rafa_epip_separate.csv')

# coefs <- summary(mod2a)$coefficients[,1]
coefs <- fixef(modmix)

bm_met2 <- bm_met %>%
    mutate(fila_prod_gCd = (coefs[2] * fila_chla_mgm2*light) * 12/32 ,
           fila_turnover = case_when(fila_gm2 > min_fila_gm2 ~ fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epi_prod_gCd = (coefs[1] * epi_chla_mgm2*light) * 12/32 ,
           epi_turnover = epi_gm2/2/epi_prod_gCd)

bm_met2 %>%
    filter(year == 2021) %>%
    select(-year, -ARf, -PAR_surface) %>%
    write_csv('data/biomass_data/2021_turnovers_for_Rafa_epip_epil_grouped.csv')

# bm_met <- bm_met %>%
#     mutate(fila_prod_gCd = (0.02394 * fila_gm2*light)*14/32 ,
#            fila_turnover = case_when(fila_prod_gCd > 0.04 ~ fila_gm2/2/fila_prod_gCd,
#                                      TRUE ~ NA_real_),
#            epil_prod_gCd = (0.435 * epil_gm2*light)*14/32 ,
#            epil_turnover = epil_gm2/2/epil_prod_gCd)


bm_met1 %>% group_by(site, year) %>%
    summarize(across(ends_with(c('gCd', 'turnover')), ~mean(., na.rm = T))) %>%
    summary()
summary(bm_met2)
bm_fila <- filter(bm_met, fila_gm2 > min_fila_gm2)

p1 <- bm_fila %>%
    pivot_longer(cols = ends_with('turnover'),
                 values_to = 'turnover', names_to = 'Biomass',
                 names_pattern = '(epil|fila)_.*') %>%
    ggplot(aes(x=turnover, group = Biomass, fill = Biomass)) +
    geom_density(adjust=1.5, alpha=.4) +
    scale_fill_manual(values = c('#1B9EC9', '#97BB43'))+
    xlim(0, 40)+
    ylab('Density')+
    xlab('Turnover time (d)')+
    theme_classic()+
    theme(legend.position = 'none',
          panel.border = element_rect(fill = NA))
p1 <- p1 + annotate(geom = 'text', x = 12, y = 0.25,
                    label="Epilithic", col = '#1B9EC9')
p1 <- p1 + annotate(geom = 'text', x = 32, y = 0.06,
                    label="Filamentous", col = '#97BB43')
p2 <- ggplot(npp, aes(fila_gm/2/(biomass_C)*100, days_bio))+
    xlab('Filamentous fraction of \ntotal biomass (%)')+
    ylab('Turnover time (d)')+
    geom_point(aes(pch = factor(fila_bloom), col = factor(fila_bloom)),
               size = 1.6)+
    theme_classic()+
    theme(legend.title = element_blank(),
          legend.position = c(0.22, 0.85),
          panel.border = element_rect(fill = NA))
png('figures/turnover_time_by_composition.png',
    width = 6.5, height = 3.5, units = 'in', res = 300)
ggpubr::ggarrange(p1, p2, common.legend = TRUE)
dev.off()

ggplot(bm_met, aes(date, fila_prod_gCd))+
    geom_line(col = 'forestgreen') +
    geom_line(aes(y = epil_prod_gCd), col = 'sienna')+
    facet_grid(site~year, scales = 'free_x')
# png('figures/NPP_vs_standing_crop_all_sites.png',
#     width = 6.5, height = 6.5,
#     res = 300, units = 'in')
#     bm_met %>%
#         mutate(doy = as.numeric(format(date, '%j')),
#                Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
#                site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
#         select(site, Date, year, fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
#         pivot_longer(cols = starts_with(c('fila', 'epil')),
#                      values_to = 'value', names_to = c('biomass', 'measure'),
#                      names_pattern = '(fila|epil)_(.*)') %>%
#         pivot_wider(names_from = measure, values_from = value) %>%
#         mutate(gCm2 = gm2/2)%>%
#     ggplot(aes(prod_gCd, gCm2, col = Date, pch = biomass))+
#         geom_point() +
#         facet_grid(site~year)+
#         theme_classic()+
#         scale_shape_manual('Biomass', values = c(1, 19))+
#         ylab(expression(paste('Biomass standing crop (gC ', m^-2, ')')))+
#         xlab(expression(paste('Biomass production (gC ', m^-2, d^-1, ')')))+
#         theme(panel.border = element_rect(fill = NA),
#               panel.spacing = unit(0, units = 'in'))
# dev.off()
# png('figures/NPP_vs_standing_crop_all_sites.png',
#     width = 6.5, height = 6.5,
#     res = 300, units = 'in')
coeff = 125
ann_text2 <- read_csv('data/metabolism/auto_sites_figure_labels.csv') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           year = factor(year),
           Date = case_when(year == 2020 ~ as.Date('2020-10-20'),
                            TRUE ~ as.Date('2020-10-14')),
           prod_gCd = rep(0.45, 12), gCm2 = rep(1, 12))

p <- bm_met %>%
    mutate(doy = as.numeric(format(date, '%j')),
           Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
           site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
    mutate(fila_gm2 = case_when(fila_gm2 < min_fila_gm2 ~ NA_real_,
                                TRUE ~ fila_gm2)) %>%
    select(site, Date, year, light,
           fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
    pivot_longer(cols = starts_with(c('fila', 'epil')),
                 values_to = 'value', names_to = c('biomass', 'measure'),
                 names_pattern = '(fila|epil)_(.*)') %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(Light = factor(rep(" ", 2*nrow(bm_met)), levels = c("1"," ")),
           gCm2 = gm2/2,
           biomass = case_when(biomass == 'fila' ~ 'Filamentous',
                               biomass == 'epil'~'Epilithic'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
    ggplot(aes(Date, prod_gCd/gCm2, col = biomass))+
    geom_area(aes(Date, gCm2/coeff, fill = biomass), color = NA, alpha = 0.4)+
    geom_line(size = 1.2)+
    geom_line(aes(y = light/2.5, lty = Light), col = 'grey20') +
    facet_grid(site~year, scales = 'free_x', ) +
    geom_text(data = ann_text2, aes(label = trophic), col = 'black') +
    scale_y_continuous(
        name = expression(paste('Production rate (', d^-1, ')')),
        sec.axis = sec_axis(~.*coeff,
                            name = expression(paste('Biomass (g C', m^-2, ')')))
    )+
    scale_fill_manual('Biomass standing stock',
                      values = c('#97BB43', '#1B9EC9')) +
    scale_color_manual('Biomass production rate',
                       values = c('#97BB43', '#1B9EC9')) +
    scale_linetype_manual('Light', values = c(2)) +
    theme_classic()+
    xlab('Date')+
    theme(strip.text.y = element_blank(),
          panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'in'),
          # legend.spacing.x = unit(1, 'cm'),
          # legend.justification = 'top')
          # legend.box.margin = margin(0, 5, 0, 5, "cm"),
          legend.position = 'top')+
    guides(color = guide_legend(title.position = "top", title.hjust = 0,
                                order = 1),
           fill = guide_legend(title.position = "top", title.hjust = 0,
                               order = 2),
           linetype = guide_legend(title.position = 'top', title.hjust = 0,
                                   order = 3))
png('figures/biomass_prod_and_turnover.png', width = 6, height =8,
    units = 'in', res = 300)
p
dev.off()

ann_text <- data.frame(Date = rep(as.Date('2020-07-16'), 6),
                       prod_gCd = rep(0.48, 6),
                       gCm2 = rep(1, 6),
                       year = rep(2020,6),
                       biomass = factor(rep('Epilithic', 6),
                                        levels = c('Filamentous','Epilithic')),
                       site = factor(c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'),
                                     levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))


p3 <- p + geom_text(data = ann_text, aes(label = site), col = 'black')

bm_sum <-  bm_met %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
    group_by(site, year) %>%
    summarize(n = n(),
              fila_Biomass = max(fila_gm2)/2,
              epil_Biomass = max(epil_gm2)/2,
              frac_fila = mean(fila_gm2)/mean(fila_gm2+epil_gm2),
              fila_cumprod = sum(fila_prod_gCd)/n*100,
              epil_cumprod = sum(epil_prod_gCd)/n*100,
              cumprod_gC = fila_cumprod + epil_cumprod) %>%
    mutate(bloom = factor(case_when(fila_Biomass > 30 ~ 'Bloom',
                                    TRUE ~ 'No Bloom')))

p2 <- ggplot(bm_sum, aes(frac_fila, cumprod_gC))+
    xlab('Filamentous fraction of \ntotal biomass (%)')+
    ylab('Turnover time (d)')+
    geom_point(aes(pch = factor(bloom)), size = 1.6, col = 'black')+
    theme_classic()+
    theme(legend.title = element_blank(),
          legend.position = c(0.22, 0.85),
          panel.border = element_rect(fill = NA))
p4 <- bm_sum %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 values_to = 'value',
                 names_to = c('biomass', 'measure'),
                 names_pattern = '(fila|epil)_(Biomass|cumprod)') %>%
    mutate(measure = case_when(measure == 'cumprod' ~ 'Cumulative \nProduction',
                               measure == 'Biomass' ~ 'Maximum \nBiomass'),
           measure = factor(measure, levels = c('Maximum \nBiomass',
                                                'Cumulative \nProduction')),
           biomass = case_when(biomass == 'epil' ~ 'Epilithic',
                               biomass == 'fila' ~ 'Filamentous'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
    ggplot(aes(x = measure, y = value, fill = biomass)) +
    geom_boxplot(outlier.shape=NA, alpha = 0.4)+
    geom_point(aes(pch = bloom, group = biomass),
               position = position_jitterdodge(jitter.width = 0.25))+
    # geom_jitter(color="black", alpha=0.9) +
    scale_fill_manual('',
                      values = c('#97BB43', '#1B9EC9')) +
    ylab(expression(paste('Mass (g C ', m^-2, ')')))+
    xlab('')+
    theme_classic()+
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.spacing.y = unit(-0.1, "cm"),
          legend.position = c(0.24, 0.79),
          legend.background = element_rect(fill = NA))
png('figures/biomass_cumulative_and_turnover.png', width = 6.5, height = 3.5,
    units = 'in', res = 300)
ggpubr::ggarrange(p4, p1, ncol = 2,
                  labels = c('A', 'B'),
                  align = 'h')

dev.off()

bm_sum %>% group_by(bloom) %>%
    summarize(across(where(is.numeric), .fns = c(mean = ~mean(.), sd = ~sd(.))))

p5 <- ggplot(data.frame(a = 1, b = 1), aes(a,b)) +
    geom_point(col = 'white') +
    xlab('') + ylab('')+
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())

png('figures/biomass_prod_and_turnover_3panel.png', width = 8, height = 8,
    units = 'in', res = 300)
ggpubr::ggarrange(p3,
                  ggpubr::ggarrange(p2, p1, p4, p5, ncol = 1,
                                    labels = c('B', 'C', 'D', ''),
                                    heights = c(5.1, 5, 5.2, 1.3),
                                    align = 'v'),
                  ncol = 2, labels = 'A', #label.y = 0.917,
                  widths = c(2,1))
dev.off()
png('figures/biomass_turnover_frac.png', width = 3.5, height = 3,
    units = 'in', res = 300)
p2
dev.off()


p2 <- bm_met %>%
    filter(year == 2021, site == 'GC') %>%

    mutate(doy = as.numeric(format(date, '%j')),
           Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
           site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
    mutate(fila_gm2 = case_when(fila_gm2 < min_fila_gm2 ~ NA_real_,
                                TRUE ~ fila_gm2)) %>%
    select(site, Date, year, light,
           fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
    pivot_longer(cols = starts_with(c('fila', 'epil')),
                 values_to = 'value', names_to = c('biomass', 'measure'),
                 names_pattern = '(fila|epil)_(.*)') %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(gCm2 = gm2/2,
           biomass = case_when(biomass == 'fila' ~ 'Filamentous',
                               biomass == 'epil'~'Epilithic'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
    filter(biomass == 'Epilithic') %>%
    ggplot(aes(Date, gCm2, fill = biomass))+
    geom_area(color = NA, alpha = 0.7)+
    ylab(expression(paste('Biomass (g C', m^-2, ')')))+
    scale_fill_manual('Biomass standing stock',
                      values = c( '#1B9EC9')) +
    theme_classic()+
    xlab('Date')+
    ylim(0,47)+
    theme(legend.position = 'none')
# p1 <-
bm_met %>%
    filter(year == 2021, site == 'GC') %>%
    mutate(doy = as.numeric(format(date, '%j')),
           Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
           fila_cum = cumsum(fila_prod_gCd),
           epil_cum = cumsum(epil_prod_gCd)) %>%
    ggplot(aes(Date, fila_cum))+
    geom_line(color = '#97BB43', size = 1.5)+
    geom_line(aes(y = epil_cum), color = '#1B9EC9', size = 1.5)+
    ylab(expression(paste('Biomass production (g C', m^-2, d^-1, ')')))+
    scale_color_manual('Biomass production rate',
                       values = c('#97BB43', '#1B9EC9')) +
    theme_classic()+
    xlab('Date')




bm_met %>%
    group_by(site, year) %>%
    mutate(fila_rate = fila_prod_gCd/(fila_gm2/2),
           epil_rate = epil_prod_gCd/(epil_gm2/2)) %>%
    filter(fila_gm2 >= min_fila_gm2) %>%
    summary()
bm_met %>%
    group_by(site, year) %>%
    filter((year == 2020 & site %in% c('GR', 'GC'))|
               (year == 2021 & site %in% c('GC', 'BN','BM')))%>%
    summarize(n = n(),
              fila_Biomass = max(fila_gm2)/2,
              epil_Biomass = max(epil_gm2)/2,
              fila_prod = mean(fila_prod_gCd),
              epil_prod = mean(epil_prod_gCd),
              fila_cumprod = sum(fila_prod_gCd),
              epil_cumprod = sum(epil_prod_gCd)) %>%
    mutate(n_epil = epil_cumprod/epil_Biomass,
           n_fila = fila_cumprod/fila_Biomass)%>%
    summary()

bm_met %>%
    filter(!(site %in% c('DL', 'PL') ))%>%
    ggplot(aes(date, fila_turnover)) +
    geom_line() +
    facet_grid(site~year, scale = 'free_x') + ylim(0,100)

ggplot(aes(Date, NPP/(fila_gm2 + epil_gm2),
           col = fila_gm2/(epil_gm2 + fila_gm2)))+
    geom_point() +
    facet_grid(site~year, scales = 'free')+
    geom_line(aes(y = light))+
    theme_classic()+
    ylab(expression(paste('Biomass standing crop (gC ', m^-2, ')')))+
    xlab(expression(paste('Biomass production (gC ', m^-2, d^-1, ')')))+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'in'))
dev.off()


bm_met %>%
    mutate(fila_gm2 = case_when(fila_gm2 < min_fila_gm2 ~ NA_real_,
                                TRUE ~ fila_gm2)) %>%
    mutate(fila_mgchlagC = fila_chla_mgm2/(fila_gm2/2),
           epil_mgchlagC = epil_chla_mgm2/(epil_gm2/2)) %>%
    ggplot(aes(date, fila_mgchlagC))+
    geom_point() +
    geom_point(aes(y = epil_mgchlagC), col = 2) +
    facet_grid(site~year, scales = 'free_x')

dd <- bm_met %>%
    select(site, year, date, epil_gm2, fila_gm2, GPP, ER, NPP, ARf,
           fila_prod_gCd, epil_prod_gCd) %>%
    mutate(GPP = GPP * 12/32,
           ER = ER *12/32,
           NPP = NPP *12/32,
           AR = GPP*ARf,
           HR = ER+AR,
           fila_gC = fila_gm2/2,
           epil_gC = epil_gm2/2)

dd$fila_diff = c(0, diff(dd$fila_gC) - dd$fila_prod_gCd[1:(nrow(dd)-1)])
dd$epil_diff = c(0, diff(dd$epil_gC) - dd$epil_prod_gCd[1:(nrow(dd)-1)])

ids <- dd %>%
    mutate(r_number = row_number()) %>%
    group_by(site, year) %>%
    summarise(index = min(r_number))
dd$epil_diff[ids$index] <- dd$fila_diff[ids$index] <- NA

dd %>%
    mutate(algal_loss = fila_diff + epil_diff - HR) %>%
    pivot_longer(cols = ends_with('diff'),
                 values_to = 'loss', names_to = 'algae') %>%
    mutate(loss = case_when(loss>0 ~ 0,
                     TRUE ~ loss))%>%
    ggplot(aes(date, loss, fill = algae))+
    geom_area()+
    geom_line(aes(y = HR))+
    facet_grid(site~year, scales = 'free_x') +
    theme_classic()

ggplot(data, aes(x=time, y=value, fill=group)) +
    geom_area()

filter(dd, site == 'DL', year == 2020)
fila_prod_gCd/(fila_prod_gCd + epil_prod_gCd),
           epil_AR = AR * epil_prod_gCd/(fila_prod_gCd + epil_prod_gCd))



