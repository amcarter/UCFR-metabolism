# fit GAMS to biomass data
# 10/2022

library(tidyverse)
library(mgcv)

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    filter(!is.na(site) & site != 'CR') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')),
           year = as.factor(lubridate::year(date)),
           site_year = paste(site, year, sep = '_'))

glimpse(biomass)


# look at the data
biomass %>%
    select(date, doy, year, site, ends_with('_gm2')) %>%
    pivot_longer(cols = -c('date', 'doy', 'year', 'site'),
                 names_to = 'bm_category',
                 values_to = 'gm2') %>%
    ggplot(aes(doy, gm2, col = year)) +
    geom_point() +
    facet_grid(site~bm_category, scales = 'free_y')


# find 2%quantile ####
# add the 2% quantile to each of the biomass categories in order to
# fit log transformed GAMs
colnames(biomass)
min_vals <- biomass %>%
    select(epilithon_gm2, epilithon_chla_mgm2,
           filamentous_gm2, filamentous_chla_mgm2,
           fila_macro_gm2) %>%
    mutate(across(.fns = ~case_when(. == 0 ~ NA_real_,
                                    TRUE ~ .))) %>%
    pivot_longer(cols = everything(), names_to = 'biomass', values_to = 'value') %>%
    group_by(biomass) %>%
    summarize(min = quantile(value, 0.01, na.rm = T),
              q2 = quantile(value, 0.02, na.rm = T))

dates <- c(seq(min(biomass$date), max(biomass$date[biomass$year == 2020]), by = 'day'),
           seq(min(biomass$date[biomass$year == 2021]), max(biomass$date), by = 'day'))


# biomass <- biomass %>%
#     mutate(epilithon_gm2 = epilithon_gm2 + min_vals$min[min_vals$biomass == 'epilithon_gm2'],
#            epilithon_chla_mgm2 = epilithon_chla_mgm2 +
#                min_vals$min[min_vals$biomass == 'epilithon_chla_mgm2'],
#            filamentous_gm2 = filamentous_gm2 + min_vals$min[min_vals$biomass == 'filamentous_gm2'],
#            filamentous_chla_mgm2 = filamentous_chla_mgm2 +
#                min_vals$min[min_vals$biomass == 'filamentous_chla_mgm2'],
#            fila_macro_gm2 = fila_macro_gm2 + min_vals$min[min_vals$biomass == 'fila_macro_gm2'])
# fit gams
preds <- data.frame(date = dates,
                    doy = as.numeric(format(dates, '%j')),
                    year = lubridate::year(dates))%>%
    mutate(year = as.factor(year))
s_preds <- data.frame()
for(site in unique(biomass$site)){
    preds$site = site
    s_preds <- bind_rows(s_preds, preds)
}
s_preds$site_year <- paste(s_preds$site, s_preds$year, sep = "_")

link_fn = 'log' # gamma link function, log or inverse?
delta = 0.5
par(mfrow = c(2,2))

hist(biomass$epilithon_gm2)
hist(log(biomass$epilithon_gm2 + min_vals$min[min_vals$biomass == 'epilithon_gm2'] ))
hist(biomass$filamentous_gm2)
hist(log(biomass$filamentous_gm2 + min_vals$min[min_vals$biomass == 'filamentous_gm2'] ))

# try it for all sites
biomass$site_year <- factor(biomass$site_year)
biomass$year <- factor(biomass$year)
biomass$site <- factor(biomass$site)
fg2_gamma <- gam(epilithon_gm2 +
                     min_vals$min[min_vals$biomass == 'epilithon_gm2'] ~ s(doy) +
                     s(doy, site_year, bs = 'fs'),
                 data = biomass, method = 'REML',
                 family = Gamma(link = link_fn))
AIC(fg2_gamma) # s(doy) = 2965.128; s(doy, year) = 2969.329, vs 2970.138
summary(fg2_gamma) #42.8%
gam.check(fg2_gamma)
sqrt(mean((fitted(fg2_gamma) - biomass$epilithon_gm2)^2)) # 8.68

# pp_lin <- mgcv::predict.gam(fg2, s_preds_lin, se.fit = TRUE)
pp_gamma <- mgcv::predict.gam(fg2_gamma, s_preds, type = 'response',
                              se.fit = TRUE)

# s_preds_lin <- mutate(s_preds_lin,
#                       epil_gm2_fit = c(pp_lin$fit),
#                       epil_gm2_se = c(pp_lin$se.fit))
s_preds <- mutate(s_preds,
                  epil_gm2_fit = c(pp_gamma$fit) -
                      min_vals$min[min_vals$biomass == 'epilithon_gm2'],
                  epil_gm2_se = c(pp_gamma$se.fit))

# filamentous
# fg2_fila <- gam(filamentous_gm2 ~ s(doy) +
#                     s(doy, site_year, bs = 'fs'),
#                 data = biomass, method = 'REML', family = 'gaussian')

# fg2_fila_gamma <- gam(filamentous_gm2 +
#                           min_vals$min[min_vals$biomass == 'filamentous_gm2'] ~
#                           s(doy, k = 20) +
#                           s(doy, site_year, bs = 'fs', k = 20),
#                       data = biomass, method = 'REML',
#                       family = Gamma(link = link_fn))
# Presence/absence model
biomass$present <- as.numeric(biomass$filamentous_gm2 > 0)
fg2_fila_binom <- gam(present ~ s(doy, site_year, bs = 'fs', k = 4),
            data = biomass, family = binomial)

AIC(fg2_fila_binom) # 315 vs 315
gam.check(fg2_fila_binom)
summary(fg2_fila_binom) #62.8%
sqrt(mean((fitted(fg2_fila_binom) - biomass$present)^2)) #0.29
# Gamma model on positives only
fg2_fila_gamma <- gam(filamentous_gm2 ~ s(doy, year, bs = 'fs') +
                           s(doy, site_year, bs = 'fs', k = 10),
            data = subset(biomass, filamentous_gm2 > 0),
            family = Gamma(link = "log"), method = 'REML')
AIC(fg2_fila_gamma) # s(doy) = 2392; s(doy, year) = 2390.851 cs 2388.598
sqrt(mean((fitted(fg2_fila_gamma) - biomass$filamentous_gm2[biomass$filamentous_gm2>0])^2)) # 55.7
gam.check(fg2_fila_gamma)
summary(fg2_fila_gamma) #55.3%
# biomass_sub$residuals <- residuals(fg2_fila_gamma)
# ggplot(biomass_sub, aes(doy, residuals, col = site_year)) +
#     geom_point() +
#     geom_line()

# pp_lin <- mgcv::predict.gam(fg2_fila, s_preds_lin, se.fit = TRUE)
pp_binom <- mgcv::predict.gam(fg2_fila_binom, newdata = s_preds,
                          type = "response", se.fit = TRUE)
pp_gamma1 <- mgcv::predict.gam(fg2_fila_gamma, type = 'response',
                              s_preds, se.fit = TRUE)
pp_gamma1$fit <- if_else(is.na(pp_gamma1$fit), 1, pp_gamma1$fit)
pp_gamma1$se.fit <- if_else(is.na(pp_gamma1$se.fit), 1, pp_gamma1$se.fit)
pp_gamma <- data.frame(fit = pp_binom$fit * pp_gamma1$fit,
                       se.fit = sqrt((pp_binom$fit^2 * pp_gamma1$se.fit^2) +
                                         (pp_gamma1$fit^2 * pp_binom$se.fit^2) +
                                         (pp_binom$se.fit^2 * pp_gamma1$se.fit^2)))

s_preds %>%
    mutate(pred_binom = pp_binom$fit) %>%
    ggplot(aes(doy, pred_binom)) +
    geom_line() +
    geom_point(data = biomass, aes(y = present), col = "green")+
    facet_grid(site~year)

s_preds <- mutate(s_preds,
                  fila_gm2_fit = c(pp_gamma$fit),
                  fila_gm2_se = c(pp_gamma$se.fit))


ggplot(s_preds, aes(doy, fila_gm2_fit)) +
    geom_line() +
    geom_point(data = biomass, aes(y = filamentous_gm2), col = "green")+
    facet_grid(site~year)




# chlorophyll
fg2_chla_gamma <- gam(epilithon_chla_mgm2 +
                          min_vals$min[min_vals$biomass == 'epilithon_chla_mgm2']~
                          s(doy, year, bs = "fs") +
                          s(doy, site_year, bs = 'fs'),
                      data = biomass, method = 'REML',
                      family = Gamma(link = link_fn))
AIC(fg2_chla_gamma) # s(doy) = 3820; s(doy, year) = 3817.509, vs 3815.476
summary(fg2_chla_gamma) #48.9%
sqrt(mean((fitted(fg2_chla_gamma) - biomass$epilithon_chla_mgm2)^2)) # 18.9

# fg2_chla <- gam(epilithon_chla_mgm2 ~ s(doy) +
#                   s(doy, site_year, bs = 'fs'),
#               data = biomass, method = 'REML', family = 'gaussian')
#
# gam.check(fg2_chla)
gam.check(fg2_chla_gamma)

# AIC(fg2_chla_gamma)
# AIC(fg2_chla)
#
pp_gamma <- mgcv::predict.gam(fg2_chla_gamma, s_preds,
                              type ='response', se.fit = TRUE)
# pp_lin <- mgcv::predict.gam(fg2_chla, s_preds_lin, se.fit = TRUE)

# s_preds_lin <- mutate(s_preds_lin,
#                       epil_chla_mgm2_fit = c(pp_lin$fit),
#                       epil_chla_mgm2_se = c(pp_lin$se.fit))
s_preds <- mutate(s_preds,
                  epil_chla_mgm2_fit = c(pp_gamma$fit) -
                      min_vals$min[min_vals$biomass == 'epilithon_chla_mgm2'],
                  epil_chla_mgm2_se = c(pp_gamma$se.fit))

# filamentous
biomass$present <- as.numeric(biomass$filamentous_chla_mgm2 > 0)
fg2_fila_chla_binom <- gam(present ~ s(doy, site_year, bs = 'fs', k = 4),
                       data = biomass, family = binomial)
AIC(fg2_fila_chla_binom)
gam.check(fg2_fila_chla_binom)
summary(fg2_fila_chla_binom) #59.5%
sqrt(mean((fitted(fg2_fila_chla_binom) - biomass$present)^2)) #0.30

# Gamma model on positives only
fg2_fila_chla_gamma <- gam(filamentous_chla_mgm2 ~ s(doy, year, bs = "fs") +
                           s(doy, site_year, bs = 'fs'),
                       data = subset(biomass, filamentous_chla_mgm2 > 0),
                       family = Gamma(link = "log"), method = 'REML')
AIC(fg2_fila_chla_gamma) # s(doy) = 2714 s(doy, year) = 2714.772 vs 2719.691
gam.check(fg2_fila_chla_gamma)
summary(fg2_fila_chla_gamma) # 49.9%
sqrt(mean((fitted(fg2_fila_chla_gamma) -
               biomass$filamentous_chla_mgm2[biomass$filamentous_chla_mgm2 >0])^2)) # 82.7

# biomass_sub$residuals <- residuals(fg2_fila_chla_gamma)
# ggplot(biomass_sub, aes(doy, residuals, col = site_year)) +
#     geom_point() +
#     geom_line()

pp_binom <- mgcv::predict.gam(fg2_fila_chla_binom, newdata = s_preds,
                          type = "response", se.fit = TRUE)
pp_gamma1 <- mgcv::predict.gam(fg2_fila_chla_gamma, type = 'response',
                               s_preds, se.fit = TRUE)
pp_gamma1$fit <- if_else(is.na(pp_gamma1$fit), 1, pp_gamma1$fit)
pp_gamma1$se.fit <- if_else(is.na(pp_gamma1$se.fit), 1, pp_gamma1$se.fit)
pp_gamma <- data.frame(fit = pp_binom$fit * pp_gamma1$fit,
                       se.fit = sqrt((pp_binom$fit^2 * pp_gamma1$se.fit^2) +
                                         (pp_gamma1$fit^2 * pp_binom$se.fit^2) +
                                         (pp_binom$se.fit^2 * pp_gamma1$se.fit^2)))

s_preds <- mutate(s_preds,
                  fila_chla_mgm2_fit = c(pp_gamma$fit),
                  fila_chla_mgm2_se = c(pp_gamma$se.fit))

s_preds %>%
    mutate(pred_binom = pp_binom$fit) %>%
    ggplot(aes(doy, pred_binom)) +
    geom_line() +
    geom_point(data = biomass, aes(y = present), col = "green")+
    facet_grid(site~year)



# Calculate fit metrics for the smoothness parameter
# gamma model
k1 <- k.check(fg2_gamma) %>% data.frame() %>%
    mutate(model = 'biofilm AFDM',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k1) <- NULL
k2 <- k.check(fg2_fila_gamma) %>% data.frame() %>%
    mutate(model = 'filamentous AFDM',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k2) <- NULL
k3 <- k.check(fg2_chla_gamma) %>% data.frame() %>%
    mutate(model = 'biofilm chla',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k3) <- NULL
k4 <- k.check(fg2_fila_chla_gamma) %>% data.frame() %>%
    mutate(model = 'filamentous chla',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k4) <- NULL
k_check_gamma <- bind_rows(k1, k2, k3, k4)

# plot diagnostics
png('figures/SI/biomass_loggammaGAMs_diagnostics.png', width = 7.5, height = 7.5,
    units = 'in', res = 300)
    par(mfrow = c(4,4),
        mar = c(3,4,2,2),
        oma = c(0,3,2,0))
    gam.check(fg2_gamma)
    mtext(expression(paste('epilithon g',m^-2)), 2, line = 45)
    gam.check(fg2_fila_gamma)
    mtext(expression(paste('filamentous g',m^-2)), 2, line = 45)
    gam.check(fg2_chla_gamma)
    mtext(expression(paste('epilithon chla mg',m^-2)), 2, line = 45)
    gam.check(fg2_fila_chla_gamma)
    mtext(expression(paste('filamentous chla mg',m^-2)), 2, line = 45)
    par(mfrow = c(1,1), new = T)
    mtext('Model fit metrics for Biomass GAMS', line = 2.75)
dev.off()
# png('figures/biomass_linGAMs_diagnostics.png', width = 7.5, height = 7.5,
#     units = 'in', res = 300)
#     par(mfrow = c(4,4),
#         mar = c(3,4,2,2),
#         oma = c(0,3,1,0))
#     gam.check(fg2)
#     mtext(expression(paste('epilithon g',m^-2)), 2, line = 45)
#     gam.check(fg2_fila)
#     mtext(expression(paste('filamentous g',m^-2)), 2, line = 45)
#     gam.check(fg2_chla)
#     mtext(expression(paste('epilithon chla mg',m^-2)), 2, line = 45)
#     gam.check(fg2_fila_chla)
#     mtext(expression(paste('filamentous chla mg',m^-2)), 2, line = 45)
# dev.off()

# qq = as_tibble(lapply(s_preds_lin, c)) %>% select(-site_year)
# write_csv(qq, 'data/biomass_data/linear_gam_fits_biomass.csv')
# write_csv(k_check, 'data/biomass_data/linear_gam_smoothness_parameter_checks.csv')

qq = as_tibble(lapply(s_preds, c)) %>% select(-site_year) %>%
    mutate(across(starts_with(c('epil', 'fila')),
                              ~case_when(. < 0 ~ 0,
                                         TRUE ~ .)))
write_csv(qq, 'data/biomass_data/log_gamma_gam_fits_biomass.csv')
write_csv(k_check_gamma, 'data/biomass_data/log_gamma_gam_smoothness_parameter_checks.csv')

qq <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv')
# plot GAMS
qq <- qq %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN'))) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
             names_to = c('biomass_type', 'units', 'stat'),
             names_pattern = '([a-z]+)_([a-z0-9_]+)_([a-z]+)',
             values_to = 'value') %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    group_by(site, biomass_type, units) %>%
    mutate(across(c('fit', 'se'), ~zoo::na.approx(., x = date, na.rm = F)))# %>%
    # mutate(fit = case_when(se > fit ~ NA_real_,
    #                        se > 500 ~ NA_real_,
    #                        TRUE ~ fit),
    #        se = case_when(is.na(fit)~NA_real_,
    #                       TRUE ~ se))
meas <- select(biomass, date, site, sample,
       epil_gm2_meas = epilithon_gm2,
       fila_gm2_meas = filamentous_gm2,
       epil_chla_mgm2_meas = epilithon_chla_mgm2,
       fila_chla_mgm2_meas = filamentous_chla_mgm2) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 names_to = c('biomass_type', 'units', 'stat'),
                 names_pattern = '([a-z]+)_([a-z0-9_]+)_([a-z]+)',
                 values_to = 'value') %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    mutate(year = lubridate::year(date))
meas_chl <- filter(meas, units == 'chla_mgm2')
meas_mass <- filter(meas, units == 'gm2')

# mm <- qq %>% filter(units == 'gm2') %>%
#     mutate(fit = fit + x,
#            fit_high = fit + se,
#            fit_low = fit - se,
#            fit_low = case_when(fit_low < 0.3 ~ 0.3,
#                                TRUE ~ fit_low)) %>%
#     ggplot(aes(date, fit, col = biomass_type)) +
#     geom_line()+
#     geom_ribbon(aes(ymax = fit_high, ymin = fit_low,
#                     fill = biomass_type), alpha = 0.4, color = NA)+
#     geom_point(data = meas_mass, aes(date, meas, col = biomass_type))+
#     facet_grid(site~year, scales = 'free_x') +
#     scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
#     scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
#     scale_y_log10(limits = c(0.3, 600))+
#     xlab('Date') +
#     ylab(expression('Algal Standing Crop (AFDM g '~ m^-2*')')) +
#     theme_bw()
# cc <- qq %>%
#     mutate(fit = fit + x,
#            fit_high = fit + se,
#            fit_low = fit - se,
#            fit_low = case_when(fit_low < 0.3 ~ 0.3,
#                                TRUE ~ fit_low)) %>%
#     filter(units == 'chla_mgm2') %>%
#     ggplot(aes(date, fit, col = biomass_type)) +
#     geom_line()+
#     geom_ribbon(aes(ymax = fit_high, ymin = fit_low,
#                     fill = biomass_type), alpha = 0.4, color = NA)+
#     geom_point(data = meas_chl, aes(date, meas, col = biomass_type))+
#     facet_grid(site~year, scales = 'free_x') +
#     scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
#     scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
#     scale_y_log10(limits = c(0.3, 1000))+
#     xlab('Date') +
#     ylab(expression('Algal Standing Crop (mg chl a '~ m^-2*')')) +
#     theme_bw()
meas_mass2 <- mutate(meas_mass,
                     biomass_type = case_when(biomass_type == 'epil'~ 'Epilithon',
                                              biomass_type == 'fila' ~ 'Filamentous'),
                     meas = case_when(meas < delta ~ delta,
                                      TRUE ~ meas))
meas_chl2 <- mutate(meas_chl,
                    biomass_type = case_when(biomass_type == 'epil'~ 'Epilithon',
                                             biomass_type == 'fila' ~ 'Filamentous'),
                    meas = case_when(meas < delta ~ delta,
                                     TRUE ~ meas))
mm <- qq %>% filter(units == 'gm2') %>%
    mutate(biomass_type = case_when(biomass_type == 'epil'~ 'Epilithon',
                                    biomass_type == 'fila' ~ 'Filamentous'),
           fit = case_when(fit < delta ~ delta,
                           TRUE ~ fit),
           fit_high = fit + se,
           fit_low = fit - se,
           fit_low = case_when(fit_low < delta ~ delta,
                               TRUE ~ fit_low)) %>%
    ggplot(aes(date, fit, col = biomass_type)) +
    geom_line()+
    geom_ribbon(aes(ymax = fit_high, ymin = fit_low,
                    fill = biomass_type), alpha = 0.4, color = NA)+
    geom_point(data = meas_mass2, aes(date, meas, col = biomass_type))+
    facet_grid(site~year, scales = 'free_x') +
    scale_color_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 600))+
    xlab('Date') +
    ylab(expression('Algal Biomass (AFDM g '~ m^-2*')')) +
    theme_classic()+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, 'line'))
cc <- qq %>%
    mutate(biomass_type = case_when(biomass_type == 'epil'~ 'Epilithon',
                                    biomass_type == 'fila' ~ 'Filamentous'),
           fit = case_when(fit < delta ~ delta,
                           TRUE ~ fit),
           fit_high = fit + se,
           fit_low = fit - se,
           fit_low = case_when(fit_low < delta ~ delta,
                               TRUE ~ fit_low)) %>%
    filter(units == 'chla_mgm2') %>%
    ggplot(aes(date, fit, col = biomass_type)) +
    geom_line()+
    geom_ribbon(aes(ymax = fit_high, ymin = fit_low,
                    fill = biomass_type), alpha = 0.4, color = NA)+
    geom_point(data = meas_chl2, aes(date, meas, col = biomass_type))+
    facet_grid(site~year, scales = 'free_x') +
    scale_color_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 1000))+
    xlab('Date') +
    ylab(expression('Algal Chlorophyll (mg chl a '~ m^-2*')')) +
    theme_classic()+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, 'line'))

png('figures/biomass_log_gamma_gams_comb_zeros.png', width = 7.5, height = 5, units = 'in',
    res = 300)

    ggpubr::ggarrange(mm, cc, nrow = 1, common.legend = TRUE,
                      labels = c('(a)', '(b)'))
dev.off()

# trim estimates so that they aren't more than 2 weeks from an actual measurement
qq %>%
  mutate(fit = case_when(se > 2*fit & se >100 ~ NA_real_,
                         fit < 0 ~ NA_real_,
                         fit > 1000 ~ NA_real_,
                         TRUE ~ fit),
         se = case_when(is.na(fit) ~ NA_real_,
                        TRUE ~ se))  %>%
    ggplot(aes(date, fit, col = units)) +
    geom_point()+
    facet_grid(site~biomass_type)


ndays = 0
mega_dates <- data.frame()
for(bt in c('epil', 'fila')){
    u_dates <- data.frame()
    for(u in c('gm2', 'chla_mgm2')){
        dates  <- meas %>%
            filter(biomass_type == bt, units == u)%>%
            select(date, site)

        cov_dates <- data.frame()
        for(s in unique(meas$site)){
            d <- filter(dates, site == s)
            dd <- unique(d$date)
            dlist <- vector()
            for(i in 1:length(dd)){

                dlist <- lubridate::as_date(c(dlist, seq((dd[i]-ndays),
                                                         (dd[i] + ndays),
                                                         by = 'day')))
            }
            cd = data.frame(date = unique(dlist)) %>%
                mutate(site = s)
            cov_dates <- bind_rows(cov_dates, cd)
        }
        cov_dates<- mutate(cov_dates, units = u)
        u_dates <- bind_rows(u_dates, cov_dates)
    }
    u_dates <- mutate(u_dates, biomass_type = bt)
    mega_dates <- bind_rows(mega_dates, u_dates)
}

mega_dates$coverage = 'good'
st_end <- mega_dates %>%
    mutate(year = factor(lubridate::year(date))) %>%
    group_by(site, year, units, biomass_type) %>%
    summarize(start = min(date),
              end = max(date)) %>%
    ungroup() %>%print(n = 48)

# qq2 <-
    left_join(qq, st_end, by = c('site','year', 'units', 'biomass_type')) %>%
    mutate(coverage = case_when(date > start & date < end ~ 'good',
                                TRUE ~ 'bad'),
           fit = case_when(coverage == 'bad' ~ NA_real_,
                           se > 2*fit ~ NA_real_,
                           fit > 1000 ~ NA_real_,
                           fit < 0 ~ NA_real_,
                           TRUE ~ fit),
           se = case_when(is.na(fit)~ NA_real_,
                          TRUE ~ se)) %>%
    ggplot(aes(date, fit, col = units)) + geom_line() +
    geom_ribbon(aes(ymin = fit-se, ymax = fit+se, fill = units), alpha = 0.3)+
    facet_grid(site~biomass_type)#, scales = 'free')
write_csv(qq2, 'data/biomass_data/invgamma_gam_fits_biomass.csv')
