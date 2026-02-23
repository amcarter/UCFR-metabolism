# fit GAMS to biomass data
# 10/2022

library(tidyverse)
library(mgcv)

biomass <- read_csv('data/biomass_data/biomass_working_data_epip.csv') %>%
    filter(!is.na(site) & site != 'CR') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')),
           year = as.factor(lubridate::year(date)),
           site_year = paste(site, year, sep = '_'),
           epi_chla_mgm2 = epiphyton_chla_mgm2 + epilithon_chla_mgm2,
           epi_gm2 = epilithon_gm2 + epiphyton_gm2,
           filaepip_chla_mgm2 = epiphyton_chla_mgm2 + filamentous_chla_mgm2,
           filaepip_gm2 = epiphyton_gm2 + filamentous_gm2)

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


# s_preds_lin <- s_preds_gamma <-
#     mutate(s_preds,
#            site_year = as.factor(paste(site, year, sep = '_'))) %>%
#     tibble()

link_fn = 'log' # gamma link function, log or inverse?
delta = 0.5
par(mfrow = c(3,2))

hist(biomass$epilithon_gm2)
hist(log(biomass$epilithon_gm2 + min_vals$min[min_vals$biomass == 'epilithon_gm2'] ))
hist(biomass$epiphyton_gm2)
hist(log(biomass$epiphyton_gm2 + min_vals$min[min_vals$biomass == 'epiphyton_gm2'] ))
hist(biomass$filamentous_gm2)
hist(log(biomass$filamentous_gm2 + min_vals$min[min_vals$biomass == 'filamentous_gm2'] ))

par(mfrow = c(2,2))


# try it for all sites
biomass$site_year <- factor(biomass$site_year)
biomass$year <- factor(biomass$year)
biomass$site <- factor(biomass$site)

# epilithon
fg2_gamma <- gam(epilithon_gm2 +
                     min_vals$min[min_vals$biomass == 'epilithon_gm2'] ~ s(doy) +
                     s(doy, site_year, bs = 'fs'),
                 data = biomass, method = 'REML',
                 family = Gamma(link = link_fn))

AIC(fg2_gamma) # 2965
summary(fg2_gamma) #42.8%
gam.check(fg2_gamma)
sqrt(mean((fitted(fg2_gamma) - biomass$epilithon_gm2)^2)) # 8.68

pp_gamma <- mgcv::predict.gam(fg2_gamma, s_preds, type = 'response',
                              se.fit = TRUE)

# s_preds_lin <- mutate(s_preds_lin,
#                       epil_gm2_fit = c(pp_lin$fit),
#                       epil_gm2_se = c(pp_lin$se.fit))
s_preds <- mutate(s_preds,
                  epil_gm2_fit = c(pp_gamma$fit) -
                      min_vals$min[min_vals$biomass == 'epilithon_gm2'],
                  epil_gm2_se = c(pp_gamma$se.fit))

# epiphyton
# Presence/absence model
biomass$epip_present <- as.numeric(biomass$epiphyton_gm2 > 0)
fg2_epip_binom <- gam(epip_present ~ s(doy) + s(doy, site_year, bs = 'fs', k = 4),
                      data = biomass, family = binomial)

AIC(fg2_epip_binom) # 373
gam.check(fg2_epip_binom)
summary(fg2_epip_binom) #56.2%
sqrt(mean((fitted(fg2_epip_binom) - biomass$epip_present)^2)) #0.31

# Gamma model on positives only
fg2_epip_gamma <- gam(epiphyton_gm2 ~ s(doy, year, bs = 'fs') +
                          s(doy, site_year, bs = 'fs', k = 10),
                      data = subset(biomass, epiphyton_gm2 > 0),
                      family = Gamma(link = "log"), method = 'REML')
AIC(fg2_epip_gamma) #  both=786
sqrt(mean((fitted(fg2_epip_gamma) - biomass$epiphyton_gm2[biomass$epiphyton_gm2>0])^2)) # 7.7
gam.check(fg2_epip_gamma)
summary(fg2_epip_gamma) #78.3%


# pp_lin <- mgcv::predict.gam(fg2_epip, s_preds_lin, se.fit = TRUE)
pp_binom <- mgcv::predict.gam(fg2_epip_binom, newdata = s_preds,
                              type = "response", se.fit = TRUE)
pp_gamma1 <- mgcv::predict.gam(fg2_epip_gamma, type = 'response',
                               s_preds, se.fit = TRUE)
pp_gamma1$fit <- if_else(is.na(pp_gamma1$fit), 1, pp_gamma1$fit)
pp_gamma1$se.fit <- if_else(is.na(pp_gamma1$se.fit), 1, pp_gamma1$se.fit)
pp_gamma <- data.frame(fit = pp_binom$fit * pp_gamma1$fit,
                       se.fit = sqrt((pp_binom$fit^2 * pp_gamma1$se.fit^2) +
                                         (pp_gamma1$fit^2 * pp_binom$se.fit^2) +
                                         (pp_binom$se.fit^2 * pp_gamma1$se.fit^2)))


s_preds <- mutate(s_preds,
                  epip_gm2_fit = c(pp_gamma$fit),
                  epip_gm2_se = c(pp_gamma$se.fit))


# filamentous
# Presence/absence model
biomass$present <- as.numeric(biomass$filamentous_gm2 > 0)
fg2_fila_binom <- gam(present ~ s(doy) + s(doy, site_year, bs = 'fs', k = 4),
                      data = biomass, family = binomial)

AIC(fg2_fila_binom) # doy=634, s(doy) = 622, s(doy, siteyear)=315, both=307
gam.check(fg2_fila_binom)
summary(fg2_fila_binom) #65.2%
sqrt(mean((fitted(fg2_fila_binom) - biomass$present)^2)) #0.28
# Gamma model on positives only
fg2_fila_gamma <- gam(filamentous_gm2 ~ s(doy, year, bs = 'fs') +
                          s(doy, site_year, bs = 'fs', k = 10),
                      data = subset(biomass, filamentous_gm2 > 0),
                      family = Gamma(link = "log"), method = 'REML')
AIC(fg2_fila_gamma) # doy=2515; s(doy)=2488; s(doy, site_year)=2388; both=2390
sqrt(mean((fitted(fg2_fila_gamma) - biomass$filamentous_gm2[biomass$filamentous_gm2>0])^2)) # 55.7
gam.check(fg2_fila_gamma)
summary(fg2_fila_gamma) #55.3%

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

s_preds <- mutate(s_preds,
                  fila_gm2_fit = c(pp_gamma$fit),
                  fila_gm2_se = c(pp_gamma$se.fit))
##########

# chlorophyll
fg2_chla_gamma <- gam(epilithon_chla_mgm2 +
                          min_vals$min[min_vals$biomass == 'epilithon_chla_mgm2'] ~ s(doy) +
                          s(doy, site_year, bs = 'fs'),
                      data = biomass, method = 'REML',
                      family = Gamma(link = link_fn))

AIC(fg2_chla_gamma) # 3820
summary(fg2_chla_gamma) #49%
gam.check(fg2_chla_gamma)
sqrt(mean((fitted(fg2_chla_gamma) - biomass$epilithon_chla_mgm2)^2)) # 18.9

pp_chla_gamma <- mgcv::predict.gam(fg2_chla_gamma, s_preds, type = 'response',
                                   se.fit = TRUE)

s_preds <- mutate(s_preds,
                  epil_chla_mgm2_fit = c(pp_chla_gamma$fit) -
                      min_vals$min[min_vals$biomass == 'epilithon_chla_mgm2'],
                  epil_chla_mgm2_se = c(pp_chla_gamma$se.fit))

# epiphyton
# Presence/absence model
biomass$epip_chla_present <- as.numeric(biomass$epiphyton_chla_mgm2 > 0)
fg2_epip_chla_binom <- gam(epip_chla_present ~ s(doy) + s(doy, site_year, bs = 'fs', k = 4),
                           data = biomass, family = binomial)

AIC(fg2_epip_chla_binom) # 227
gam.check(fg2_epip_chla_binom)
summary(fg2_epip_chla_binom) #69.9%
sqrt(mean((fitted(fg2_epip_chla_binom) - biomass$epip_chla_present)^2)) #0.23

# Gamma model on positives only
fg2_epip_chla_gamma <- gam(epiphyton_chla_mgm2 ~ s(doy, year, bs = 'fs') +
                               s(doy, site_year, bs = 'fs', k = 10),
                           data = subset(biomass, epiphyton_chla_mgm2 > 0),
                           family = Gamma(link = "log"), method = 'REML')
AIC(fg2_epip_chla_gamma) # 691
sqrt(mean((fitted(fg2_epip_chla_gamma) - biomass$epiphyton_chla_mgm2[biomass$epiphyton_chla_mgm2>0])^2)) # 12.3
gam.check(fg2_epip_chla_gamma)
summary(fg2_epip_chla_gamma) #66.5%


pp_binom <- mgcv::predict.gam(fg2_epip_chla_binom, newdata = s_preds,
                              type = "response", se.fit = TRUE)
pp_gamma1 <- mgcv::predict.gam(fg2_epip_chla_gamma, type = 'response',
                               s_preds, se.fit = TRUE)
pp_gamma1$fit <- if_else(is.na(pp_gamma1$fit), 1, pp_gamma1$fit)
pp_gamma1$se.fit <- if_else(is.na(pp_gamma1$se.fit), 1, pp_gamma1$se.fit)
pp_gamma <- data.frame(fit = pp_binom$fit * pp_gamma1$fit,
                       se.fit = sqrt((pp_binom$fit^2 * pp_gamma1$se.fit^2) +
                                         (pp_gamma1$fit^2 * pp_binom$se.fit^2) +
                                         (pp_binom$se.fit^2 * pp_gamma1$se.fit^2)))


s_preds <- mutate(s_preds,
                  epip_chla_mgm2_fit = c(pp_gamma$fit),
                  epip_chla_mgm2_se = c(pp_gamma$se.fit))


# filamentous
# Presence/absence model
biomass$chla_present <- as.numeric(biomass$filamentous_chla_mgm2 > 0)
fg2_fila_chla_binom <- gam(chla_present ~ s(doy) + s(doy, site_year, bs = 'fs', k = 4),
                           data = biomass, family = binomial)

AIC(fg2_fila_chla_binom) # 313
gam.check(fg2_fila_chla_binom)
summary(fg2_fila_chla_binom) #65.2%
sqrt(mean((fitted(fg2_fila_chla_binom) - biomass$chla_present)^2)) #0.28
# Gamma model on positives only
fg2_fila_chla_gamma <- gam(filamentous_chla_mgm2 ~ s(doy, year, bs = 'fs') +
                               s(doy, site_year, bs = 'fs', k = 10),
                           data = subset(biomass, filamentous_chla_mgm2 > 0),
                           family = Gamma(link = "log"), method = 'REML')
AIC(fg2_fila_chla_gamma) # 2714
sqrt(mean((fitted(fg2_fila_chla_gamma) - biomass$filamentous_chla_mgm2[biomass$filamentous_chla_mgm2>0])^2)) # 82.7
gam.check(fg2_fila_chla_gamma)
summary(fg2_fila_chla_gamma) #49.9%

# pp_lin <- mgcv::predict.gam(fg2_fila, s_preds_lin, se.fit = TRUE)
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
################################################################################

ggplot(biomass, aes(date, epilithon_gm2))+
    geom_point() +
    geom_line(data = s_preds, aes(y = epil_gm2_fit))+
    facet_grid(site~year, scales = 'free_x')
ggplot(biomass, aes(date, epiphyton_gm2))+
    geom_point() +
    geom_line(data = s_preds, aes(y = epip_gm2_fit))+
    facet_grid(site~year, scales = 'free_x')
ggplot(biomass, aes(date, filamentous_gm2))+
    geom_point() +
    geom_line(data = s_preds, aes(y = fila_gm2_fit))+
    facet_grid(site~year, scales = 'free_x')


# Calculate fit metrics for the smoothness parameter

# gamma model
k1 <- k.check(fg2_gamma) %>% data.frame() %>%
    mutate(model = 'epilithon AFDM',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k1) <- NULL
k2 <- k.check(fg2_epip_gamma) %>% data.frame() %>%
    mutate(model = 'epiphyton AFDM',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k2) <- NULL
k3 <- k.check(fg2_fila_gamma) %>% data.frame() %>%
    mutate(model = 'filamentous AFDM',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k3) <- NULL
k4 <- k.check(fg2_chla_gamma) %>% data.frame() %>%
    mutate(model = 'epilithon chla',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k4) <- NULL
k5 <- k.check(fg2_epip_chla_gamma) %>% data.frame() %>%
    mutate(model = 'epiphyton chla',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k5) <- NULL
k6 <- k.check(fg2_fila_chla_gamma) %>% data.frame() %>%
    mutate(model = 'filamentous chla',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k6) <- NULL
k_check_gamma <- bind_rows(k1, k2, k3, k4, k5, k6)

# plot diagnostics
png('figures/SI/biomass_loggammaGAMs_diagnostics_epip.png', width = 7.5, height = 10,
    units = 'in', res = 300)
par(mfrow = c(6,4),
    mar = c(3,4,2,2),
    oma = c(0,3,2,0))
gam.check(fg2_gamma)
mtext(expression(paste('epil g',m^-2)), 2, line = 45)
gam.check(fg2_epip_gamma)
mtext(expression(paste('epip g',m^-2)), 2, line = 45)
gam.check(fg2_fila_gamma)
mtext(expression(paste('fila g',m^-2)), 2, line = 45)
gam.check(fg2_chla_gamma)
mtext(expression(paste('epil chla mg',m^-2)), 2, line = 45)
gam.check(fg2_epip_chla_gamma)
mtext(expression(paste('epip chla mg',m^-2)), 2, line = 45)
gam.check(fg2_fila_chla_gamma)
mtext(expression(paste('fila chla mg',m^-2)), 2, line = 45)
par(mfrow = c(1,1), new = T)
mtext('Model fit metrics for Biomass GAMS', line = 2.75)
dev.off()

qq = as_tibble(lapply(s_preds, c)) %>% select(-site_year) %>%
    mutate(across(starts_with(c('epi', 'fila')),
                  ~case_when(. < 0 ~ 0,
                             TRUE ~ .)))
write_csv(qq, 'data/biomass_data/log_gamma_gam_fits_biomass_epip.csv')
write_csv(k_check_gamma, 'data/biomass_data/log_gamma_gam_smoothness_parameter_checks_epip.csv')

qq <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass_epip.csv')
# plot GAMS
qq <- qq %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN'))) %>%
    pivot_longer(cols = starts_with(c('epi_','epil', 'epip', 'fila')),
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
               epip_gm2_meas = epiphyton_gm2,
               epi_gm2_meas = epi_gm2,
               fila_gm2_meas = filamentous_gm2,
               filaepip_gm2_meas = filaepip_gm2,
               epil_chla_mgm2_meas = epilithon_chla_mgm2,
               epip_chla_mgm2_meas = epiphyton_chla_mgm2,
               epi_chla_mgm2_meas = epi_chla_mgm2,
               filaepip_chla_mgm2_meas = filaepip_chla_mgm2,
               fila_chla_mgm2_meas = filamentous_chla_mgm2) %>%
    pivot_longer(cols = starts_with(c('epi', 'fila')),
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
                     meas = case_when(meas < delta ~ delta,
                                      TRUE ~ meas))
meas_chl2 <- mutate(meas_chl,
                    meas = case_when(meas < delta ~ delta,
                                     TRUE ~ meas))
mm <- qq %>% filter(units == 'gm2') %>%
    mutate(fit = case_when(fit < delta ~ delta,
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
    scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 600))+
    xlab('Date') +
    ylab(expression('Algal Standing Crop (AFDM g '~ m^-2*')')) +
    theme_classic()+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, 'line'))
cc <- qq %>%
    mutate(fit = case_when(fit < delta ~ delta,
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
    scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 1000))+
    xlab('Date') +
    ylab(expression('Algal Standing Crop (mg chl a '~ m^-2*')')) +
    theme_classic()+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, 'line'))

png('figures/biomass_log_gamma_gams_comb_zeros_epip.png', width = 7.5, height = 5, units = 'in',
    res = 300)

ggpubr::ggarrange(mm, cc, nrow = 1, common.legend = TRUE,
                  labels = c('a', 'b'))
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
