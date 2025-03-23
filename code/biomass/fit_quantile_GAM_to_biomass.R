# fit GAMS to biomass data
# 10/2022

library(tidyverse)
library(qgam)
library(MASS)

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site),
           doy = as.numeric(format(date, '%j')))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           year = as.factor(year)) %>%
    filter(!is.na(site))

glimpse(biomass)

# There is a lot of missing data for the epiphyton and macrophyte categories

# combine filamentous algae and epiphyton for this analysis:
biomass <- biomass %>%
    mutate(site_year = as.factor(paste(site, year, sep = '_')),
           fila.om.area.g.m2 = case_when(is.na(epip.om.area.g.m2) ~ fila.om.area.g.m2,
                                         TRUE ~ fila.om.area.g.m2 + epip.om.area.g.m2),
           fila.chla.mg.m2.ritchie = case_when(is.na(epip.chla.mg.m2.ritchie) ~ fila.chla.mg.m2.ritchie,
                                         TRUE ~ fila.chla.mg.m2.ritchie + epip.chla.mg.m2.ritchie))
dates <- c(seq(min(biomass$date), max(biomass$date[biomass$year == 2020]), by = 'day'),
           seq(min(biomass$date[biomass$year == 2021]), max(biomass$date), by = 'day'))

preds <- data.frame(date = dates,
                    doy = as.numeric(format(dates, '%j')),
                    year = lubridate::year(dates))%>%
    mutate(year = as.factor(year))
s_preds <- data.frame()
for(site in unique(biomass$site)){
    preds$site = site
    s_preds <- bind_rows(s_preds, preds)
}
s_preds <- mutate(s_preds,
                  site_year = as.factor(paste(site, year, sep = '_'))) %>%
    tibble()
k_check <- data.frame()


# try it for all sites
if( suppressWarnings(require(RhpcBLASctl)) ){ blas_set_num_threads(1) } # Optional

fg2_log <- qgam(log(epil.om.area.g.m2) ~
                s(doy + site_year,  k = 20, bs = 'fs'),
            data = biomass,
            qu = 0.5)
# gam.check(fg2)

AIC(fg2_log)
gam.check(fg2_log)
kk <- k.check(fg2_log) %>% data.frame() %>%
    mutate(model = 'biofilm AFDM',
           smooth_par = c('s(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(kk) <- NULL
k_check <- bind_rows(k_check, kk)
pp <- mgcv::predict.gam(fg2_log, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                epil_gm2_fit = c(pp$fit),
                epil_gm2_se = c(pp$se.fit))

# filamentous
fg2_fila_log <- qgam(log(fila.om.area.g.m2) ~
                         s(doy + site_year,  k = 20, bs = 'fs'),
                     data = biomass,
                     qu = 0.5)

AIC(fg2_fila_log)
gam.check(fg2_fila_log)
kk <- k.check(fg2_fila_log) %>% data.frame() %>%
    mutate(model = 'filamentous AFDM',
           smooth_par = c( 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(kk) <- NULL
k_check <- bind_rows(k_check, kk)

pp <- mgcv::predict.gam(fg2_fila_log, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                fila_gm2_fit = pp$fit,
                fila_gm2_se = pp$se.fit)


# chlorophyll
# try it for all sites
fg2_chla_log <- qgam(log(epil.chla.mg.m2.ritchie) ~
                          s(doy + site_year,  k = 20, bs = 'fs'),
                      data = biomass,
                      qu = 0.5)
AIC(fg2_chla_log)

gam.check(fg2_chla_log)
kk <- k.check(fg2_chla_log) %>% data.frame() %>%
    mutate(model = 'biofilm chla',
           smooth_par = c( 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(kk) <- NULL
k_check <- bind_rows(k_check, kk)

pp <- mgcv::predict.gam(fg2_chla_log, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                epil_chla_mgm2_fit = pp$fit,
                epil_chla_mgm2_se = pp$se.fit)

# filamentous
fg2_fila_chla_log <- qgam(log(fila.chla.mg.m2.ritchie) ~
                          s(doy + site_year,  k = 20, bs = 'fs'),
                      data = biomass,
                      qu = 0.5)

AIC(fg2_fila_chla_log)
gam.check(fg2_fila_chla_log)
kk <- k.check(fg2_fila_chla_log) %>% data.frame() %>%
    mutate(model = 'filamentous chla',
           smooth_par = c( 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(kk) <- NULL
k_check <- bind_rows(k_check, kk)

pp <- mgcv::predict.gam(fg2_fila_chla_log, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                fila_chla_mgm2_fit = pp$fit,
                fila_chla_mgm2_se = pp$se.fit)

# plot diagnostics
png('figures/biomass_logGAMs_diagnostics.png', width = 7.5, height = 7.5,
    units = 'in', res = 300)
    par(mfrow = c(4,4),
        mar = c(3,4,2,2),
        oma = c(0,3,1,0))
    gam.check(fg2_log)
    mtext(expression(paste('biofilm g',m^-2)), 2, line = 45)
    gam.check(fg2_fila_log)
    mtext(expression(paste('fila. g',m^-2)), 2, line = 45)
    gam.check(fg2_chla_log)
    mtext(expression(paste('biofilm chla mg',m^-2)), 2, line = 45)
    gam.check(fg2_fila_chla_log)
    mtext(expression(paste('fila. chla mg',m^-2)), 2, line = 45)
dev.off()
png('figures/biomass_GAMs_diagnostics.png', width = 7.5, height = 7.5,
    units = 'in', res = 300)
    par(mfrow = c(4,4),
        mar = c(3,4,2,2),
        oma = c(0,3,1,0))
    gam.check(fg2)
    mtext(expression(paste('biofilm g',m^-2)), 2, line = 45)
    gam.check(fg2_fila)
    mtext(expression(paste('fila. g',m^-2)), 2, line = 45)
    gam.check(fg2_chla)
    mtext(expression(paste('biofilm chla mg',m^-2)), 2, line = 45)
    gam.check(fg2_fila_chla)
    mtext(expression(paste('fila. chla mg',m^-2)), 2, line = 45)
dev.off()

# plot GAMS
ddd <- s_preds %>%
    mutate(across(where(is.array), c),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    dplyr::select(date, doy, year, site, starts_with(c('epil_chla', 'fila_chla'))) %>%
    full_join(dplyr::select(biomass, date, site, sample, epil.chla.mg.m2.ritchie,
                     fila.chla.mg.m2.ritchie), by = c('date', 'site')) %>%
    rename(epil_chla_mgm2_meas = epil.chla.mg.m2.ritchie,
           fila_chla_mgm2_meas = fila.chla.mg.m2.ritchie) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 names_to = c('biomass_type', 'stat'),
                 names_pattern = '([a-z]+)_chla_mgm2_([a-z]+)',
                 values_to = 'value') %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    mutate(fit.low = exp(fit - se),
           fit.high = exp(fit + se),
           fit = exp(fit))

dpoly <- bind_rows(ddd, ddd[nrow(ddd):1,])
dpoly$fit[1:nrow(ddd)] <- dpoly$fit.high[1:nrow(ddd)]
dpoly$fit[(nrow(ddd)+1):nrow(dpoly)] <- dpoly$fit.low[(nrow(ddd)+1):nrow(dpoly)]
# png('figures/biomass_chla_gams.png', width = 5, height = 6, units = 'in',
#     res = 300)
chl <- ggplot(ddd, aes(date, fit, col = biomass_type)) +
    geom_line() +
    geom_point(aes(y = meas), size = 1.2) +
    # geom_line(aes(y = fit + se), lty = 2)+
    # geom_line(aes(y = fit - se), lty = 2)+
    geom_polygon(data = dpoly, aes(x = date, y = fit, fill = biomass_type),
                 alpha = 0.3, linetype = 0)+
    scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    facet_grid(site~year, scales = 'free_x')+
    xlab('Date') +
    ylab(expression('Algal Standing Crop (mg chl a '~ m^-2*')')) +
    scale_y_log10()+
    theme_bw() +
    theme(legend.position = 'none')

# dev.off()

ddd <-
    s_preds %>%
    mutate(across(where(is.array), c),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    dplyr::select(date, doy, year, site, starts_with(c('epil_gm2', 'fila_gm2'))) %>%
    full_join(dplyr::select(biomass, date, site, sample, epil.om.area.g.m2,
                     fila.om.area.g.m2), by = c('date', 'site')) %>%
    rename(epil_gm2_meas = epil.om.area.g.m2,
           fila_gm2_meas = fila.om.area.g.m2) %>%
    mutate(fila_gm2_meas  = ifelse(fila_gm2_meas > 300, NA, fila_gm2_meas))%>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 names_to = c('biomass_type', 'stat'),
                 names_pattern = '([a-z]+)_gm2_([a-z]+)',
                 values_to = 'value') %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    group_by(site, biomass_type) %>%
    mutate(across(c('fit', 'se'), ~zoo::na.approx(., x = date, na.rm = F))) %>%
    mutate(fit.low = exp(fit - se),
           fit.high = exp(fit + se),
           fit = exp(fit))

dpoly <- bind_rows(ddd, ddd[nrow(ddd):1,])
dpoly$fit[1:nrow(ddd)] <- dpoly$fit.high[1:nrow(ddd)]
dpoly$fit[(nrow(ddd)+1):nrow(dpoly)] <- dpoly$fit.low[(nrow(ddd)+1):nrow(dpoly)]
# dpoly <- bind_rows(ddd, ddd[nrow(ddd):1,])
# dpoly$fit[1:nrow(ddd)] <- dpoly$fit[1:nrow(ddd)] + dpoly$se[1:nrow(ddd)]
# dpoly$fit[(nrow(ddd)+1):nrow(dpoly)] <- dpoly$fit[(nrow(ddd)+1):nrow(dpoly)] -
#     dpoly$se[(nrow(ddd)+1):nrow(dpoly)]
# png('figures/biomass_gm2_gams.png', width = 5, height = 6, units = 'in',
#     res = 300)
mass <- ggplot(ddd, aes(date, fit, col = biomass_type)) +
    geom_line() +
    geom_point(aes(y = meas), size = 1.2) +
    # geom_line(aes(y = fit + se), lty = 2)+
    # geom_line(aes(y = fit - se), lty = 2)+
    geom_polygon(data = dpoly, aes(x = date, y = fit, fill = biomass_type),
                 alpha = 0.3, linetype = 0)+
    scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    facet_grid(site~year, scales = 'free_x')+
    scale_y_log10()+
    xlab('Date') +
    ylab(expression('Algal Standing Crop (AFDM g '~ m^-2*')')) +
    theme_bw() +
    theme(legend.position = 'none')

# dev.off()

png('figures/biomass_log_gams_comb.png', width = 7.5, height = 5, units = 'in',
    res = 300)
ggpubr::ggarrange(mass, chl, nrow = 1, ncol = 2, common.legend = TRUE,
                  labels = c('a', 'b'))
dev.off()



s_preds <- select(s_preds, -site_year)
# attributes(s_preds) <- NULL
as_tibble(s_preds)
qq = as_tibble(lapply(s_preds, c))
write_csv(qq, 'data/biomass_data/log_gam_fits_biomass.csv')
write_csv(k_check, 'data/biomass_data/log_gam_smoothness_parameter_checks.csv')
