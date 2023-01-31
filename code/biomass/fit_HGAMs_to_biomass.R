# fit GAMS to biomass data
# 10/2022

library(tidyverse)
library(mgcv)

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site),
           doy = as.numeric(format(date, '%j')))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           year = as.factor(year)) %>%
    filter(!is.na(site))

glimpse(biomass)


# look at the data
biomass %>%
    select(date, site, ends_with('.om.area.g.m2')) %>%
    select(-cpom.om.area.g.m2)%>%
    filter(!is.na(site))%>%
    pivot_longer(cols = -c('date', 'site'), names_to = 'bm_category',
                 values_to = 'g_m2') %>%
    mutate(doy = format(date, "%j"),
           year = as.factor(substr(date, 1,4))) %>%
    ggplot(aes(doy, g_m2, col = year)) +
    geom_point() +
    facet_grid(site~bm_category, scales = 'free_y')


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
s_preds_lin <- s_preds_gamma <-
    mutate(s_preds,
           site_year = as.factor(paste(site, year, sep = '_'))) %>%
    tibble()

link_fn = 'log' # gamma link function, log or inverse?
par(mfrow = c(2,2))

# try it for all sites
fg2_gamma <- gam(epil.om.area.g.m2 ~ s(doy) +
                  s(doy, site_year, bs = 'fs'),
              data = biomass, method = 'REML',
              family = Gamma(link = link_fn))

fg2 <- gam(epil.om.area.g.m2 ~ s(doy) +
               s(doy, site_year, bs = 'fs'),
           data = biomass, method = 'REML', family = 'gaussian')

gam.check(fg2)
gam.check(fg2_gamma)
AIC(fg2)
AIC(fg2_gamma)

pp_lin <- mgcv::predict.gam(fg2, s_preds_lin, se.fit = TRUE)
pp_gamma <- mgcv::predict.gam(fg2_gamma, s_preds_gamma, type = 'response',
                        se.fit = TRUE)

s_preds_lin <- mutate(s_preds_lin,
                      epil_gm2_fit = c(pp_lin$fit),
                      epil_gm2_se = c(pp_lin$se.fit))
s_preds_gamma <- mutate(s_preds_gamma,
                        epil_gm2_fit = c(pp_gamma$fit),
                        epil_gm2_se = c(pp_gamma$se.fit))

# filamentous
fg2_fila_gamma <- gam(fila.om.area.g.m2 ~ s(doy) +
                          s(doy, site_year, bs = 'fs'),
                      data = biomass, method = 'REML',
                      family = Gamma(link = link_fn))

fg2_fila <- gam(fila.om.area.g.m2 ~ s(doy) +
                    s(doy, site_year, bs = 'fs'),
                data = biomass, method = 'REML', family = 'gaussian')

AIC(fg2_fila)
AIC(fg2_fila_gamma)

gam.check(fg2_fila)
gam.check(fg2_fila_gamma)


pp_lin <- mgcv::predict.gam(fg2_fila, s_preds_lin, se.fit = TRUE)
pp_gamma <- mgcv::predict.gam(fg2_fila_gamma, type = 'response',
                              s_preds_gamma, se.fit = TRUE)

s_preds_lin <- mutate(s_preds_lin,
                      fila_gm2_fit = c(pp_lin$fit),
                      fila_gm2_se = c(pp_lin$se.fit))
s_preds_gamma <- mutate(s_preds_gamma,
                        fila_gm2_fit = c(pp_gamma$fit),
                        fila_gm2_se = c(pp_gamma$se.fit))


# chlorophyll
fg2_chla_gamma <- gam(epil.chla.mg.m2.ritchie ~ s(doy) +
                          s(doy, site_year, bs = 'fs'),
                      data = biomass, method = 'REML',
                      family = Gamma(link = link_fn))
fg2_chla <- gam(epil.chla.mg.m2.ritchie ~ s(doy) +
                  s(doy, site_year, bs = 'fs'),
              data = biomass, method = 'REML', family = 'gaussian')

gam.check(fg2_chla)
gam.check(fg2_chla_gamma)

AIC(fg2_chla_gamma)
AIC(fg2_chla)

pp_gamma <- mgcv::predict.gam(fg2_chla_gamma, s_preds_gamma,
                              type ='response', se.fit = TRUE)
pp_lin <- mgcv::predict.gam(fg2_chla, s_preds_lin, se.fit = TRUE)

s_preds_lin <- mutate(s_preds_lin,
                      epil_chla_mgm2_fit = c(pp_lin$fit),
                      epil_chla_mgm2_se = c(pp_lin$se.fit))
s_preds_gamma <- mutate(s_preds_gamma,
                        epil_chla_mgm2_fit = c(pp_gamma$fit),
                        epil_chla_mgm2_se = c(pp_gamma$se.fit))

# filamentous
fg2_fila_chla_gamma <- gam(fila.chla.mg.m2.ritchie ~ s(doy) +
                               s(doy, site_year, bs = 'fs'),
                           data = biomass, method = 'REML',
                           family = Gamma(link = link_fn))
fg2_fila_chla <- gam(fila.chla.mg.m2.ritchie ~ s(doy) +
                        s(doy, site_year, bs = 'fs'),
                     data = biomass, method = 'REML', family = 'gaussian')

gam.check(fg2_fila_chla)
gam.check(fg2_fila_chla_gamma)
AIC(fg2_fila_chla)
AIC(fg2_fila_chla_gamma)


pp_gamma <- mgcv::predict.gam(fg2_fila_chla_gamma, s_preds_gamma,
                              type ='response', se.fit = TRUE)
pp_lin <- mgcv::predict.gam(fg2_fila_chla, s_preds_lin, se.fit = TRUE)

s_preds_lin <- mutate(s_preds_lin,
                      fila_chla_mgm2_fit = c(pp_lin$fit),
                      fila_chla_mgm2_se = c(pp_lin$se.fit))
s_preds_gamma <- mutate(s_preds_gamma,
                        fila_chla_mgm2_fit = c(pp_gamma$fit),
                        fila_chla_mgm2_se = c(pp_gamma$se.fit))

# Calculate fit metrics for the smoothness parameter
k1 <- k.check(fg2) %>% data.frame() %>%
    mutate(model = 'biofilm AFDM',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k1) <- NULL
k2 <- k.check(fg2_fila) %>% data.frame() %>%
    mutate(model = 'filamentous AFDM',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k2) <- NULL
k3 <- k.check(fg2_chla) %>% data.frame() %>%
    mutate(model = 'biofilm chla',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k3) <- NULL
k4 <- k.check(fg2_fila_chla) %>% data.frame() %>%
    mutate(model = 'filamentous chla',
           smooth_par = c('s(doy)', 's(doy, site_year)')) %>%
    relocate(model, smooth_par)
rownames(k4) <- NULL
k_check <- bind_rows(k1, k2, k3, k4)

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
png('figures/biomass_idGAMs_diagnostics.png', width = 7.5, height = 7.5,
    units = 'in', res = 300)
    par(mfrow = c(4,4),
        mar = c(3,4,2,2),
        oma = c(0,3,1,0))
    gam.check(fg2_gamma)
    mtext(expression(paste('biofilm g',m^-2)), 2, line = 45)
    gam.check(fg2_fila_gamma)
    mtext(expression(paste('fila. g',m^-2)), 2, line = 45)
    gam.check(fg2_chla_gamma)
    mtext(expression(paste('biofilm chla mg',m^-2)), 2, line = 45)
    gam.check(fg2_fila_chla_gamma)
    mtext(expression(paste('fila. chla mg',m^-2)), 2, line = 45)
dev.off()
png('figures/biomass_linGAMs_diagnostics.png', width = 7.5, height = 7.5,
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

qq = as_tibble(lapply(s_preds_lin, c)) %>% select(-site_year)
write_csv(qq, 'data/biomass_data/linear_gam_fits_biomass.csv')
write_csv(k_check, 'data/biomass_data/linear_gam_smoothness_parameter_checks.csv')

qq = as_tibble(lapply(s_preds_gamma, c)) %>% select(-site_year)
write_csv(qq, 'data/biomass_data/log_gamma_gam_fits_biomass.csv')
write_csv(k_check_gamma, 'data/biomass_data/log_gamma_gam_smoothness_parameter_checks.csv')

# plot GAMS
qq <- qq %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
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
       epil_gm2_meas = epil.om.area.g.m2,
       fila_gm2_meas = fila.om.area.g.m2,
       epil_chla_mgm2_meas = epil.chla.mg.m2.ritchie,
       fila_chla_mgm2_meas = fila.chla.mg.m2.ritchie) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 names_to = c('biomass_type', 'units', 'stat'),
                 names_pattern = '([a-z]+)_([a-z0-9_]+)_([a-z]+)',
                 values_to = 'value') %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    mutate(year = lubridate::year(date))
meas_chl <- filter(meas, units == 'chla_mgm2')
meas_mass <- filter(meas, units == 'gm2')
x = 10
mm <- qq %>% filter(units == 'gm2') %>%
    mutate(fit = fit + x,
           fit_high = fit + se,
           fit_low = fit - se,
           fit_low = case_when(fit_low < 0.3 ~ 0.3,
                               TRUE ~ fit_low)) %>%
    ggplot(aes(date, fit, col = biomass_type)) +
    geom_line()+
    geom_ribbon(aes(ymax = fit_high, ymin = fit_low,
                    fill = biomass_type), alpha = 0.4, color = NA)+
    geom_point(data = meas_mass, aes(date, meas, col = biomass_type))+
    facet_grid(site~year, scales = 'free_x') +
    scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 600))+
    xlab('Date') +
    ylab(expression('Algal Standing Crop (AFDM g '~ m^-2*')')) +
    theme_bw()
cc <- qq %>%
    mutate(fit = fit + x,
           fit_high = fit + se,
           fit_low = fit - se,
           fit_low = case_when(fit_low < 0.3 ~ 0.3,
                               TRUE ~ fit_low)) %>%
    filter(units == 'chla_mgm2') %>%
    ggplot(aes(date, fit, col = biomass_type)) +
    geom_line()+
    geom_ribbon(aes(ymax = fit_high, ymin = fit_low,
                    fill = biomass_type), alpha = 0.4, color = NA)+
    geom_point(data = meas_chl, aes(date, meas, col = biomass_type))+
    facet_grid(site~year, scales = 'free_x') +
    scale_color_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 1000))+
    xlab('Date') +
    ylab(expression('Algal Standing Crop (mg chl a '~ m^-2*')')) +
    theme_bw()

# png('figures/biomass_log_gamma_gams_comb.png', width = 7.5, height = 5, units = 'in',
#     res = 300)

    ggpubr::ggarrange(mm, cc, nrow = 1, common.legend = TRUE,
                      labels = c('a', 'b'))
# dev.off()

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
