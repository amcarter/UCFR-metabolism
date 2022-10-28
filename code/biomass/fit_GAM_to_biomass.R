# fit GAMS to biomass data
# 10/2022

library(tidyverse)
library(bayesGAM)


biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site),
           doy = as.numeric(format(date, '%j')))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           year = as.factor(year))

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

# test bayesGAM on the 2020 epilition data from deer lodge

dl_2020 <- biomass %>%
    mutate(year = as.factor(year))%>%
    filter(site == 'DL', year == 2020)

bg <- bayesGAM::bayesGAM(epil.om.area.g.m2 ~ doy, data = dl_2020)

preds <- data.frame(doy = seq(min(dl_2020$doy), max(dl_2020$doy), by = 1))

post <- bayesGAM::predict(bg, preds, draws = 200)

y <- data.frame(median = apply(post@pp$ypred, 2, FUN = median),
                y_2.5 = apply(post@pp$ypred, 2, FUN = quantile, probs = 0.025),
                y_97.5 = apply(post@pp$ypred, 2, FUN = quantile, probs = 0.975))

preds <- bind_cols(preds, y)
# try a non bayesian approach
library(mgcv)
fg <- gam(epil.om.area.g.m2 ~ s(doy),# by = year),
          data = dl_2020, method = "REML")

par(mfrow = c(2,1),
    mar = c(2,5,1,2), oma = c(3,0,1,0))
ylim = range(c(preds$y_2.5, preds$y_97.5))
plot(dl_2020$doy, dl_2020$epil.om.area.g.m2, ylim = ylim, ylab = 'BayesGAM fit')
lines(preds$doy, preds$median)
polygon(x = c(preds$doy, rev(preds$doy)),
        y = c(preds$y_2.5, rev(preds$y_97.5)),
        col = alpha('grey', .5), border = NA)

# par(mfrow = c(2,1))
plot(fg,
     residuals = TRUE,
     rug = FALSE,
     shade = TRUE,
     shift = coef(fg)[1],
     seWithMean = TRUE,
     ylab = 'mgcvGAM fit')

points(dl_2020$doy, dl_2020$epil.om.area.g.m2)

preds_2 <- mgcv::predict.gam(fg, preds, se.fit = TRUE)

preds$f_med <- preds_2$fit
preds$f_se <- preds_2$se.fit
ylim = range(c(preds$f_med + 1.96 * sqrt(nrow(dl_2020)) *preds$f_se,
               preds$f_med - 1.96* sqrt(nrow(dl_2020)) * preds$f_se))
plot(dl_2020$doy, dl_2020$epil.om.area.g.m2, ylim = ylim, ylab = 'mgcvGAM fit')

lines(preds$doy, preds$f_med)

summary(fg)
polygon(x = c(preds$doy, rev(preds$doy)),
        y = c(preds$f_med + 1.96 * sqrt(nrow(dl_2020)) *preds$f_se,
              rev(preds$f_med - 1.96 * sqrt(nrow(dl_2020)) * preds$f_se)),
        col = alpha('grey', .5), border = NA)
polygon(x = c(preds$doy, rev(preds$doy)),
        y = c(preds$f_med + preds$f_se,
              rev(preds$f_med -  preds$f_se)),
        col = alpha('grey', .5), border = NA)

# wy is this working better?
# fitting it with just time as a predictor is clearly not working very well
acf(preds$f_med)
biomass <- biomass %>%
    filter(!is.na(site)) %>%
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

# try it for all sites
fg2 <- gam(epil.om.area.g.m2 ~ s(doy) +
                  s(doy, site_year, bs = 'fs'),
              data = biomass, method = 'REML')

pp <- mgcv::predict.gam(fg2, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                epil_gm2_fit = c(pp$fit),
                epil_gm2_se = c(pp$se.fit))

# filamentous
fg2_fila <- gam(fila.om.area.g.m2 ~ s(doy) +
                  s(doy, site_year, bs = 'fs'),
              data = biomass, method = 'REML')

pp <- mgcv::predict.gam(fg2_fila, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                fila_gm2_fit = pp$fit,
                fila_gm2_se = pp$se.fit)


epi <- ggplot(s_preds, aes(doy, epil_gm2_fit, col = year)) +
    geom_line() +
    geom_line(aes(y = epil_gm2_fit + epil_gm2_se), lty = 2)+
    geom_line(aes(y = epil_gm2_fit - epil_gm2_se), lty = 2)+
    geom_point(aes(doy, epil.om.area.g.m2), data = biomass) +
    facet_wrap(.~site, scales = 'free_y', ncol = 1, strip.position = 'right')+
    ylab('Epilithion (g/m2)')

fila <- ggplot(s_preds, aes(doy, fila_gm2_fit, col = year)) +
    geom_line() +
    geom_line(aes(y = fila_gm2_fit + fila_gm2_se), lty = 2)+
    geom_line(aes(y = fila_gm2_fit - fila_gm2_se), lty = 2)+
    geom_point(aes(doy, fila.om.area.g.m2 ), data = biomass)+
    facet_wrap(.~site, scales = 'free_y', ncol = 1, strip.position = 'right')+
    ylab('Filamentous algae (g/m2)')
# chlorophyll
# try it for all sites
fg2_chla <- gam(epil.chla.mg.m2.ritchie ~ s(doy) +
                  s(doy, site_year, bs = 'fs'),
              data = biomass, method = 'REML')

pp <- mgcv::predict.gam(fg2_chla, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                epil_chla_mgm2_fit = pp$fit,
                epil_chla_mgm2_se = pp$se.fit)

# filamentous
fg2_fila_chla <- gam(fila.chla.mg.m2.ritchie ~ s(doy) +
                  s(doy, site_year, bs = 'fs'),
              data = biomass, method = 'REML')

pp <- mgcv::predict.gam(fg2_fila_chla, s_preds, se.fit = TRUE)

s_preds <- mutate(s_preds,
                fila_chla_mgm2_fit = pp$fit,
                fila_chla_mgm2_se = pp$se.fit)


epi_chla <- ggplot(s_preds, aes(doy, epil_chla_mgm2_fit, col = year)) +
    geom_line() +
    geom_line(aes(y = epil_chla_mgm2_fit + epil_chla_mgm2_se), lty = 2)+
    geom_line(aes(y = epil_chla_mgm2_fit - epil_chla_mgm2_se), lty = 2)+
    geom_point(aes(doy, epil.chla.mg.m2.ritchie), data = biomass) +
    facet_wrap(.~site, scales = 'free_y', ncol = 1, strip.position = 'right')+
    ylab('Epilithion Chla (mg/m2)')

fila_chla <- ggplot(s_preds, aes(doy, fila_chla_mgm2_fit, col = year)) +
    geom_line() +
    geom_line(aes(y = fila_chla_mgm2_fit + fila_chla_mgm2_se), lty = 2)+
    geom_line(aes(y = fila_chla_mgm2_fit - fila_chla_mgm2_se), lty = 2)+
    geom_point(aes(doy, fila.chla.mg.m2.ritchie ), data = biomass)+
    facet_wrap(.~site, scales = 'free_y', ncol = 1, strip.position = 'right')+
    ylab('Filamentous algae Chla (mg/m2)')


png('figures/biomass_chla_gams.png', height = 6, width = 8, res = 300, units = 'in')
    ggpubr::ggarrange( epi_chla, fila_chla,
                      nrow = 1, common.legend = TRUE)
dev.off()

png('figures/biomass_gams.png', height = 6, width = 8, res = 300, units = 'in')
    ggpubr::ggarrange(epi, fila,
                      nrow = 1, common.legend = TRUE)
dev.off()


s_preds <- select(s_preds, -site_year)
# attributes(s_preds) <- NULL
as_tibble(s_preds)
qq=as_tibble(lapply(s_preds, c))
write_csv(qq, 'data/biomass_data/gam_fits_biomass.csv')

