# linear model(s) of GPP as a function of biomass

# setup ####

# library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(brms)

met <- read_csv('data/metabolism_compiled_all_sites.csv')
light <- read_csv('data/sw_radiation_all_sites') %>%
    group_by(sitecode, date) %>%
    summarize(light = sum(SW)) %>%
    ungroup() %>%
    rename(site = sitecode)
q <- read_csv('data/discharge_UCFRsites_2020.csv')

met <- left_join(met, light, by = c('site', 'date')) %>%
    left_join(q, by = c('site', 'date'))

biomass <- read_csv('data/biomass_data/gam_fits_biomass.csv')

met <- left_join(met, biomass, by = c('site', 'date'))

m <- met %>%
    select(site, date, GPP, light, q.cms, ends_with(c('2_fit'))) %>%
    mutate(log_q = log(q.cms),
           across(-c(site, date), function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)))


# Fit BRMS models

epil <- brms::brm(GPP ~ epil_gm2_fit + light + (1|site), data = m)
fila <- brms::brm(GPP ~ fila_gm2_fit + light + (1|site), data = m)
epil_fila <- brms::brm(GPP ~ fila_gm2_fit + epil_gm2_fit +
                           light + (1|site), data = m)

epil_chl <- brms::brm(GPP ~ epil_chla_mgm2_fit + light + (1|site), data = m)
fila_chl <- brms::brm(GPP ~ fila_chla_mgm2_fit + light + (1|site), data = m)
epil_fila_chl <- brms::brm(GPP ~ fila_chla_mgm2_fit +
                               epil_chla_mgm2_fit + light + (1|site), data = m)

brms_mods <- list(epil = epil,
                  fila = fila,
                  epil_fila = epil_fila,
                  epil_chl = epil_chl,
                  fila_chl = fila_chl,
                  epil_fila_chl = epil_fila_chl)

saveRDS(brms_mods, 'data/brms_gpp_models.rds')

bmods <- readRDS('data/brms_gpp_models.rds')
plot(bmods$epil)
plot(conditional_effects(bmods$epil), points = TRUE)
pp_check(bmods$epil)

plot(bmods$fila)
plot(conditional_effects(bmods$fila), points = TRUE)
pp_check(bmods$fila)

plot(bmods$epil_fila)
plot(conditional_effects(bmods$epil_fila), points = TRUE)
pp_check(bmods$epil_fila)

plot(bmods$epil_chl)
plot(conditional_effects(bmods$epil_chl), points = TRUE)
pp_check(bmods$epil_chl)

plot(bmods$fila_chl)
plot(conditional_effects(bmods$fila_chl), points = TRUE)
pp_check(bmods$fila_chl)

plot(bmods$epil_fila_chl)
plot(conditional_effects(bmods$epil_fila_chl), points = TRUE)
pp_check(bmods$epil_fila_chl)

loo(bmods$epil_fila_chl, bmods$fila_chl)
