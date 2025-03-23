# linear model(s) of GPP as a function of biomass

# setup ####

# library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(brms)

m <- read_csv('data/model_fits/biomass_metab_model_data.csv')

# Lm models

BM <- filter(m, site == 'BM' & !is.na(GPP) & !is.na(epil_chla_mgm2_fit))
BM$year <- year(BM$date)
mod <- lm(GPP ~ light + epil_chla_mgm2_fit + fila_chla_mgm2_fit,
          data = BM)

summary(mod)
BM$mod_resid <- residuals(mod)
ggplot(BM, aes(date, mod_resid)) + geom_point() + facet_wrap(.~year, scales = 'free')
acf(BM$GPP)
acf(BM$light)
acf(BM$fila_chla_mgm2_fit)
acf(BM$mod_resid)
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

saveRDS(brms_mods, 'data/model_fits/brms_gpp_models.rds')

bmods <- readRDS('data/model_fits/brms_gpp_models.rds')
plot(bmods$epil)
plot(conditional_effects(bmods$epil), points = TRUE)
pp_check(bmods$epil)
rr <- data.frame(residuals(bmods$epil))
plot(rr$Estimate)
bind_cols(m, rr)

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

