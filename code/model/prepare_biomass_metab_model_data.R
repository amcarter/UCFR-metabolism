# compile data for running metabolism-biomass models
library(tidyverse)
library(lubridate)

met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK.csv')# %>%
metse <- read_csv('data/metabolism/metabolism_compiled_all_sites_kb_oipi_05_bdr.csv') %>%
    mutate(GPP.plus = GPP.upper - GPP,
           GPP.minus = GPP - GPP.lower,
           ER.plus = ER.upper - ER,
           ER.minus = ER - ER.lower) %>%
    select(site, date, ends_with(c("plus", "minus")))

met <- left_join(met, metse, by = c("site", "date")) %>%
    mutate(GPP.upper = GPP + GPP.plus,
           GPP.lower = GPP - GPP.minus,
           ER.upper = ER + ER.plus,
           ER.lower = ER - ER.minus) %>%
    select(-ends_with(c("plus", "minus")))

write_csv(met, "data/metabolism/metabolism_compiled_all_sites_mle_fixedK_correctedSE.csv")
#     mutate(across(starts_with(c('GPP', 'ER', 'K600'), ignore.case = FALSE),
#                   ~case_when(DO_fit == 'bad' ~ NA_real_,
#                             TRUE ~ .)))
light <- read_csv('data/site_data/daily_modeled_light_all_sites.csv') %>%
    mutate(light_rel = PAR_surface/max(PAR_surface),
           mean_PAR_umolm2s = PAR_surface/light_hrs)

q <- read_csv('data/site_data/discharge_UCFRsites_2020.csv')


met <- left_join(met, light, by = c('site', 'date')) %>%
    left_join(q, by = c('site', 'date'))

biomass <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv') %>%
# biomass <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass_epip.csv') %>%
    mutate(across(starts_with(c('fila', 'epil')), ~case_when(. < 0 ~ 0,
                                                             TRUE ~ .)))
met <- left_join(met, biomass, by = c('site', 'date', 'doy', 'year'))

m <- met %>%
    select(site, date, GPP, GPP.lower, GPP.upper,
           ER, ER.lower, ER.upper, K600, #K600.lower, K600.upper,
           light = light_rel, mean_PAR_umolm2s, light_hrs, q.cms,
           ends_with(c('2_fit', '2_se'))) %>%
    mutate(log_q = log(q.cms),
           year = year(date))

write_csv(m, 'data/model_fits/biomass_metab_model_data.csv')
# write_csv(m, 'data/model_fits/biomass_metab_model_data_epip.csv')
