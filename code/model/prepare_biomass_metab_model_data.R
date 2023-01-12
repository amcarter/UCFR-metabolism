# compile data for running metabolism-biomass models
library(tidyverse)

met <- read_csv('data/metabolism/metabolism_compiled_all_sites.csv') %>%
    mutate(across(starts_with(c('GPP', 'ER', 'K600'), ignore.case = FALSE),
                  ~case_when(DO_fit == 'bad' ~ NA_real_,
                             TRUE ~ .)))
light <- read_csv('data/site_data/sw_radiation_all_sites.csv') %>%
    group_by(sitecode, date) %>%
    summarize(light = sum(SW)) %>%
    ungroup() %>%
    rename(site = sitecode) %>%
    mutate(light_rel = light/max(light))

q <- read_csv('data/site_data/discharge_UCFRsites_2020.csv')


met <- left_join(met, light, by = c('site', 'date')) %>%
    left_join(q, by = c('site', 'date'))

biomass <- read_csv('data/biomass_data/gam_fits_biomass.csv')

met <- left_join(met, biomass, by = c('site', 'date'))

m <- met %>%
    select(site, date, GPP, GPP.lower, GPP.upper,
           ER, ER.lower, ER.upper, K600, K600.lower, K600.upper,
           light = light_rel, q.cms,
           ends_with(c('2_fit', '2_se'))) %>%
    mutate(log_q = log(q.cms))

write_csv(m, 'data/model_fits/biomass_metab_model_data.csv')
