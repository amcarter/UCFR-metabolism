# Calculate biomass numbers for results
library(tidyverse)

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    filter(!is.na(site) & site != 'CR') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           year = as.factor(lubridate::year(date)),
           site_year = paste(site, year, sep = '_'))

glimpse(biomass)


bb <- biomass %>% group_by(site, date, year) %>%
    summarize(across(ends_with('gm2'), ~mean(., na.rm = T)))

# Filamentous
min(bb$filamentous_gm2)
max(bb$filamentous_gm2)
mean(bb$filamentous_gm2)
sd(bb$filamentous_gm2)
mean(bb$filamentous_chla_mgm2)
sd(bb$filamentous_chla_mgm2, na.rm = T)
# Epilithon
min(bb$epilithon_gm2)
max(bb$epilithon_gm2)
mean(bb$epilithon_gm2)
sd(bb$epilithon_gm2)
mean(bb$epilithon_chla_mgm2)
sd(bb$epilithon_chla_mgm2)

bb %>% group_by(site, year) %>%
    summarize(across(ends_with('gm2'), ~mean(., na.rm = T))) %>%
    ungroup() %>%
    mutate(bloom = c(0,0,0,0,1,0,1,1,0,1,0,1)) %>%
    group_by(bloom) %>%
    summarize(across(ends_with('_gm2'), .fns = c(mean = ~mean(., na.rm = T),
                                                 sd = ~sd(., na.rm = T))))

bb %>%
    mutate(fila = filamentous_chla_mgm2/filamentous_gm2,
           epil = epilithon_chla_mgm2/epilithon_gm2) %>%
    ungroup() %>%
    summarize(across(c('fila', 'epil'), .fns = c(mean = ~mean(., na.rm = T),
                                                 sd = ~sd(., na.rm = T))))


bb %>% ungroup() %>%
    group_by(site, year) %>%
    summarize(across(ends_with('gm2'), ~mean(., na.rm = T))) %>%
    ungroup() %>%
    mutate(bloom = c(0,0,0,0,1,0,1,1,0,1,0,1),
           per_fila = filamentous_gm2/(filamentous_gm2 + epilithon_gm2)) %>%
    filter(bloom == 0) %>%
    select(per_fila) %>% summary()
