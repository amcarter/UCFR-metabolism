# Compare UCFR metabolism data to biomass data for 2020 and 2021
# A Carter 8/2022

library(tidyverse)
library(lubridate)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

# Biomass data:
biomass <- read_csv('data/biomass_data/INTEGRATED_BIOMASS_FOR_AC_RFL.csv')
mdat <- read_csv('data/biomass_data/BIOMASS_DatasetAttributes.csv')
glimpse(biomass)

colnames(biomass) <- tolower(colnames(biomass))
unique(biomass$substrate)
bm <- biomass %>%
    mutate(date = as.Date(date, format = '%m/%d/%Y'),
           cobble = case_when(substrate == 'COBBLE' ~ 1,
                              TRUE ~ 0),
           sand = case_when(substrate == 'SAND' ~ 1,
                              TRUE ~ 0),
           gravel = case_when(substrate == 'GRAVEL' ~ 1,
                              TRUE ~ 0),
           pebble = case_when(substrate == 'PEBBLE' ~ 1,
                              TRUE ~ 0)) %>%
    select(date, site, average.depth,
           cobble, sand, gravel, pebble,
           ends_with('.om.area.g.m2'),
           ends_with('.chla.mg.m2.ritchie'),
           total.algal.biomass.g.m2) %>%
    group_by(date, site) %>%
    summarize(n_samp = n(),
              across(-n_samp, .fns = list(mean = function(x) mean(x, na.rm = T),
                                           sd = function(x) sd(x, na.rm = T)))) %>%
    select(-cobble_sd, -sand_sd, -gravel_sd, -pebble_sd) %>%
    rename(pct_cobble = cobble_mean, pct_sand = sand_mean,
           pct_gravel = gravel_mean, pct_pebble = pebble_mean) %>%
    ungroup() %>%
    mutate(year = year(date))

ggplot(bm, aes(date, total.algal.biomass.g.m2_mean, color = site)) +
    geom_point() +
    geom_line() +
    facet_wrap(.~year, scales = 'free_x')

write_csv(bm, 'data/biomass_data/biomass_working_data.csv')
