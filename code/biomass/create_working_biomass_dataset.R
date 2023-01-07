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
    select(date, site, sample, average.depth,
           cobble, sand, gravel, pebble,
           ends_with('om.area.g.m2'),
           ends_with('.chla.mg.m2.ritchie'),
           matches('phicocyanin'),
           total.algal.biomass.g.m2) %>%
    mutate(year = year(date))

bm_sum <- bm %>%
    group_by(date, site) %>%
    summarize(n_samp = n(),
              across(-n_samp, .fns = list(mean = function(x) mean(x, na.rm = T),
                                           sd = function(x) sd(x, na.rm = T)))) %>%
    select(-starts_with('sample'),
           -cobble_sd, -sand_sd, -gravel_sd, -pebble_sd) %>%
    rename(pct_cobble = cobble_mean, pct_sand = sand_mean,
           pct_gravel = gravel_mean, pct_pebble = pebble_mean) %>%
    ungroup() %>%
    mutate(year = year(date))

pairwise_sum = function(x, y){
    z = mapply(function(xx, yy) sum(xx, yy, na.rm=T), x, y)
    z[is.na(x) & is.na(y)] = NA_real_
    return(z)
}

bm_sum$filamentous_gm2 = pairwise_sum(bm_sum$epip.om.area.g.m2_mean,
                                      bm_sum$fila.om.area.g.m2_mean)
bm_sum$filamentous_chla_mgm2 = pairwise_sum(bm_sum$epip.chla.mg.m2.ritchie_mean,
                                      bm_sum$fila.chla.mg.m2.ritchie_mean)
bm_sum$filamentous_gm2_sd = sqrt(pairwise_sum(bm_sum$epip.om.area.g.m2_sd^2,
                                      bm_sum$fila.om.area.g.m2_sd^2))
bm_sum$filamentous_chla_mgm2_sd = sqrt(pairwise_sum(bm_sum$epip.chla.mg.m2.ritchie_sd^2,
                                      bm_sum$fila.chla.mg.m2.ritchie_sd^2))

bm_sum <- bm_sum %>%
    filter(site != 'CR' & !is.na(site)) %>%
    mutate(site = case_when(site == 'BG' ~'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site),
           doy = as.numeric(format(date, '%j'))) %>%
    rename(biofilm_gm2 = epil.om.area.g.m2_mean,
           biofilm_gm2_sd = epil.om.area.g.m2_sd,
           biofilm_chla_mgm2 = epil.chla.mg.m2.ritchie_mean,
           biofilm_chla_mgm2_sd = epil.chla.mg.m2.ritchie_sd) %>%
    select(-ends_with(c('ritchie_mean', 'ritchie_sd', '.m2_sd', '.m2_mean')),
           -year_mean, -year_sd)

ggplot(bm_sum, aes(doy, (biofilm_chla_mgm2), col = as.factor(year))) +
    geom_point() +
    geom_point(aes(y = (filamentous_chla_mgm2)), pch = 11) +
    facet_wrap(.~site)

ggplot(bm_sum, aes(date, total.algal.biomass.g.m2_mean, color = site)) +
    geom_point() +
    geom_line() +
    facet_wrap(.~year, scales = 'free_x')

write_csv(bm, 'data/biomass_data/biomass_working_data.csv')
write_csv(bm_sum, 'data/biomass_data/biomass_working_data_summary.csv')



# examine anomalously high chla values for fbom:

bm %>%
    filter(!is.na(site))%>%
ggplot(aes(fbom.area.g.m2*1000, fbom.chla.mg.m2.ritchie,
               col = factor(year)))+
    geom_point() +
    # geom_abline(slope = 1, intercept = 0) +
    facet_wrap(.~site, scales = 'free')
bm %>%
    filter(!is.na(site))%>%
ggplot(aes(epil.om.area.g.m2*1000, epil.chla.mg.m2.ritchie,
               col = factor(year)))+
    geom_point() +
    geom_abline(slope = .1, intercept = 0) +
    facet_wrap(.~site, scales = 'free')
bm %>%
    filter(!is.na(site))%>%
ggplot(aes(epip.om.area.g.m2*1000, epip.chla.mg.m2.ritchie,
               col = factor(year)))+
    geom_point() +
    geom_abline(slope = .1, intercept = 0) +
    facet_wrap(.~site, scales = 'free')
bm %>%
    filter(!is.na(site))%>%
ggplot(aes(fila.om.area.g.m2*1000, fila.chla.mg.m2.ritchie,
               col = factor(year)))+
    geom_point() +
    geom_abline(slope = .1, intercept = 0) +
    facet_wrap(.~site, scales = 'free')
