# Compare UCFR metabolism data to biomass data for 2020 and 2021
# A Carter 8/2022

library(tidyverse)
library(lubridate)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

bm <- read_csv('data/biomass_data/biomass_working_data.csv')
met <- read_csv('data/metab_fits/metabolism_estimates_2020-21.csv')

glimpse(bm)
substrate <- bm %>%
    group_by(site) %>%
    summarize(across(starts_with('pct'), mean))

bm <- bm %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site)) %>%
    select(date, site, n_samp,
           average.depth = average.depth_mean, average.depth_sd,
           chla_mgm2 = total.chla.mg.m2.ritchie_mean,
           chla_mgm2_sd = total.chla.mg.m2.ritchie_sd,
           biomass_gm2 = total.algal.biomass.g.m2_mean,
           biomass_gm2_sd = total.algal.biomass.g.m2_sd)

d <- full_join(met, bm, by = c('site', 'date')) %>%
    select(-msgs.fit, -warnings, -errors) %>%
    mutate(year = year(date),
           doy = as.numeric(format(date, '%j')))%>%
    relocate(site, year, date, doy) %>%
    filter(site != 'CR', !is.na(site))
png('figures/GPP_biomass_visual_comp.png', width = 600 )
d %>%
    # filter(biomass_gm2 < 500) %>%
ggplot(aes(doy, GPP, col = factor(year))) +
    geom_line() +
    geom_point(aes(y = biomass_gm2/5)) +
    # geom_point(aes(y = chla_mgm2/10)) +
    facet_wrap(.~site)
dev.off()

d %>%
    rename(Biomass = biomass_gm2)%>%
    dplyr::filter(site == 'GR') %>%
    ggplot(aes(date, Biomass)) +
    geom_line(col = 'black') +
    # geom_point(aes(y = biomass_gm2), col = 'aquamarine3', size = 2) +
    facet_wrap(.~year, scales = 'free_x')+
    theme_bw() +
    labs(ylab = 'Biomass')
ggplot(d, aes(biomass_gm2, GPP, col = site)) +
    geom_point() +
    facet_wrap(.~site, scales = 'free')
