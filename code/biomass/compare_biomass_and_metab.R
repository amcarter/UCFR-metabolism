# Compare UCFR metabolism data to biomass data for 2020 and 2021
# A Carter 8/2022

# setup ####
library(tidyverse)
library(lubridate)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))
met <- read_csv('data/metab_fits/metabolism_estimates_2020-21.csv')

glimpse(biomass)
substrate <- biomass %>%
    group_by(site) %>%
    summarize(across(starts_with('pct'), mean))


# Biomass categories ####
bm_class <- biomass %>%
    select(date, site, starts_with('phico'),
           matches('\\.om\\.area')) %>%
    rename_with(~sub('\\.om\\.area\\.g\\.m2', '', .))%>%
    filter(site != 'CR',
           !is.na(site))

bm_class[is.na(bm_class)] <- 0

bm_class <- bm_class %>%
    mutate(epilithion_mean = epil_mean + epip_mean,
           epilithion_sd = epil_sd + epip_sd,
           filamentous_mean = fila_mean + macro_mean,
           filamentous_sd = fila_sd + macro_sd)%>%
    select(site, date, starts_with(c('phico', 'epili', 'filam'))) %>%
    # rename_with(function(x) {str_match(x, '^([^\\.]+)')[, 2]}, ends_with('_mean')) %>%
    pivot_longer(any_of(ends_with(c('mean', 'sd'))),
                 names_to = c('biomass_type', 'stat'),
                 names_sep = '_',
                 values_to = 'mass_gm2') %>%
    pivot_wider(id_cols = c('date', 'site', 'biomass_type'),
                values_from = 'mass_gm2',
                names_from = 'stat')%>%
    rename(biomass_gm2 = mean)

png('figures/biomass_categories_by_site.png', width = 10, height = 6,
    units = 'in', res = 300)
    bm_class %>%
        mutate(year = year(date),
               doy = as.numeric(format(date, format = '%j'))) %>%
    ggplot( aes(doy, biomass_gm2, col = factor(year))) +
        geom_line() +
        facet_grid(biomass_type~site, scale = 'free') +
        theme_bw()
dev.off()

# Biomass and Metabolism ####
bm <- biomass %>%
    select(date, site, n_samp,
           average.depth = average.depth_mean, average.depth_sd,
           chla_mgm2 = total.chla.mg.m2.ritchie_mean,
           chla_mgm2_sd = total.chla.mg.m2.ritchie_sd,
           biomass_gm2 = total.algal.biomass.g.m2_mean,
           biomass_gm2_sd = total.algal.biomass.g.m2_sd,
           phicocyanin_mgm2 = phicocyanin.mg.m2_mean,
           phicocyanin_mgm2_sd = phicocyanin.mg.m2_sd)

d <- full_join(met, bm, by = c('site', 'date')) %>%
    mutate(year = year(date),
           doy = as.numeric(format(date, '%j')))%>%
    relocate(site, year, date, doy) %>%
    filter(site != 'CR', !is.na(site))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

png('figures/GPP_biomass_visual_comp.png', width = 800, height = 600 )
    d %>%
        # filter(chla_mgm2 < 500) %>%
    ggplot(aes(doy, temp.water, col = factor(year))) +
        geom_line() +
        # geom_point(aes(y = biomass_gm2/5), size = 2) +
        # geom_point(aes(y = chla_mgm2/10)) +
        facet_wrap(.~site, ncol = 2) +
        # ylim(0, 50)
        theme_bw()
dev.off()

ggplot(dat, aes(log(biomass_gm2), log(chla_mgm2))) +
    geom_point()
    ylim(-20,500)


png('figures/GPP_biomass_relationship.png', width = 800, height = 600 )
    ggplot(d, aes(biomass_gm2, GPP+ER, col = doy)) +
        geom_point(size = 3) +
        geom_smooth(method = 'lm', se = FALSE, lty = 2, col = 'grey')+
        scale_color_continuous(type = 'viridis')+
        facet_wrap(.~site, scales = 'free')+
        theme_bw()
dev.off()

    ggplot(d, aes(biomass_gm2, GPP, col = doy)) +
        geom_point(size = 3) +
        scale_color_continuous(type = 'viridis')+
        facet_wrap(.~site, scales = 'free')

write_csv(d, 'data/metabolism_biomass_compiled_all_sites.csv')
