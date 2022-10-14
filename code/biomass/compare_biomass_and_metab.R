# Compare UCFR metabolism data to biomass data for 2020 and 2021
# A Carter 8/2022

# setup ####
library(tidyverse)
library(lubridate)
library(mgcv)

setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))
met <- read_csv('data/metabolism_compiled_all_sites.csv')
light <- read_csv('data/sw_radiation_all_sites') %>%
    group_by(sitecode, date) %>%
    summarize(light = sum(SW)) %>%
    ungroup() %>%
    rename(site = sitecode)
q <- read_csv('data/discharge_UCFRsites_2020.csv')

met <- left_join(met, light, by = c('site', 'date')) %>%
    left_join(q, by = c('site', 'date'))

glimpse(biomass)
substrate <- biomass %>%
    group_by(site) %>%
    summarize(across(starts_with('pct'), mean))


# Biomass categories ####
bm_class <- biomass %>%
    select(date, site, sample, starts_with('phico'),
           matches('\\.om\\.area'), matches('\\.chla\\.mg\\.m2\\.ritchie')) %>%
    rename_with(~sub('\\.om\\.area\\.g\\.m2', '_gm2', .))%>%
    rename_with(~sub('\\.chla\\.mg\\.m2\\.ritchie', '_chla_mgm2', .))%>%
    filter(site != 'CR',
           !is.na(site)) %>%
    mutate(year = year(date))

apply(bm_class, 2, function(x) sum(is.na(x)))

ggplot(bm_class, aes(date, epil_gm2)) +
    geom_point(col = 'forestgreen') +
    geom_point(aes(y = epip_gm2), col = 'aquamarine') +
    facet_grid(site~year, scales = 'free') +
    theme_bw()
bm_class[is.na(bm_class)] <- 0

bm <- bm_class %>%
    mutate(epilitheon_gm2 = epil_gm2 + epip_gm2,
           filamentous_gm2 = fila_gm2 + macro_gm2,
           epilitheon_chla.mgm2 = epil_chla_mgm2 + epip_chla_mgm2,
           filamentous_chla.mgm2 = fila_chla_mgm2 )%>%
    select(site, date, sample, starts_with(c( 'epili', 'filam', 'phico')))%>%
    mutate(year = year(date),
           doy = as.numeric(format(date, format = '%j')))

write_csv(bm, 'data/biomass_data/epil_fila_compiled_dataset.csv')

# Fit HGAMS to biomass data
fit_hgam <- function(variable , data, k = 6){

    data <- rename(data, x = !!variable)
    modGS <- gam(x ~ s(doy, sp = k, bs = 'fs') +
                     s(doy, site, sp = k, bs = "fs"),
                 data = data, method = "REML")

    newdata <- data.frame()
    for(s in unique(data$site)){
        nd <- data.frame(doy = seq(min(data$doy), max(data$doy)),
                         site = s)
        newdata <- bind_rows(newdata, nd)
    }
    pg <- predict.gam(modGS, se.fit = T, newdata = newdata)
    newdata[[paste0(variable, '_mean')]] <- unname(c(pg$fit))
    newdata[[paste0(variable, '_se')]] <- unname(c(pg$se.fit))

    return(newdata)
}

bm_2020 <- filter(bm, year == 2020)
bm_2021 <- filter(bm, year == 2021)

k = 8
gam_2020 <- fit_hgam('epilitheon_gm2', bm_2020, k)
gam_2020 <- fit_hgam('epilitheon_chla.mgm2', bm_2020, k) %>%
    full_join(gam_2020, by = c('doy', 'site'))
gam_2020 <- fit_hgam('filamentous_gm2', bm_2020, k)%>%
    full_join(gam_2020, by = c('doy', 'site'))
gam_2020 <- fit_hgam('filamentous_chla.mgm2', bm_2020, k)%>%
    full_join(gam_2020, by = c('doy', 'site'))
gam_2020 <- fit_hgam('phicocyanin.mg.m2', bm_2020, k)%>%
    full_join(gam_2020, by = c('doy', 'site'))
gam_2021 <- fit_hgam('epilitheon_gm2', bm_2021, k)
gam_2021 <- fit_hgam('epilitheon_chla.mgm2', bm_2021, k)%>%
    full_join(gam_2021, by = c('doy', 'site'))
gam_2021 <- fit_hgam('filamentous_gm2', bm_2021, k)%>%
    full_join(gam_2021, by = c('doy', 'site'))
gam_2021 <- fit_hgam('filamentous_chla.mgm2', bm_2021, k)%>%
    full_join(gam_2021, by = c('doy', 'site'))
gam_2021 <- fit_hgam('phicocyanin.mg.m2', bm_2021, k)%>%
    full_join(gam_2021, by = c('doy', 'site'))
gam_2020$year <- 2020
gam_2021$year <- 2021

bm_gams <- bind_rows(gam_2020,gam_2021)

write_csv(bm_gams, 'data/biomass_data/gam_fits_biomass.csv')
bm_gams <- read_csv('data/biomass_data/gam_fits_biomass.csv')

# Plot biomass data:
ggplot(bm_gams, aes(doy, epilitheon_gm2_mean, col = site)) +
    geom_line() +
    geom_point(data = bm, aes(doy, epilitheon_gm2, col = site))+
    facet_grid(site~year, scales = 'free')



bm_class <- bm %>%
    # rename(phicocyanin_mgm2 = phicocyanin.mg.m2)
    # rename_with(function(x) {str_match(x, '^([^\\.]+)')[, 2]}, ends_with('_mean')) %>%
    pivot_longer(any_of(ends_with('gm2')),
                 names_to = c('biomass_type', 'meas'),
                 names_sep = '_',
                 values_to = 'value') %>%
    pivot_wider(id_cols = c('date', 'doy', 'site', 'sample', 'biomass_type'),
                values_from = 'value',
                names_from = 'meas')%>%
    rename(biomass.gm2 = gm2) %>%
    mutate(year = substr(date, 1,4))

png('figures/biomass_categories_by_site.png', width = 10, height = 6,
    units = 'in', res = 300)
    bm_class %>%
    ggplot(aes(doy, biomass.gm2, col = factor(year))) +
        geom_line() +
        facet_grid(biomass_type~site, scale = 'free') +
        theme_bw()
dev.off()
png('figures/biomass_categories_by_site.png', width = 10, height = 6,
    units = 'in', res = 300)
    bm_class %>%
    ggplot( aes(doy, chla.mgm2, col = factor(year))) +
        geom_point() +
        geom_smooth(se = FALSE)+
        scale_y_log10()+
        facet_grid(biomass_type~site, scale = 'free') +
        theme_bw()
dev.off()

# Biomass and Metabolism ####
biomass <- read_csv('data/biomass_data/biomass_working_data_summary.csv') %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))
# bm <- biomass %>%
#     select(date, site, n_samp,
#            average.depth = average.depth_mean, average.depth_sd,
#            chla_mgm2 = total.chla.mg.m2.ritchie_mean,
#            chla_mgm2_sd = total.chla.mg.m2.ritchie_sd,
#            biomass_gm2 = total.algal.biomass.g.m2_mean,
#            biomass_gm2_sd = total.algal.biomass.g.m2_sd,
#            phicocyanin_mgm2 = phicocyanin.mg.m2_mean,
#            phicocyanin_mgm2_sd = phicocyanin.mg.m2_sd)
bm_nona <- mutate(biomass, across(where(is.numeric), ~ifelse(is.na(.), 0, .)))

bm <- bm_nona %>%
    mutate(epil_gm2 = epil.om.area.g.m2_mean +
               epip.om.area.g.m2_mean,
           epil_gm2_sd = sqrt(epil.om.area.g.m2_sd^2 +
               epip.om.area.g.m2_sd^2),
           epil_chla_mgm2 = epil.chla.mg.m2.ritchie_mean +
               epip.chla.mg.m2.ritchie_mean,
           epil_chla_mgm2_sd = sqrt(epil.chla.mg.m2.ritchie_sd^2 +
               epip.chla.mg.m2.ritchie_sd^2),
           fila_gm2 = fila.om.area.g.m2_mean +
               macro.om.area.g.m2_mean,
           fila_gm2_sd = sqrt(fila.om.area.g.m2_sd^2 +
               macro.om.area.g.m2_sd^2),
           fila_chla_mgm2 = fila.chla.mg.m2.ritchie_mean,
           fila_chla_mgm2_sd = fila.chla.mg.m2.ritchie_sd) %>%
    select(date, site,
           average.depth = average.depth_mean, average.depth_sd,
           ends_with(c('_gm2', 'gm2_sd', 'mgm2','mgm2_sd')))

bm %>%
    dplyr::rename(biofilm = epil_gm2, algae = fila_gm2) %>%
    pivot_longer(cols = any_of(c('biofilm', 'algae')),
                 names_to = 'Cat', values_to = 'AFDM_gm2') %>%
    mutate(year = year(date),
           doy = as.numeric(format(date, '%j'))) %>%
    filter(!is.na(site) )%>%
ggplot( aes(x = doy,y=AFDM_gm2, fill=Cat)) +
    geom_area() +
    facet_grid(year~site,  scales = 'free_x')+
    scale_fill_manual(values = c('#5B8853', '#A68A87')) +
    xlab('Day of Year') + ylab('Biomass (AFDM, gm-2)') +
    theme_bw()



d <- full_join(met, bm, by = c('site', 'date')) %>%
    mutate(year = year(date),
           doy = as.numeric(format(date, '%j')))%>%
    relocate(site, year, date, doy) %>%
    filter(site != 'CR', !is.na(site))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(site))

png('figures/GPP_biomass_visual_comp.png', width = 800, height = 600 )
    d %>%
        filter(!is.na(site)) %>%
    ggplot(aes(doy, GPP, col = factor(year))) +
        geom_line(col = '#007929') +
        geom_point(col = '#007929') +
        geom_line(aes(y = ER), col = '#A64B00') +
        geom_point(aes(y = ER), col = '#A64B00') +
        geom_hline(yintercept = 0, col = 'grey70')+
        # geom_point(aes(y = (epil_gm2 + fila_gm2)/5), size = 2) +
        # geom_point(aes(y = chla_mgm2/10)) +
        facet_grid(year~site) +
        # ylim(0,30)+
        xlab('Day of year')+
        ylab('GPP and biomass')+
        theme_bw()
dev.off()

gc <- d %>% filter(site == 'GC', year == 2020)
par(mfrow = c(2,1))
plot(gc$date, gc$GPP, col = '#6B8B3B', pch = 16, xlab = 'date', ylab = 'GPP')
plot(gc$date, gc$q.cms, col = '#2B6785', type = 'l', xlab = 'date',
     ylab = 'Discharge', log = 'y', lwd = 3)


# png('figures/GPP_biomass_relationship.png', width = 800, height = 600 )
    d %>%
        group_by(site)%>%
        # mutate(across(ends_with('gm2'),
        #               ~zoo::na.approx(., x = date, na.rm = F))) %>%
        ggplot( aes(epil_gm2 + fila_gm2, GPP, col = light)) +
        geom_point(size = 2) +
        geom_smooth(method = 'lm', se = FALSE, lty = 2, col = 'grey')+
        scale_color_continuous(type = 'viridis')+
        facet_wrap(.~site, scales = 'free')+
        theme_bw()
    d %>%
        group_by(site)%>%
        # mutate(across(ends_with('gm2'),
        #               ~zoo::na.approx(., x = date, na.rm = F))) %>%
        ggplot( aes(epil_chla_mgm2 + fila_chla_mgm2, GPP, col = light)) +
        geom_point(size = 2) +
        geom_smooth(method = 'lm', se = FALSE, lty = 2, col = 'grey')+
        scale_color_continuous(type = 'viridis')+
        facet_wrap(.~site, scales = 'free')+
        theme_bw()
    d %>%
        group_by(site)%>%
        mutate(across(ends_with('gm2'),
                      ~zoo::na.approx(., x = date, na.rm = F))) %>%
        ggplot( aes(light, GPP, col = log(epil_gm2 + fila_gm2))) +
        geom_point(size = 2) +
        geom_smooth(method = 'lm', se = FALSE, lty = 2, col = 'grey')+
        scale_color_continuous(type = 'viridis')+
        facet_wrap(.~site, scales = 'free')+
        labs(color = 'biomass')+
        theme_bw()
    d %>%
        filter(site == 'BN')%>%
        mutate(across(ends_with('gm2'),
                      ~zoo::na.approx(., x = date, na.rm = F))) %>%
        ggplot( aes(light, GPP, col = log(epil_gm2 + fila_gm2))) +
        geom_point(size = 2) +
        scale_color_continuous(type = 'viridis') +
        labs(color = 'biomass')+
        theme_bw()
# dev.off()

    ggplot(d, aes(biomass.gm2, GPP, col = doy)) +
        geom_point(size = 3) +
        scale_color_continuous(type = 'viridis')+
        facet_wrap(.~site, scales = 'free')

write_csv(d, 'data/metabolism_biomass_compiled_all_sites.csv')



# explore basic biomass models
par(mfrow =c(2,3))
d2 <- data.frame()
for(i in 1:6){
    s <- filter(d, site == unique(d$site)[i])
    acf(zoo::na.approx(s$GPP),  main = unique(d$site)[i])
    s$GPP_pre <- c(s$GPP[1], s$GPP[1:(nrow(s)-1)])
    # m <- lm(GPP ~ GPP_pre + light + fila_chla_mgm2 + epil_chla_mgm2, data = s)
    # summary(m)
    d2 <- bind_rows(d2, s)


}

for(i in 1:6){
    s <- filter(d, site == unique(d$site)[i],
                !is.na(fila_chla_mgm2))
    acf(zoo::na.approx(s$fila_chla_mgm2),  main = unique(d$site)[i])
}

for(i in 1:6){
    s <- filter(d, site == unique(d$site)[i],
                !is.na(epil_chla_mgm2))
    acf(zoo::na.approx(s$epil_chla_mgm2),  main = unique(d$site)[i])
}
m <- lme4::lmer(GPP ~ GPP_pre + light + fila_chla_mgm2 + epil_chla_mgm2 + (1|site),
         data = d2)

