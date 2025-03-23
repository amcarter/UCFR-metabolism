# Calculate biomass numbers for results
library(tidyverse)
source('code/UCFR_helpers.R')


met <- read_csv('data/metabolism/metab_for_results.csv') %>%
    mutate(year = year(date),
           NEP = GPP + ER,
           PR = -GPP/ER,
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           site_year = paste(site, year, sep = '_'))%>%
    filter(!is.na(year))

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    filter(!is.na(site) & site != 'CR') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           year = as.factor(lubridate::year(date)),
           site_year = paste(site, year, sep = '_'))


head(biomass)


bb <- biomass %>% group_by(site, date, year) %>%
    summarize(across(ends_with('gm2'), ~mean(., na.rm = T)))

# Filamentous
min(bb$filamentous_gm2)
max(bb$filamentous_gm2)
mean(bb$filamentous_gm2)
calculate_ts_mean_se(bb$filamentous_gm2)
sd(bb$filamentous_gm2)
mean(bb$filamentous_chla_mgm2)
calculate_ts_mean_se(bb$filamentous_chla_mgm2)
sd(bb$filamentous_chla_mgm2, na.rm = T)
mean(bb$filamentous_chla_mgm2/bb$filamentous_gm2, na.rm = T)
calculate_ts_mean_se(bb$filamentous_chla_mgm2/bb$filamentous_gm2)
# Epilithon
min(bb$epilithon_gm2)
max(bb$epilithon_gm2)
mean(bb$epilithon_gm2)
calculate_ts_mean_se(bb$epilithon_gm2)
sd(bb$epilithon_gm2)
mean(bb$epilithon_chla_mgm2)
calculate_ts_mean_se(bb$epilithon_chla_mgm2)
sd(bb$epilithon_chla_mgm2)

mean(bb$epilithon_chla_mgm2/bb$epilithon_gm2)
calculate_ts_mean_se(bb$epilithon_chla_mgm2/bb$epilithon_gm2)


bio <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv')
bio <- left_join(select(met, 'site', 'date'), bio, by = c('site', 'date')) %>%
    filter(!is.na(year))

bio$site_year <- paste(bio$site, bio$year, sep = '_')
year_sum <- data.frame(site_year = unique(bio$site_year),
                       max_epil = NA_real_,
                       max_fila = NA_real_,
                       mean_epil = NA_real_,
                       mean_epil.se = NA_real_,
                       mean_fila = NA_real_,
                       mean_fila.se = NA_real_)

for(i in 1:nrow(year_sum)){
    mm <- filter(bio, site_year == year_sum$site_year[i])
    year_sum$max_epil[i] <- max(mm$epil_gm2_fit, na.rm = T)
    year_sum$max_fila[i] <- max(mm$fila_gm2_fit, na.rm = T)
    year_sum$mean_epil[i] <- mean(mm$epil_gm2_fit, na.rm = T)
    year_sum$mean_epil.se[i] <- calculate_ts_mean_se(mm$epil_gm2_fit)
    year_sum$mean_fila[i] <- mean(mm$fila_gm2_fit, na.rm = T)
    year_sum$mean_fila.se[i] <- calculate_ts_mean_se(mm$fila_gm2_fit)
}

summary(year_sum)
sd(year_sum$max_epil)/sqrt(12)
sd(year_sum$mean_epil)/sqrt(12)
sd(year_sum$max_fila)/sqrt(12)
sd(year_sum$mean_fila)/sqrt(12)

bloom <- filter(year_sum, site_year %in% c('GR_2020', 'GC_2020', 'GC_2021',
                                           'BG_2021', 'BN_2021'))
not_bloom <- filter(year_sum, !(site_year %in% c('GR_2020', 'GC_2020', 'GC_2021',
                                           'BG_2021', 'BN_2021')))
summary(bloom)
sd(bloom$max_fila)/sqrt(5)
sd(bloom$mean_fila)/sqrt(5)

summary(not_bloom)
sd(not_bloom$max_fila)/sqrt(7)
sd(not_bloom$mean_fila)/sqrt(7)


