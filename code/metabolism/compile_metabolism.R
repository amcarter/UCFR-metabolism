# Compile metabolism estimates from UCFR sites
# A Carter 8/2022

library(tidyverse)
library(lubridate)
library(streamMetabolizer)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

site_dat <- read_csv('data/site_data.csv')
bad_days <- read_csv('data/metab_fits/bad_fits_2020.csv')
bad_days <- read_csv('data/metab_fits/bad_fits_2021.csv') %>%
    bind_rows(bad_days)

# load in metabolism fits and compile into one dataframe
metab20 <- data.frame()
metab21 <- data.frame()
# for(year in c(2020, 2021)){
year = 2020
    for(s in site_dat$sitecode){
        fit <- readRDS(paste0('data/metab_fits/', s, '_', year, '_kn_oipi.rds'))
        met <- predict_metab(fit)
        met$site <- s
        met$year <- year
        metab20 <- bind_rows(metab20, met)
    }
write_csv(metab20, 'data/metab_fits/metab20_compiled.csv')

year = 2021
    for(s in site_dat$sitecode){
        print(s)
        fit <- readRDS(paste0('data/metab_fits/', s, '_', year, '_kn_oipi.rds'))
        met <- predict_metab(fit)
        rm(fit)
        met$site <- s
        met$year <- year
        metab21 <- bind_rows(metab21, met)
        gc()
    }
write_csv(metab21, 'data/metab_fits/metab21_partial.csv')

metab21 <- read_csv('data/metab_fits/metab21_partial.csv')
metab20 <- read_csv('data/metab_fits/metab20_compiled.csv')

# }

met <- bind_rows(metab20, metab21) %>%
    left_join(bad_days, by = c('site', 'date')) %>%
    select(-msgs.fit, -warnings, -errors) %>%
    mutate(across(starts_with(c('GPP', 'ER')),
                  ~ case_when(fit == 'bad' ~ NA_real_,
                              TRUE ~ .)))

daily <- data.frame()
for(year in c(2020, 2021)){
    for(s in site_dat$sitecode){
        dat <- read_csv(paste0('data/prepared_data/', s, '_', year, '.csv'))
        dd <- dat %>%
            mutate(date = as.Date(solar.time)) %>%
            group_by(date) %>%
            summarize(across(-solar.time, mean, na.rm = T)) %>%
            ungroup()
        dd$site <- s
        dd$year <- year
        daily <- bind_rows(daily, dd)
    }
}

light <- read_csv('data/sw_radiation_all_sites') %>%
    rename(site = sitecode) %>%
    group_by(site, date) %>%
    summarize(SW = sum(SW)) %>%
    ungroup()

dat <- left_join(met, daily, by = c('site', 'year', 'date')) %>%
    left_join(light, by = c('site', 'date')) %>%
    relocate(site, year, date)

write_csv(dat, 'data/metab_fits/metabolism_estimates_2020-21.csv')

