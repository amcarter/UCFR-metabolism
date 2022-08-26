# Compile metabolism estimates from UCFR sites
# A Carter 8/2022

library(tidyverse)
library(lubridate)
library(streamMetabolizer)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

site_dat <- read_csv('data/site_data.csv')

# load in metabolism fits and compile into one dataframe
metab <- data.frame()
for(year in c(2020, 2021)){
    for(s in site_dat$sitecode){
        fit <- readRDS(paste0('data/metab_fits/', s, '_', year, '_kn_oipi.rds'))
        met <- predict_metab(fit)
        met$site <- s
        met$year <- year
        metab <- bind_rows(metab, met)
    }
}


write_csv(metab, 'data/metab_fits/metabolism_estimates_2020-21.csv')

