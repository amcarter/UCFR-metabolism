# Run streamMetabolizer on prepared datasets for UCFR

#Reinstall unitted and streamMetabolizer if needed
# remotes::install_github('appling/unitted', force = TRUE)
# remotes::install_github("GLEON/LakeMetabolizer", force = TRUE)
# remotes::install_github("USGS-R/streamMetabolizer", force = TRUE)

#load all packages
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(streamMetabolizer)
library(lubridate)
# library(shinystan)
library(tidyverse)
library(dygraphs)

# load datasets: ####
# setwd('~/Desktop/donkey/ucfr/')
# setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

PL<-read_csv('data/prepared_data/PL2020_2021.csv') %>% select(-discharge)
DL<-read_csv('data/prepared_data/DL2020_2021.csv') %>% select(-discharge)
GR<-read_csv('data/prepared_data/GR2020_2021.csv') %>% select(-discharge)
GC<-read_csv('data/prepared_data/GC2020_2021.csv') %>% select(-discharge)
BM<-read_csv('data/prepared_data/BM2020_2021.csv') %>% select(-discharge)
BN<-read_csv('data/prepared_data/BN2020_2021.csv') %>% select(-discharge)


#Run models

PL_daily <- PL %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = 13.6) %>%
    select(-discharge.daily)
PL_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=PL,
                       data_daily = PL_daily)
saveRDS(PL_fit, 'data/metabolism/metab_fits/PL_knorm_mle_bdr_0.rds')

DL_daily <- DL %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = 7.6) %>%
    select(-discharge.daily)
DL_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=DL,
                       data_daily = DL_daily)
saveRDS(DL_fit, 'data/metabolism/metab_fits/DL_knorm_mle_bdr_0.rds')

GR_daily <- GR %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = 18.1) %>%
    select(-discharge.daily)
GR_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=GR,
                       data_daily = GR_daily)
saveRDS(GR_fit, 'data/metabolism/metab_fits/GR_knorm_mle_bdr_0.rds')

GC_daily <- GC %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = 10.4) %>%
    select(-discharge.daily)
GC_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=GC,
                       data_daily = GC_daily)
saveRDS(GC_fit, 'data/metabolism/metab_fits/GC_knorm_mle_bdr_0.rds')

BM_daily <- BM %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = 9.7) %>%
    select(-discharge.daily)
BM_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=BM,
                       data_daily = BM_daily)
saveRDS(BM_fit, 'data/metabolism/metab_fits/BM_knorm_mle_bdr_0.rds')


BN_daily <- BN %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = 10.8) %>%
    select(-discharge.daily)
BN_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=BN,
                       data_daily = BN_daily)
saveRDS(BN_fit, 'data/metabolism/metab_fits/BN_knorm_mle_bdr_0.rds')
