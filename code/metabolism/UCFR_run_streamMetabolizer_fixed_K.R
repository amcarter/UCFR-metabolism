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



# Compile metabolism estimates
source('code/metabolism/functions_examine_SM_output.R')

site_dat <- read_csv('data/site_data/site_data.csv') %>%
    mutate(site = factor(sitecode,
                         levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(sitecode)) %>%
    rename(distance_dwnstrm_km = 'Distance downstream (km)')

# all_bad_days <- data.frame()
compiled_metab <- data.frame()
# Perkins ####
fit <- readRDS('data/metabolism/metab_fits/PL_knorm_mle_bdr_0.rds')
dat <- read_csv('data/prepared_data/PL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
met <- extract_metab(fit, sitecode = 'PL', mle = TRUE)#, bad_days = bad_days)
compiled_metab <- bind_rows(compiled_metab, met)

# Deer Lodge ####
fit <- readRDS('data/metabolism/metab_fits/DL_knorm_mle_bdr_0.rds')
dat <- read_csv('data/prepared_data/DL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
met <- extract_metab(fit, sitecode = 'DL', mle = TRUE)#, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)

# Garrison 2020####
fit <- readRDS('data/metabolism/metab_fits/GR_knorm_mle_bdr_0.rds')
dat <- read_csv('data/prepared_data/GR2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
met <- extract_metab(fit, sitecode = 'GR', mle = TRUE)#, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)


# Gold Creek 2020####
fit <- readRDS('data/metabolism/metab_fits/GC_knorm_mle_bdr_0.rds')
dat <- read_csv('data/prepared_data/GC2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
met <- extract_metab(fit, sitecode = 'GC', mle = TRUE)
compiled_metab <- bind_rows(compiled_metab, met)

# BearMouth 2020####
fit <- readRDS('data/metabolism/metab_fits/BM_knorm_mle_bdr_0.rds')
dat <- read_csv('data/prepared_data/BM2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
met <- extract_metab(fit, sitecode = 'BM', mle = TRUE)#, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)


# Bonita 2020####
fit <- readRDS('data/metabolism/metab_fits/BN_knorm_mle_bdr_0.rds')
dat <- read_csv('data/prepared_data/BN2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
met <- extract_metab(fit, sitecode = 'BN', mle = TRUE)#, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)


compiled_metab <- compiled_metab %>%
    mutate(year = factor(year(date)),
           doy = as.numeric(format(date, '%j')),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    left_join(select(site_dat, site, distance_dwnstrm_km))


write_csv(compiled_metab, 'data/metabolism/metabolism_compiled_all_sites_mle_fixedK.csv')
# write_csv(all_bad_days, 'data/days_with_poor_DO_fits.csv')

met <- compiled_metab %>%
    select(-errors) %>%
    mutate(across(starts_with(c('GPP', 'ER', 'K600')),
                  ~case_when((!is.na(DO_fit) & DO_fit == 'bad') ~ NA_real_,
                             TRUE ~ .))) %>%
    select(-DO_fit)

dat <- read_csv('data/prepared_data/compiled_prepared_data.csv')
dd <- left_join(dat, met, by = c('site', 'date')) %>%
    select(-msgs.fit, -warnings, ends_with('Rhat') )

write_csv(dd, 'data/metabolism/metab_for_results.csv')

