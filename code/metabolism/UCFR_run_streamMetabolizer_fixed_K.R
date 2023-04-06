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


file_ext <- '_kfull_0.rds'

#Run models

K_meds <- read_csv('data/metabolism/median_K600s.csv')

PL_daily <- PL %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = K_meds$K600[K_meds$site == 'PL']) %>%
    select(-discharge.daily)
PL_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=PL,
                       data_daily = PL_daily)
saveRDS(PL_fit, paste0('data/metabolism/metab_fits/PL', file_ext))

DL_daily <- DL %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = K_meds$K600[K_meds$site == 'DL']) %>%
    select(-discharge.daily)
DL_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=DL,
                       data_daily = DL_daily)
saveRDS(DL_fit, paste0('data/metabolism/metab_fits/DL', file_ext))

GR_daily <- GR %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = K_meds$K600[K_meds$site == 'GR']) %>%
    select(-discharge.daily)
GR_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=GR,
                       data_daily = GR_daily)
saveRDS(GR_fit, paste0('data/metabolism/metab_fits/GR', file_ext))

GC_daily <- GC %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = K_meds$K600[K_meds$site == 'GC']) %>%
    select(-discharge.daily)
GC_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=GC,
                       data_daily = GC_daily)
saveRDS(GC_fit, paste0('data/metabolism/metab_fits/GC', file_ext))

BM_daily <- BM %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = K_meds$K600[K_meds$site == 'BM']) %>%
    select(-discharge.daily)
BM_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=BM,
                       data_daily = BM_daily)
saveRDS(BM_fit, paste0('data/metabolism/metab_fits/BM', file_ext))


BN_daily <- BN %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(depth, na.rm = T)) %>%
    ungroup() %>%
    mutate(K600.daily = K_meds$K600[K_meds$site == 'BN']) %>%
    select(-discharge.daily)
BN_fit <- metab_mle(specs = specs(mm_name(type = 'mle')), data=BN,
                       data_daily = BN_daily)
saveRDS(BN_fit, paste0('data/metabolism/metab_fits/BN', file_ext))



# Compile metabolism estimates
source('code/metabolism/functions_examine_SM_output.R')

site_dat <- read_csv('data/site_data/site_data.csv') %>%
    mutate(site = factor(sitecode,
                         levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(sitecode)) %>%
    rename(distance_dwnstrm_km = 'Distance downstream (km)')

compiled_metab <- data.frame()
# Perkins ####
fit <- readRDS(paste0('data/metabolism/metab_fits/PL', file_ext))
dat <- read_csv('data/prepared_data/PL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-10-23', '2021-07-01', '2021-07-21',
                      '2021-08-08', '2021-08-18'))
met <- extract_metab(fit, sitecode = 'PL', mle = TRUE, bad_days = bad_days)
compiled_metab <- bind_rows(compiled_metab, met)

# Deer Lodge ####
fit <- readRDS(paste0('data/metabolism/metab_fits/DL', file_ext))
dat <- read_csv('data/prepared_data/DL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-09-19', '2021-06-26', '2021-08-05', '2021-08-08',
                      '2021-09-14', '2021-09-21', '2021-09-22'))
met <- extract_metab(fit, sitecode = 'DL', mle = TRUE, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)

# Garrison####
fit <- readRDS(paste0('data/metabolism/metab_fits/GR', file_ext))
dat <- read_csv('data/prepared_data/GR2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-07-14', '2021-07-15',
                      '2021-08-08', '2021-08-18', '2021-09-22'))
met <- extract_metab(fit, sitecode = 'GR', mle = TRUE, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)


# Gold Creek ####
fit <- readRDS(paste0('data/metabolism/metab_fits/GC', file_ext))
dat <- read_csv('data/prepared_data/GC2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-27', '2021-08-08'))
met <- extract_metab(fit, sitecode = 'GC', mle = TRUE, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)

# BearMouth ####
fit <- readRDS(paste0('data/metabolism/metab_fits/BM', file_ext))
dat <- read_csv('data/prepared_data/BM2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-18', '2020-08-24', '2021-07-20', '2021-07-21'))
met <- extract_metab(fit, sitecode = 'BM', mle = TRUE, bad_days)
compiled_metab <- bind_rows(compiled_metab, met)


# Bonita 2020####
fit <- readRDS(paste0('data/metabolism/metab_fits/BN', file_ext))
dat <- read_csv('data/prepared_data/BN2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2021-08-01', '2021-08-02'))
met <- extract_metab(fit, sitecode = 'BN', mle = TRUE, bad_days)
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

