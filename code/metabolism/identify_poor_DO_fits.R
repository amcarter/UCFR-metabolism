# Compile metabolism estimates from UCFR sites
# A Carter 8/2022
# setwd('~/Desktop/donkey/ucfr/')
# source('functions_examine_SM_output.R')
library(tidyverse)
library(lubridate)
library(streamMetabolizer)
source('code/metabolism/functions_examine_SM_output.R')

site_dat <- read_csv('data/site_data/site_data.csv') %>%
    mutate(site = factor(sitecode,
                         levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(sitecode)) %>%
    rename(distance_dwnstrm_km = 'Distance downstream (km)')

all_bad_days <- data.frame()
# Perkins ####
# fit <- readRDS('data/metabolism/metab_fits/initial runs/PL_knorm_oipi_1000iter.rds')
fit <- readRDS('data/metabolism/metab_fits/binned_ksig0.05/PL_kb_oipi_05.rds')
dat <- read_csv('data/prepared_data/PL2020_2021.csv')
# examine DO fit for bad days
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-07-01', '2021-07-07',
                      '2021-07-20', '2021-08-05', '2021-08-08'))
all_bad_days = data.frame(site = rep('PL', length(bad_days)),
                          bad_days = bad_days) %>%
    bind_rows(all_bad_days)

# Deer Lodge ####
# fit <- readRDS('data/metabolism/metab_fits/initial runs/DL_knorm_oipi_1000iter.rds')
fit <- readRDS('data/metabolism/metab_fits/binned_ksig0.05/DL_kb_oipi_05.rds')
dat <- read_csv('data/prepared_data/DL2020_2021.csv')
# examine DO fit for bad days
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2021-06-26', '2021-08-05', '2021-08-08'))
all_bad_days = data.frame(site = rep('DL', length(bad_days)),
                          bad_days = bad_days) %>%
    bind_rows(all_bad_days)

# Garrison 2020####
# fit <- readRDS('data/metabolism/metab_fits/initial runs/GR_knorm_oipi_1000iter.rds')
fit <- readRDS('data/metabolism/metab_fits/binned_ksig0.05/GR_kb_oipi_05.rds')
dat <- read_csv('data/prepared_data/GR2020_2021.csv')
# examine DO fit for bad days
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2020-09-24',
                      '2021-07-12', '2021-07-13', '2021-07-14',
                      '2021-07-15', '2021-07-20', '2021-08-05', '2021-08-08',
                      '2021-08-10'))
all_bad_days = data.frame(site = rep('GR', length(bad_days)),
                          bad_days = bad_days) %>%
    bind_rows(all_bad_days)

# Gold Creek 2020####
# fit <- readRDS('data/metabolism/metab_fits/initial runs/GC_knorm_oipi_1000iter.rds')
fit <- readRDS('data/metabolism/metab_fits/binned_ksig0.05/GC_kb_oipi_05.rds')
dat <- read_csv('data/prepared_data/GC2020_2021.csv')
# examine DO fit for bad days
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-07-01',
                      '2021-07-19', '2021-07-20', '2021-07-27',
                      '2021-08-04', '2021-08-05','2021-08-08'))
all_bad_days = data.frame(site = rep('GC', length(bad_days)),
                          bad_days = bad_days) %>%
    bind_rows(all_bad_days)

# BearMouth 2020####
# fit <- readRDS('data/metabolism/metab_fits/initial runs/BM_knorm_oipi_1000iter.rds')
fit <- readRDS('data/metabolism/metab_fits/binned_ksig0.05/BM_kb_oipi_05.rds')
dat <- read_csv('data/prepared_data/BM2020_2021.csv')
# examine DO fit for bad days
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-07-19', '2021-07-20',
                      '2021-07-21', '2021-08-08'))
all_bad_days = data.frame(site = rep('BM', length(bad_days)),
                          bad_days = bad_days) %>%
    bind_rows(all_bad_days)

# Bonita 2020####
fit <- readRDS('data/metabolism/metab_fits/binned_ksig0.05/BN_kb_oipi_05.rds')
# fit <- readRDS('data/metabolism/metab_fits/initial runs/BN_knorm_oipi_1000iter.rds')
dat <- read_csv('data/prepared_data/BN2020_2021.csv')
# examine DO fit for bad days
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-08-01', '2021-08-02',
                      '2021-08-08'))
all_bad_days = data.frame(site = rep('BN', length(bad_days)),
                          bad_days = bad_days) %>%
    bind_rows(all_bad_days)

write_csv(all_bad_days, 'data/metabolism/days_with_poor_DO_fits_initial_run.csv')
