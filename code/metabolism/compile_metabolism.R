# Compile metabolism estimates from UCFR sites
# A Carter 8/2022

library(tidyverse)
library(lubridate)
library(streamMetabolizer)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')
source('code/metabolism/functions_examine_SM_output.R')

site_dat <- read_csv('data/site_data.csv') %>%
    filter(!is.na(sitecode))

compiled_metab <- data.frame()
# 2020 sites ####
# Perkins 2020####
fit <- readRDS('data/metab_fits/PL_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/PL_2020.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
# bad_Rhats <- get_bad_Rhats(fit)
bad_days <- c(as.Date(c('2020-08-24', '2020-08-25')))#,
              # bad_Rhats)
met <- extract_metab(fit, sitecode = 'PL', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# Deer Lodge 2020####
fit <- readRDS('data/metab_fits/DL_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/DL_2020.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
met <- extract_metab(fit, sitecode = 'DL')

plot_KxER(met)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# Gold Creek 2020####
fit <- readRDS('data/metab_fits/GC_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/GC_2020.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit)
bad_days <- unique(c(as.Date(c('2020-08-24', '2020-08-25')),
              bad_Rhats))
met <- extract_metab(fit, sitecode = 'GC', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# Garrison 2020####
fit <- readRDS('data/metab_fits/GR_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/GR_2020.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- as.Date(c('2020-08-24', '2020-08-25'))
met <- extract_metab(fit, sitecode = 'GR', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# BearMouth 2020####
fit <- readRDS('data/metab_fits/BM_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/BM_2020.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit)
bad_days <- unique(bad_Rhats)
met <- extract_metab(fit, sitecode = 'BM', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# Bonita 2020####
fit <- readRDS('data/metab_fits/BN_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/BN_2020.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit)
bad_days <- unique(c(as.Date(c('2020-08-24', '2020-08-25')),
                     bad_Rhats))
met <- extract_metab(fit, sitecode = 'BN', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)


# 2021 ####
# Perkins 2021####
fit <- readRDS('data/metab_fits/PL_2021_kn_oipi.rds')
dat <- read_csv('data/prepared_data/PL_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
bad_Rhats <- get_bad_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- unique(c(as.Date(c('2021-07-01', '2021-07-20',
                               '2021-08-08')),
                     bad_Rhats))
met <- extract_metab(fit, sitecode = 'PL', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)


# Deer Lodge 2021####
fit <- readRDS('data/metab_fits/DL_2021_kn_oipi.rds')
dat <- read_csv('data/prepared_data/DL_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
bad_Rhats <- get_bad_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- unique(c(as.Date(c('2021-06-26', '2021-08-08')),
                     bad_Rhats))
met <- extract_metab(fit, sitecode = 'DL', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# Gold Creek 2021####
fit <- readRDS('data/metab_fits/GC_2021_kn_oipi.rds')
dat <- read_csv('data/prepared_data/GC_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
bad_Rhats <- get_bad_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- unique(c(as.Date(c('2021-07-12', '2021-07-13',
                               '2021-07-14', '2021-07-15',
                               '2021-08-05', '2021-08-08')),
                     bad_Rhats))
met <- extract_metab(fit, sitecode = 'GC', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# Garrison 2021####
fit <- readRDS('data/metab_fits/GR_2021_kn_oipi.rds')
dat <- read_csv('data/prepared_data/GR_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
bad_Rhats <- get_bad_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- unique(c(as.Date(c('2021-07-01', '2021-07-19',
                               '2021-07-20', '2021-07-27',
                               '2021-08-04',
                               '2021-08-05','2021-08-08')),
                     bad_Rhats))
met <- extract_metab(fit, sitecode = 'GR', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# BearMouth 2021####
fit <- readRDS('data/metab_fits/BM_2021_kn_oipi.rds')
dat <- read_csv('data/prepared_data/BM_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
bad_Rhats <- get_bad_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- unique(c(as.Date(c('2021-08-08'))))#,
                     # bad_Rhats))
met <- extract_metab(fit, sitecode = 'BM', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

# Bonita 2021####
fit <- readRDS('data/metab_fits/BN_2021_kn_oipi.rds')
dat <- read_csv('data/prepared_data/BN_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
bad_Rhats <- get_bad_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_days <- unique(c(as.Date(c('2021-08-01', '2021-08-02', '2021-08-08'))))#,
                     # bad_Rhats))
met <- extract_metab(fit, sitecode = 'BN', bad_days)

plot_KxER(met, rm.bds = TRUE)
plot_KxQ(met, dat)
plot_KxQ_bins(fit)

get_fit(fit)$overall %>%
    select(ends_with('Rhat'))
compiled_metab <- bind_rows(compiled_metab, met)

write_csv(compiled_metab, 'data/metabolism_compiled_all_sites.csv')
