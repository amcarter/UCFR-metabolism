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

# all_bad_days <- data.frame()
compiled_metab <- data.frame()
# Perkins ####
fit <- readRDS('data/metabolism/metab_fits/PL_knorm_oipi_2000iter_bdr_005.rds')
dat <- read_csv('data/prepared_data/PL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.05)
bad_days <- unique(as.Date(bad_Rhats))
met <- extract_metab(fit, sitecode = 'PL', bad_days = bad_days)

plot(met$K600, met$ER)
# identify(met$K600, met$ER, labels = met$date)
a <- plot_KxER(met, rm.bds = FALSE)
b <- plot_KxQ(met, dat)
# plot_KxQ_bins(fit)

png('figures/model_fits/PL_K_diagnostics.png', width = 6, height = 4,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE)
dev.off()

get_fit(fit)$overall %>% select(ends_with('Rhat')) %>%glimpse()
get_fit(fit)$KQ_overall %>% glimpse()
compiled_metab <- bind_rows(compiled_metab, met)


# Deer Lodge ####
fit <- readRDS('data/metabolism/metab_fits/DL_knorm_oipi_2000iter_bdr_05.rds')
dat <- read_csv('data/prepared_data/DL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(bad_Rhats)
met <- extract_metab(fit, sitecode = 'DL')#, bad_days)

a <- plot_KxER(met)#, rm.bds = TRUE)
b <- plot_KxQ(met, dat)
# plot_KxQ_bins(fit)

png('figures/model_fits/DL_K_diagnostics.png', width = 6, height = 4,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE)
dev.off()

get_fit(fit)$overall %>% glimpse()
get_fit(fit)$KQ_overall %>% glimpse()

compiled_metab <- bind_rows(compiled_metab, met)

# Garrison 2020####
fit <- readRDS('data/metabolism/metab_fits/GR_knorm_oipi_2000iter_bdr_01.rds')
dat <- read_csv('data/prepared_data/GR2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(bad_Rhats)
met <- extract_metab(fit, sitecode = 'GR')#, bad_days)

a <- plot_KxER(met, rm.bds = TRUE)
b <- plot_KxQ(met, dat)
# plot_KxQ_bins(fit)

png('figures/model_fits/GR_K_diagnostics.png', width = 6, height = 4,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE)
dev.off()

get_fit(fit)$overall %>% glimpse()
get_fit(fit)$KQ_overall %>% glimpse()

compiled_metab <- bind_rows(compiled_metab, met)


# Gold Creek 2020####
fit <- readRDS('data/metabolism/metab_fits/GC_knorm_oipi_2000iter_bdr_005.rds')
dat <- read_csv('data/prepared_data/GC2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(c(as.Date(c('2020-08-24')), bad_Rhats))
met <- extract_metab(fit, sitecode = 'GC', bad_days)

plot(met$K600, met$ER)
a <- plot_KxER(met, rm.bds = TRUE)
b <- plot_KxQ(met, dat)
# plot_KxQ_bins(fit)

png('figures/model_fits/GC_K_diagnostics.png', width = 6, height = 4,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE)
dev.off()

get_fit(fit)$overall %>% glimpse()
get_fit(fit)$KQ_overall %>% glimpse()

compiled_metab <- bind_rows(compiled_metab, met)



# BearMouth 2020####
fit <- readRDS('data/metabolism/metab_fits/BM_knorm_oipi_2000iter_bdr_05.rds')
dat <- read_csv('data/prepared_data/BM2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(bad_Rhats)
met <- extract_metab(fit, sitecode = 'BM')#, bad_days)

a <- plot_KxER(met, rm.bds = TRUE)
b <- plot_KxQ(met, dat)
# plot_KxQ_bins(fit)

png('figures/model_fits/BM_K_diagnostics.png', width = 6, height = 4,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE)
dev.off()

get_fit(fit)$overall %>% glimpse()
get_fit(fit)$KQ_overall %>% glimpse()

compiled_metab <- bind_rows(compiled_metab, met)


# Bonita 2020####
fit <- readRDS('data/metabolism/metab_fits/BN_knorm_oipi_2000iter_bdr_005.rds')
dat <- read_csv('data/prepared_data/BN2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.05)
bad_days <- unique(bad_Rhats)
met <- extract_metab(fit, sitecode = 'BN')#, bad_days)

plot(met$K600, met$ER)
a <- plot_KxER(met, rm.bds = TRUE)
b <- plot_KxQ(met, dat)
# plot_KxQ_bins(fit)

png('figures/model_fits/BN_K_diagnostics.png', width = 6, height = 4,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE)
dev.off()

get_fit(fit)$overall %>% glimpse()
get_fit(fit)$KQ_overall %>% glimpse()

compiled_metab <- bind_rows(compiled_metab, met)
compiled_metab <- compiled_metab %>%
    mutate(year = factor(year(date)),
           doy = as.numeric(format(date, '%j')),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    left_join(select(site_dat, site, distance_dwnstrm_km))


write_csv(compiled_metab, 'data/metabolism/metabolism_compiled_all_sites_2000iter_bdr_kss005.csv')
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
