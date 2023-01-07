# Compile metabolism estimates from UCFR sites
# A Carter 8/2022
# setwd('~/Desktop/donkey/ucfr/')
# source('functions_examine_SM_output.R')
library(tidyverse)
library(lubridate)
library(streamMetabolizer)
source('code/metabolism/functions_examine_SM_output.R')

site_dat <- read_csv('data/site_data.csv') %>%
    mutate(site = factor(sitecode,
                         levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(sitecode)) %>%
    rename(distance_dwnstrm_km = 'Distance downstream (km)')

# all_bad_days <- data.frame()
compiled_metab <- data.frame()
# Perkins ####
fit <- readRDS('data/metab_fits/PL_knorm_oipi.rds')
dat <- read_csv('data/prepared_data/PL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(c(as.Date(c('2020-10-16')), bad_Rhats))
# bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-07-01', '2021-07-07',
#                       '2021-07-20', '2021-08-02', '2021-08-05', '2021-08-08'))
# all_bad_days = data.frame(site = rep('PL', length(bad_days)),
#                           bad_days = bad_days) %>%
#     bind_rows(all_bad_days)
met <- extract_metab(fit, sitecode = 'PL', bad_days)

# plot(met$K600, met$ER)
# identify(met$K600, met$ER, labels = met$date)
a <- plot_KxER(met, rm.bds = TRUE)
b <- plot_KxQ(met, dat)
# plot_KxQ_bins(fit)

png('figures/model_fits/PL_K_diagnostics.png', width = 6, height = 4,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE)
dev.off()

get_fit(fit)$overall %>% glimpse()
get_fit(fit)$KQ_overall %>% glimpse()
compiled_metab <- bind_rows(compiled_metab, met)


# Deer Lodge ####
fit <- readRDS('data/metab_fits/DL_knorm_oipi.rds')
dat <- read_csv('data/prepared_data/DL2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(bad_Rhats)
# bad_days <- as.Date(c('2021-06-26', '2021-08-05', '2021-08-08'))
# all_bad_days = data.frame(site = rep('DL', length(bad_days)),
#                           bad_days = bad_days) %>%
#     bind_rows(all_bad_days)
met <- extract_metab(fit, sitecode = 'DL')#, bad_days)

a <- plot_KxER(met, rm.bds = TRUE)
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
fit <- readRDS('data/metab_fits/GR_knorm_oipi.rds')
dat <- read_csv('data/prepared_data/GR2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(bad_Rhats)
# bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2020-09-07', '2021-07-01',
#                       '2021-07-19', '2021-07-20', '2021-07-25', '2021-07-27',
#                       '2021-07-28', '2021-08-04', '2021-08-05','2021-08-08'))
# all_bad_days = data.frame(site = rep('GR', length(bad_days)),
#                           bad_days = bad_days) %>%
#     bind_rows(all_bad_days)
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
fit <- readRDS('data/metab_fits/GC_knorm_oipi.rds')
dat <- read_csv('data/prepared_data/GC2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(c(as.Date(c('2020-08-24', '2020-09-24')), bad_Rhats))
# bad_days <- as.Date(c('2020-07-22, 2020-08-24', '2020-08-25', '2020-09-22',
#                       '2021-06-30', '2021-07-12', '2021-07-13', '2021-07-14',
#                       '2021-07-15', '2021-07-20', '2021-08-05', '2021-08-08',
#                       '2021-08-10'))
# all_bad_days = data.frame(site = rep('GC', length(bad_days)),
#                           bad_days = bad_days) %>%
#     bind_rows(all_bad_days)
met <- extract_metab(fit, sitecode = 'GC', bad_days)

plot(met$K600, met$ER)
identify(x = met$K600, y = met$ER, labels = met$date)
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
fit <- readRDS('data/metab_fits/BM_knorm_oipi.rds')
dat <- read_csv('data/prepared_data/BM2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(bad_Rhats)
# bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-07-19', '2021-07-20',
#                       '2021-07-21', '2021-08-08'))
# all_bad_days = data.frame(site = rep('BM', length(bad_days)),
#                           bad_days = bad_days) %>%
#     bind_rows(all_bad_days)
met <- extract_metab(fit, sitecode = 'BM', bad_days)

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
fit <- readRDS('data/metab_fits/BN_knorm_oipi.rds')
dat <- read_csv('data/prepared_data/BN2020_2021.csv')
# examine DO fit for bad days
plot_metab_preds(fit)
plot_Rhats(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs')
bad_Rhats <- get_bad_Rhats(fit, threshold = 1.1)
bad_days <- unique(bad_Rhats)
# bad_days <- as.Date(c('2020-08-24', '2020-08-25', '2021-08-01', '2021-08-02',
#                       '2021-08-08'))
# all_bad_days = data.frame(site = rep('BN', length(bad_days)),
#                           bad_days = bad_days) %>%
#     bind_rows(all_bad_days)

met <- extract_metab(fit, sitecode = 'BN', bad_days)

plot(met$K600, met$ER)
identify(met$K600, met$ER, labels = met$date)
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


write_csv(compiled_metab, 'data/metabolism_compiled_all_sites.csv')
# write_csv(all_bad_days, 'data/days_with_poor_DO_fits.csv')
# png('figures/model_fits/QxK_across_sites.png', width = 8, height = 6.5,
png('figures/model_fits/Rhat_across_sites_normal_pool.png', width = 8, height = 6.5,
    res = 300, units = 'in')

    par(mfrow = c(2,3),
        mar = c(1,1,0,1),
        oma = c(4,4,3,0))
    fit <- readRDS('data/metab_fits/PL_knorm_oipi.rds')
    # plot_KxQ_bins(fit, labs = FALSE, legend = FALSE)
    plot_Rhats(fit)
    mtext('PL', 3, line = -2.5, adj = 0.9)
    rm(fit)
    fit <- readRDS('data/metab_fits/DL_knorm_oipi.rds')
    plot_Rhats(fit)
    # plot_KxQ_bins(fit, labs = FALSE, legend = FALSE)
    mtext('DL', 3, line = -2.5, adj = 0.9)
    rm(fit)
    fit <- readRDS('data/metab_fits/GR_knorm_oipi.rds')
    # plot_KxQ_bins(fit, labs = FALSE, legend = FALSE)
    plot_Rhats(fit)
    mtext('GR', 3, line = -2.5, adj = 0.9)
    rm(fit)
    fit <- readRDS('data/metab_fits/GC_knorm_oipi.rds')
    # plot_KxQ_bins(fit, labs = FALSE, legend = FALSE)
    plot_Rhats(fit)
    mtext('GC', 3, line = -2.5, adj = 0.9)
    rm(fit)
    fit <- readRDS('data/metab_fits/BM_knorm_oipi.rds')
    # plot_KxQ_bins(fit, labs = FALSE, legend = FALSE)
    plot_Rhats(fit)
    mtext('BM', 3, line = -2.5, adj = 0.9)
    rm(fit)
    fit <- readRDS('data/metab_fits/BN_knorm_oipi.rds')
    # plot_KxQ_bins(fit, labs = FALSE, legend = FALSE)
    plot_Rhats(fit)
    mtext('BN', 3, line = -2.5, adj = 0.9)

    par(mfrow = c(1,1), oma = c(0,0,2.5,0), new = T)
    # mtext(expression(paste("log discharge (m"^"3"~"s"^"-1"*")")), 1, adj = .55)
    # mtext(expression(paste("K600 (d"^"-1"*")")), 2, line = 0)
    mtext('Date', 1, adj = .55)
    mtext('Rhat', 2, line = 0)
    # legend('top', legend = c("prior", "posterior", "data"),
    #        col = c( "brown3", "brown3", "grey25"), xpd = NA, inset = c(0, -0.12),
    #        pch = c(1, 19, 20), bty = 'n', ncol = 3)
dev.off()
