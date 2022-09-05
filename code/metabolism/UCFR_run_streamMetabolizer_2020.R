# Run streamMetabolizer on prepared datasets for UCFR

#Reinstall unitted and streamMetabolizer if needed
# remotes::install_github('appling/unitted', force = TRUE)
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
setwd('~/Desktop/donkey/ucfr/')
# setwd('C:/Users/alice.carter/git/UCFR-metabolism/')
depth_Q <- read_csv('data/depth_discharge_relationships_allsites.csv')

PL<-read_csv('data/prepared_data/PL_2020.csv')
DL<-read_csv('data/prepared_data/DL_2020.csv')
GR<-read_csv('data/prepared_data/GR_2020.csv')
GC<-read_csv('data/prepared_data/GC_2020.csv')
BM<-read_csv('data/prepared_data/BM_2020.csv')
BN<-read_csv('data/prepared_data/BN_2020.csv')


# Visualize the data #####
# dat <- DL
#
# dat %>% unitted::v() %>%
#   mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
#   select(solar.time, starts_with('DO')) %>%
#   gather(type, DO.value, starts_with('DO')) %>%
#   mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
#   ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')
#
# labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)',
#             light='PAR\n(umol m^-2 s^-1)', discharge='Q\n(cms)')
# dat %>% unitted::v() %>%
#   mutate(discharge = log(discharge)) %>%
#   select(solar.time, depth, temp.water, light, discharge) %>%
#   gather(type, value, depth, temp.water, light, discharge) %>%
#   mutate(
#     type=ordered(type, levels=c('depth','temp.water','light','discharge')),
#     units=ordered(labels[type], unname(labels))) %>%
#   ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
#   facet_grid(units ~ ., scale='free_y') + theme_bw() +
#   scale_color_discrete('variable')

#sM model specs--------------- ####
bayes_name <- mm_name(type='bayes', pool_K600='binned',
                      err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_specs = specs(bayes_name,
                    burnin_steps = 2000,
                    saved_steps = 1000,
                    verbose = T, n_cores = 4)
# bayes_specs <- specs(bayes_name, #K600_daily_meanlog_meanlog=1.986474396, # where did this number for K600 come from?
#                   K600_daily_meanlog_sdlog=0.75,
#                   K600_daily_sdlog_sigma=0.5,
#                   burnin_steps=500,
#                   saved_steps=500,
#                   n_cores=4, verbose=T)

# Only relevant for running binned SM:
# set the nodes to reflect the range of discharge at this site:
set_Q_nodes <- function(bayes_specs, discharge){
    Qrange = quantile(log(discharge),
                      probs = c(0.1, 0.9), na.rm = T)
    n = 3
    delta = (Qrange[2]-Qrange[1])/n
    while(delta > 1){
        n = n + 1
        delta <- (Qrange[2]-Qrange[1])/n
    }
    nodes <- seq(Qrange[1], Qrange[2], length.out = n)
    bayes_specs$K600_lnQ_nodes_centers <- nodes
    bayes_specs$K600_lnQ_nodes_meanlog <- rep(2.484906649788, n)
    bayes_specs$K600_lnQ_nodes_sdlog=0.75
    bayes_specs$K600_daily_sigma_sigma = 0.5
    return(bayes_specs)
}

#Run models
bayes_specs <- set_Q_nodes(bayes_specs, PL$discharge)
PL_fit <- metab(bayes_specs, data=PL)
saveRDS(PL_fit, 'data/metab_fits/PL_2020_kn_oipi.rds')
rm(PL_fit)
gc()

bayes_specs <- set_Q_nodes(bayes_specs,  DL$discharge)
DL_fit <- metab(bayes_specs, data=DL)
saveRDS(DL_fit, 'data/metab_fits/DL_2020_kn_oipi.rds')
rm(DL_fit)
gc()

bayes_specs <- set_Q_nodes(bayes_specs, GR$discharge)
GR_fit <- metab(bayes_specs, data=GR)
saveRDS(GR_fit, 'data/metab_fits/GR_2020_kn_oipi.rds')
rm(GR_fit)
gc()

bayes_specs <- set_Q_nodes(bayes_specs, GC$discharge)
GC_fit <- metab(bayes_specs, data=GC)
saveRDS(GC_fit, 'data/metab_fits/GC_2020_kn_oipi.rds')
rm(GC_fit)
gc()

bayes_specs <- set_Q_nodes(bayes_specs, BM$discharge)
BM_fit <- metab(bayes_specs, data=BM)
saveRDS(BM_fit, 'data/metab_fits/BM_2020_kn_oipi.rds')
rm(BM_fit)
gc()

bayes_specs <- set_Q_nodes(bayes_specs, BN$discharge)
BN_fit <- metab(bayes_specs, data=BN)
saveRDS(BN_fit, 'data/metab_fits/BN_2020_kn_oipi.rds')
rm(BN_fit)
gc()
#Check model output: ####
bds <- data.frame()
# Perkins ####
fit <- readRDS('data/metab_fits/PL_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/PL_2020.csv')
daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(across(-solar.time, mean, na.rm = T)) %>%
    ungroup()

met <- fit@fit$daily
plot_metab_preds(fit)
plot(met$K600_daily_mean, met$ER_mean)

get_params(fit)
get_fit(fit)$overall %>%
  select(ends_with('Rhat'))
get_fit(fit)$overall %>%
  select('err_proc_iid_sigma_mean')

# predict_DO(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs',
              y_lim = list(conc = c(NA, NA),
                           pctsat = c(NA, NA), ddodt = c(-50, 50)))

# Dates with poor DO fits:
bad_days <- data.frame(site = 'PL',
                       date = as.Date(c('2020-08-24', '2020-08-25')),
                       fit = 'bad')
met %>%
    left_join(bad_days) %>%
    ggplot(aes(K600_daily_mean, ER_mean, col = fit)) +
    geom_point(size = 2)
daily %>%
    left_join(met, by = 'date') %>%
    left_join(bad_days, by = 'date') %>%
    ggplot(aes(log(discharge), K600_daily_mean, col = fit)) +
    geom_point(size = 2)

# these points don't look like they are influencing the overall KxER relationship
# Just remove these days but no need to rerun for now
bds <- bind_rows(bds, bad_days)



# Deer Lodge ####
fit <- readRDS('data/metab_fits/DL_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/DL_2020.csv')
daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(across(-solar.time, mean, na.rm = T)) %>%
    ungroup()

met <- fit@fit$daily
plot_metab_preds(fit)
plot(met$K600_daily_mean, met$ER_mean)

get_params(fit)
get_fit(fit)$overall %>%
  select(ends_with('Rhat'))
get_fit(fit)$overall %>%
  select('err_proc_iid_sigma_mean')

# predict_DO(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs',
              y_lim = list(conc = c(NA, NA),
                           pctsat = c(NA, NA), ddodt = c(-50, 50)))

# Dates with poor DO fits: everything looks good!

# Garrison ####
fit <- readRDS('data/metab_fits/GR_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/GR_2020.csv')
daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(across(-solar.time, mean, na.rm = T)) %>%
    ungroup()

met <- fit@fit$daily
plot_metab_preds(fit)
plot(met$K600_daily_mean, met$ER_mean)

get_params(fit)
get_fit(fit)$overall %>%
  select(ends_with('Rhat'))
get_fit(fit)$overall %>%
  select('err_proc_iid_sigma_mean')

# predict_DO(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs',
              y_lim = list(conc = c(NA, NA),
                           pctsat = c(NA, NA), ddodt = c(-50, 50)))

# Dates with poor DO fits:
bad_days <- data.frame(site = 'GR',
                       date = as.Date(c('2020-08-02', '2020-08-03',
                                        '2020-10-16')),
                       fit = 'bad')
met %>%
    left_join(bad_days) %>%
    ggplot(aes(K600_daily_mean, ER_mean, col = fit)) +
    geom_point(size = 2)
daily %>%
    left_join(met, by = 'date') %>%
    left_join(bad_days, by = 'date') %>%
    ggplot(aes(log(discharge), K600_daily_mean, col = fit)) +
    geom_point(size = 2)

bds <- bind_rows(bds, bad_days)
# Gold Creek ####
fit <- readRDS('data/metab_fits/GC_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/GC_2020.csv')
daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(across(-solar.time, mean, na.rm = T)) %>%
    ungroup()

met <- fit@fit$daily
plot_metab_preds(fit)
plot(met$K600_daily_mean, met$ER_mean)

get_params(fit)
get_fit(fit)$overall %>%
  select(ends_with('Rhat'))
get_fit(fit)$overall %>%
  select('err_proc_iid_sigma_mean')

# predict_DO(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs',
              y_lim = list(conc = c(NA, NA),
                           pctsat = c(NA, NA), ddodt = c(-50, 50)))

# Dates with poor DO fits:
bad_days <- data.frame(site = 'GR',
                       date = as.Date(c('2020-08-24', '2020-08-27')),
                       fit = 'bad')
met %>%
    left_join(bad_days) %>%
    ggplot(aes(K600_daily_mean, ER_mean, col = fit)) +
    geom_point(size = 2)
daily %>%
    left_join(met, by = 'date') %>%
    left_join(bad_days, by = 'date') %>%
    ggplot(aes(log(discharge), K600_daily_mean, col = fit)) +
    geom_point(size = 2)

bds <- bind_rows(bds, bad_days)



# Bear Gulch ####
fit <- readRDS('data/metab_fits/BM_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/BM_2020.csv')
daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(across(-solar.time, mean, na.rm = T)) %>%
    ungroup()

met <- fit@fit$daily
plot_metab_preds(fit)
plot(met$K600_daily_mean, met$ER_mean)

get_params(fit)
get_fit(fit)$overall %>%
  select(ends_with('Rhat'))
get_fit(fit)$overall %>%
  select('err_proc_iid_sigma_mean')

# predict_DO(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs',
              y_lim = list(conc = c(NA, NA),
                           pctsat = c(NA, NA), ddodt = c(-50, 50)))



# Bonita ####
fit <- readRDS('data/metab_fits/BN_2020_kn_oipi.rds')
dat <- read_csv('data/prepared_data/BN_2020.csv')
daily <- dat %>%
    mutate(date = as.Date(solar.time)) %>%
    group_by(date) %>%
    summarize(across(-solar.time, mean, na.rm = T)) %>%
    ungroup()

met <- fit@fit$daily
plot_metab_preds(fit)
plot(met$K600_daily_mean, met$ER_mean)
identify(met$K600_daily_mean, met$ER_mean, met$date)
get_params(fit)
get_fit(fit)$overall %>%
  select(ends_with('Rhat'))
get_fit(fit)$overall %>%
  select('err_proc_iid_sigma_mean')

# predict_DO(fit)
plot_DO_preds(fit, y_var=c( "pctsat"), style='dygraphs',
              y_lim = list(conc = c(NA, NA),
                           pctsat = c(NA, NA), ddodt = c(-50, 50)))



write_csv(bds, 'data/metab_fits/bad_fits_2020.csv')
