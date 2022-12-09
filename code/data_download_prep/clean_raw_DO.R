# Clean oxygen timeseries for UCFR to flag bad or questionable data
# A Carter
# 7/2022

# NOTE#############################################################
# The raw sensor data for each site has both a UTC and MST datetime columns
# However, the MST column is actually in MT (ie, MST in winter, MDT in summer)
# This is an important distinction for how these times are converted to solartime
# For this reason, the UTC and unix time columns are the correct times.
####################################################################

# setup
library(tidyverse)
library(lubridate)
library(dygraphs)
library(xts)

setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

raw_files <- list.files('data/raw_DO/')
max_gap = 3 # maximum number of hours to fill in
# 2020 data ####
# Perkins Lane ####
ff <- grep('PERKINS_2020', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))

dat %>% select(water.temp, DO.obs, DO.sat) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dyRangeSelector()
flag_points <- data.frame(MST = c(ymd_hms(c(
    '2020-08-19 14:36:00', '2020-08-22 01:36:00',
    '2020-08-23 03:41:00', '2020-08-24 06:56:00',
    '2020-08-24 07:01:00', '2020-08-26 11:06:00',
    '2020-08-26 11:11:00',
    '2020-09-21 15:00:00', '2020-09-21 15:05:00',
    '2020-09-28 17:30:00', '2020-09-28 17:35:00',
    '2020-09-28 18:00:00', '2020-10-27 11:45:00',
    '2020-10-27 11:50:00')),
    seq(ymd_hms('2020-08-27 11:01:00'),
        ymd_hms('2020-08-27 11:29:00'),
        by = 'min'),
    seq(ymd_hms('2020-09-05 03:49:00'),
        ymd_hms('2020-09-05 09:29:00'),
        by = 'min'),
    seq(ymd_hms('2020-10-10 17:55:00'),
        ymd_hms('2020-10-10 18:05:00'),
        by = 'min')),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()

# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Deer Lodge ####
ff <- grep('DL_2020', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(ymd_hms(c(
    '2020-07-22 15:18:00', '2020-08-14 12:30:00',
    '2020-08-14 12:40:00', '2020-08-17 11:50:00',
    '2020-08-17 12:24:00')),
    seq(ymd_hms('2020-08-23 13:54:00'),
        ymd_hms('2020-08-26 11:41:00'),
        by = 'min')),
    flag = 'bad data')
dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:

dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '10 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Gold Creek ####
ff <- grep('GC_2020', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()
# looks good!

# align samples to consistent intervals:
dt <- as.numeric(median(diff(dat$UTC)))
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat) %>%
    mutate(across(c(-UTC, -MST),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC')%>%
    select(-DO.sat, -MST, -Q)

# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Garrison ####
ff <- grep('Gar_2020', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(
    ymd_hms(c('2020-08-14 14:42:00', '2020-08-10 13:52:00',
    '2020-09-02 12:23:00')),
    seq(ymd_hms('2020-07-16 10:40:00'),
        ymd_hms('2020-07-16 11:00:00'),
        by = 'min'),
    seq(ymd_hms('2020-07-22 13:51:00'),
        ymd_hms('2020-07-22 14:21:00'),
        by = 'min'),
    seq(ymd_hms('2020-08-02 16:51:00'),
        ymd_hms('2020-08-03 12:31:00'),
        by = 'min'),
    seq(ymd_hms('2020-09-01 09:27:00'),
        ymd_hms('2020-09-01 13:37:00'),
        by = 'min')),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '10 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Bear Gulch (Mouth?) ####
ff <- grep('BEARGULCH_2020', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(
    ymd_hms('2020-09-23 15:06:00'),
    seq(ymd_hms('2020-08-27 14:23:00'),
        ymd_hms('2020-08-27 14:36:00'),
        by = 'min')),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Bonita ####
ff <- grep('BONITA_2020', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(
    ymd_hms('2020-08-31 15:21:00'),
    seq(ymd_hms('2020-08-27 15:11:00'),
        ymd_hms('2020-08-27 15:21:00'),
        by = 'min'),
    seq(ymd_hms('2020-09-23 11:16:00'),
        ymd_hms('2020-09-23 14:28:00'),
        by = 'min')),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))


# 2021 data ####
# Perkins Lane ####
ff <- grep('PERKINS_2021', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))
DO_comp <- dat %>%
    select(UTC, DO_PL = DO.obs)

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(ymd_hms(c(
    '2021-06-23 18:03:00', '2021-06-25 04:48:00',
    '2021-06-25 05:33:00', '2021-06-25 06:48:00',
    '2021-07-08 03:53:00', '2021-07-17 06:08:00',
    '2021-08-26 17:58:00')),
    seq(ymd_hms('2021-06-14 15:42:00'),
        ymd_hms('2021-06-15 10:48:00'),
        by = 'min'),
    seq(ymd_hms('2021-06-23 22:03:00'),
        ymd_hms('2021-06-23 23:48:00'),
        by = 'min'),
    seq(ymd_hms('2021-06-24 23:48:00'),
        ymd_hms('2021-06-25 00:03:00'),
        by = 'min'),
    seq(ymd_hms('2021-07-26 21:08:00'),
        ymd_hms('2021-07-26 22:53:00'),
        by = 'min'),
    seq(ymd_hms('2021-11-15 10:40:00'),
        ymd_hms('2021-11-29 09:25:00'),
        by = 'min')
    ),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt,
                                x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()

# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Deer Lodge ####
ff <- grep('DL_2021', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))
DO_comp <- dat %>%
    select(UTC, DO_DL = DO.obs) %>%
    bind_rows(DO_comp)

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(
    seq(ymd_hms('2021-06-14 15:36:00'),
        ymd_hms('2021-06-15 11:34:00'),
        by = 'min'),
    seq(ymd_hms('2021-06-16 18:04:00'),
        ymd_hms('2021-06-18 20:34:00'),
        by = 'min'),
    seq(ymd_hms('2021-06-27 12:04:00'),
        ymd_hms('2021-06-30 13:53:00'),
        by = 'min'),
    seq(ymd_hms('2021-07-26 20:38:00'),
        ymd_hms('2021-08-02 05:53:00'),
        by = 'min'),
    seq(ymd_hms('2021-11-15 10:56:00'),
        ymd_hms('2021-11-29 09:11:00'),
        by = 'min')
    ),
    flag = 'bad data')
dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:

dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '10 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Gold Creek ####
ff <- grep('GC_2021', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))
DO_comp <- dat %>%
    select(UTC, DO_GC = DO.obs) %>%
    bind_rows(DO_comp)

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()
dat %>% select(DO.sat) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.sat', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(ymd_hms(
    c('2021-07-20 05:47:00', '2021-07-30 13:02:00', '2021-07-31 03:17:00',
      '2021-07-31 10:02:00', '2021-07-31 10:17:00', '2021-08-01 09:32:00',
      '2021-08-01 11:47:00', '2021-08-03 16:32:00', '2021-08-03 18:47:00',
      '2021-08-09 17:02:00', '2021-08-09 17:17:00', '2021-08-24 15:17:00',
      '2021-08-29 16:02:00', '2021-08-29 16:17:00', '2021-08-31 13:47:00',
      '2021-08-31 15:17:00', '2021-09-11 19:46:00', '2021-09-11 20:16:00',
      '2021-09-13 13:46:00', '2021-09-16 11:31:00', '2021-09-19 14:31:00',
      '2021-09-25 00:31:00', '2021-09-25 01:16:00', '2021-09-26 19:16:00',
      '2021-10-01 13:01:00', '2021-10-01 13:46:00', '2021-10-01 14:01:00',
      '2021-10-02 11:31:00', '2021-10-02 12:46:00', '2021-10-02 13:16:00',
      '2021-10-03 08:16:00', '2021-10-03 15:01:00', '2021-10-03 19:46:00',
      '2021-10-04 05:01:00', '2021-10-04 05:16:00', '2021-10-04 17:31:00',
      '2021-10-06 12:01:00', '2021-10-06 12:16:00', '2021-10-07 05:16:00',
      '2021-10-07 05:31:00', '2021-10-07 10:16:00', '2021-10-07 14:16:00',
      '2021-10-09 20:01:00', '2021-10-10 06:46:00', '2021-10-10 15:01:00',
      '2021-10-10 15:16:00', '2021-10-12 01:46:00', '2021-10-15 22:31:00',
      '2021-10-16 10:16:00', '2021-10-16 10:31:00', '2021-10-16 10:46:00',
      '2021-10-16 17:01:00', '2021-10-22 03:46:00', '2021-10-22 05:01:00',
      '2021-10-22 16:31:00', '2021-10-23 07:46:00', '2021-10-24 06:01:00',
      '2021-10-24 13:16:00', '2021-10-24 14:31:00', '2021-10-24 17:01:00',
      '2021-10-24 17:16:00', '2021-10-29 06:46:00')),
    seq(ymd_hms('2021-07-19 17:47:00'), ymd_hms('2021-07-19 18:47:00'), by = 'min'),
    seq(ymd_hms('2021-07-25 13:47:00'), ymd_hms('2021-07-25 14:32:00'), by = 'min'),
    seq(ymd_hms('2021-07-29 01:47:00'), ymd_hms('2021-07-29 03:47:00'), by = 'min'),
    seq(ymd_hms('2021-08-02 14:17:00'), ymd_hms('2021-08-02 16:47:00'), by = 'min'),
    seq(ymd_hms('2021-09-08 00:31:00'), ymd_hms('2021-09-08 01:31:00'), by = 'min'),
    seq(ymd_hms('2021-09-08 05:46:00'), ymd_hms('2021-09-08 07:16:00'), by = 'min'),
    seq(ymd_hms('2021-09-13 16:31:00'), ymd_hms('2021-09-13 17:16:00'), by = 'min'),
    seq(ymd_hms('2021-09-24 21:16:00'), ymd_hms('2021-09-25 00:01:00'), by = 'min'),
    seq(ymd_hms('2021-09-25 02:46:00'), ymd_hms('2021-09-25 03:16:00'), by = 'min'),
    seq(ymd_hms('2021-09-27 03:31:00'), ymd_hms('2021-09-27 05:31:00'), by = 'min'),
    seq(ymd_hms('2021-10-02 01:31:00'), ymd_hms('2021-10-02 09:31:00'), by = 'min'),
    seq(ymd_hms('2021-10-03 17:46:00'), ymd_hms('2021-10-03 18:31:00'), by = 'min'),
    seq(ymd_hms('2021-10-06 14:31:00'), ymd_hms('2021-10-06 15:46:00'), by = 'min'),
    seq(ymd_hms('2021-10-09 21:01:00'), ymd_hms('2021-10-09 23:46:00'), by = 'min'),
    seq(ymd_hms('2021-10-10 04:31:00'), ymd_hms('2021-10-10 05:31:00'), by = 'min'),
    seq(ymd_hms('2021-10-11 03:46:00'), ymd_hms('2021-10-11 10:16:00'), by = 'min'),
    seq(ymd_hms('2021-10-16 06:16:00'), ymd_hms('2021-10-16 06:46:00'), by = 'min'),
    seq(ymd_hms('2021-10-25 04:31:00'), ymd_hms('2021-10-25 06:31:00'), by = 'min'),
    seq(ymd_hms('2021-10-28 04:31:00'), ymd_hms('2021-10-28 12:31:00'), by = 'min'),
    seq(ymd_hms('2021-10-28 16:46:00'), ymd_hms('2021-11-29 09:31:00'), by = 'min')
),
flag = 'bad data')
dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
           DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                              TRUE ~ DO.obs),
           DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:

dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '10 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Garrison ####
ff <- grep('Gar_2021', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))
DO_comp <- dat %>%
    select(UTC, DO_GR = DO.obs) %>%
    bind_rows(DO_comp)

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(
    ymd_hms(c('2021-06-28 22:47:00', '2021-07-20 08:30:00',
    '2021-08-14 22:03:00', '2021-09-26 09:14:00',
    '2021-09-26 10:44:00', '2021-10-05 10:59:00',
    '2021-10-05 11:14:00', '2021-10-09 10:59:00',
    '2021-10-10 11:44:00', '2021-10-12 09:59:00',
    '2021-10-13 09:29:00', '2021-10-17 09:59:00'
    )),
    seq(ymd_hms('2021-06-30 14:47:00'),
        ymd_hms('2021-06-30 15:32:00'), by = 'min'),
    seq(ymd_hms('2021-07-13 08:15:00'),
        ymd_hms('2021-07-13 09:45:00'), by = 'min'),
    seq(ymd_hms('2021-09-27 09:29:00'),
    #     ymd_hms('2021-10-06 10:59:00'), by = 'min'),
    # seq(ymd_hms('2021-10-13 10:14:00'),
    #     ymd_hms('2021-10-13 10:44:00'), by = 'min'),
    # seq(ymd_hms('2021-11-13 11:14:00'),
        ymd_hms('2021-11-24 15:59:00'), by = 'min')
    ),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F)) %>%
    distinct()

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '10 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Bear Gulch (Mouth?) ####
ff <- grep('BEARMOUTH_2021', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))
DO_comp <- dat %>%
    select(UTC, DO_BM = DO.obs) %>%
    bind_rows(DO_comp)

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(
    ymd_hms('2021-09-17 20:14:00','2021-10-13 00:14:00'),
    seq(ymd_hms('2021-06-14 15:38:00'),
        ymd_hms('2021-06-15 13:53:00'), by = 'min'),
    seq(ymd_hms('2021-09-12 21:29:00'),
        ymd_hms('2021-09-12 22:44:00'), by = 'min'),
    seq(ymd_hms('2021-10-08 05:59:00'),
        ymd_hms('2021-10-08 20:59:00'), by = 'min'),
    seq(ymd_hms('2021-10-13 08:14:00'),
        ymd_hms('2021-10-13 19:14:00'), by = 'min'),
    seq(ymd_hms('2021-10-20 03:44:00'),
        ymd_hms('2021-10-20 05:29:00'), by = 'min'),
    seq(ymd_hms('2021-11-15 12:29:00'),
        ymd_hms('2021-11-24 15:44:00'), by = 'min')),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Bonita ####
ff <- grep('BONITA_2021', raw_files, value = TRUE)
dat <- read_csv(paste0('data/raw_DO/', ff))
DO_comp <- dat %>%
    select(UTC, DO_BN = DO.obs) %>%
    bind_rows(DO_comp)

dat %>% select(DO.obs) %>%
    xts(order.by = dat$UTC) %>%
    dygraph() %>%
    dySeries('DO.obs', drawPoints = TRUE) %>%
    dyRangeSelector()

flag_points <- data.frame(MST = c(
    ymd_hms('2020-08-31 15:21:00'),
    seq(ymd_hms('2021-06-14 15:49:00'),
        ymd_hms('2021-06-15 14:30:00'),
        by = 'min'),
    seq(ymd_hms('2021-09-06 04:26:00'),
        ymd_hms('2021-09-10 19:56:00'),
        by = 'min'),
    seq(ymd_hms('2021-11-04 12:56:00'),
        ymd_hms('2021-11-24 15:41:00'),
        by = 'min'),
    seq(ymd_hms('2021-08-30 02:41:00'),
        ymd_hms('2021-08-30 04:41:00'),
        by = 'min')),
    flag = 'bad data')

dt <- as.numeric(median(diff(dat$UTC)))

dat_clean <- left_join(dat, flag_points) %>%
    mutate(DO.raw = DO.obs,
        DO.obs = case_when(flag == 'bad data' ~ NA_real_,
                                   TRUE ~ DO.obs),
        DO.obs = zoo::na.approx(DO.obs, maxgap = max_gap * 60/dt, x = UTC, na.rm = F))

dat_clean %>%
    select(DO.obs, DO.raw) %>%
    xts(order.by = dat_clean$UTC) %>%
    dygraph() %>%
    dySeries('DO.raw') %>%
    dySeries('DO.obs', drawPoints = TRUE, strokeWidth = 0) %>%
    dyRangeSelector()
# align samples to consistent intervals:
dat_filled <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = 'min')) %>%
    left_join(dat_clean) %>%
    mutate(across(c(-UTC, -MST, -flag),
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -DO.sat, -flag, -MST, -Q)
# plot(dat$UTC, dat$DO.obs, type = 'l')
# points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))


# DO cross site comparison ####
DO_comp <-
    arrange(DO_comp, UTC) %>%
    mutate(across(-UTC, zoo::na.approx, na.rm = F))

ggplot(DO_comp, aes(DO_PL, DO_DL, col = UTC))+
    geom_point() +
    geom_abline(slope = 1, intercept = 0)
ggplot(DO_comp, aes(UTC, DO_PL))+
    geom_point() +
    geom_point(aes(y = DO_BM), col = 2)
ggplot(DO_comp, aes(UTC, DO_PL)) +
    geom_line() +
    geom_line(aes(y = DO_DL), col = 2)
ggplot(DO_comp, aes(DO_PL, DO_GC, col = UTC))+
    geom_point() +
    geom_abline(slope = 1, intercept = 0)
ggplot(DO_comp, aes(UTC, DO_DL)) +
    geom_line() +
    geom_line(aes(y = DO_GC), col = 2)
ggplot(DO_comp, aes(DO_GC, DO_GR, col = UTC))+
    geom_point() +
    geom_abline(slope = 1, intercept = 0)+
    ylim(0,20)
ggplot(DO_comp, aes(UTC, DO_GC)) +
    geom_line() +
    geom_line(aes(y = DO_GR), col = 2)+
    ylim(0,20)
ggplot(DO_comp, aes(DO_BM, DO_BN, col = UTC))+
    geom_point() +
    geom_abline(slope = 1, intercept = 0)
ggplot(DO_comp, aes(UTC, DO_BM)) +
    geom_line() +
    geom_line(aes(y = DO_BN), col = 2)
