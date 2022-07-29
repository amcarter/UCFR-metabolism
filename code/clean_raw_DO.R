# Clean oxygen timeseries for UCFR to flag bad or questionable data
# A Carter
# 7/2022

# setup
library(tidyverse)
library(dygraphs)
library(xts)

setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

raw_files <- list.files('data/raw_DO/')
max_gap = 4 # maximum number of hours to fill in
# Perkins Lane ####
ff <- raw_files[7]
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
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)),
           MST = with_tz(UTC, 'MST'))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -flag)
plot(dat$UTC, dat$DO.obs, type = 'l')
points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Deer Lodge ####
ff <- raw_files[4]
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
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)),
           MST = with_tz(UTC, 'MST'))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '10 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -flag)
plot(dat$UTC, dat$DO.obs, type = 'l')
points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Gold Creek ####
ff <- raw_files[6]
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
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)),
           MST = with_tz(UTC, 'MST'))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC')

plot(dat$UTC, dat$DO.obs, type = 'l')
points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Garrison ####
ff <- raw_files[5]
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
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)),
           MST = with_tz(UTC, 'MST'))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '10 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -flag)
plot(dat$UTC, dat$DO.obs, type = 'l')
points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Bear Gulch (Mouth?) ####
ff <- raw_files[1]
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
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)),
           MST = with_tz(UTC, 'MST'))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -flag)
plot(dat$UTC, dat$DO.obs, type = 'l')
points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

# Bonita ####
ff <- raw_files[2]
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
                  function(x) zoo::na.approx(x, maxgap = max_gap * 60, x = UTC, na.rm = F)),
           MST = with_tz(UTC, 'MST'))

dat_clean <- data.frame(UTC = seq(min(dat$UTC), max(dat$UTC), by = '5 min')) %>%
    left_join(dat_filled, by = 'UTC') %>%
    select(-DO.raw, -flag)
plot(dat$UTC, dat$DO.obs, type = 'l')
points(dat_clean$UTC, dat_clean$DO.obs, pch = 20, col = 4)

write_csv(dat_clean, paste0('data/prepared_data/cleaned_data/cleaned_', ff))

