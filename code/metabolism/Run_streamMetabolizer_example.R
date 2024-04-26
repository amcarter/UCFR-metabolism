# Run streamMetabolizer on prepared datasets

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
library(tidyverse)
library(dygraphs)
library(dataRetrieval)
library(zoo)

# prep raw O2 data from UCFR for running metabolism
# - add air pressure, temperature, and light data
# - download discharge from nearby USGS gages

##set working directory

##Save USGS gage number and site lat, long for downloading data
usgs.gage<-'12323800' # Perkins Ln. USGS gage
Latitude <- 46.2087
Longitude <- -112.7700

# load in air pressure data (downloaded from local airport)
air_press_temp <- read_csv('air_pressure_temp.csv')


# prep raw data file for running metabolism -------------------------------####
dat <- read_csv('cleaned_DO_data.csv')

# convert to solar time
time<-lubridate::as_datetime(dat$UTC)
dat$solar.time <- convert_UTC_to_solartime(
    time,
    Longitude,
    time.type = c("mean solar"))

# clean up data
metab <- select(dat, solar.time, DO.obs, temp.water = water.temp) %>%
    mutate(date = as.Date(solar.time))

# Generate light
metab$light<-calc_light(metab$solar.time,
                        latitude = Latitude,
                        longitude= Longitude)

#Add pressure and air temp data from DL to metab data frame
#pressure data is Denver time zone; Convert to solar time first.
air_press_temp$solar.time <- convert_UTC_to_solartime(
    air_press_temp$datetime_UTC,
    longitude = Longitude)
#interpolate to sensor time
press.dl.solar<-approx(x = air_press_temp$solar.time,
                       y = air_press_temp$press.mb,
                       xout = metab$solar.time)
temp.dl.solar<-approx(x = air_press_temp$solar.time,
                      y = air_press_temp$air_temp,
                      xout = metab$solar.time)

metab<-tibble(metab, press.dl.solar$y, temp.dl.solar$y)
metab<-rename(metab, "press.mb"="press.dl.solar$y",
              "air_temp"="temp.dl.solar$y")

# Calculate DO concentration at saturation
metab$DO.sat <- calc_DO_sat(
    temp=metab$temp.water,
    press=metab$press.mb,
    salinity.water = 0,)

# Discharge from USGS gage
instFlow <- readNWISdata(sites = usgs.gage,
                         service = "dv",
                         parameterCd = "00060",
                         startDate = as.Date(dat$solar.time[1]),
                         endDate = as.Date(dat$solar.time[nrow(dat)]))

instFlow$dateTime <- as.Date(instFlow$dateTime)
# convert cfs to cms
instFlow$q.m3s<-instFlow$X_00060_00003/35.31
names(instFlow)<-c("agency", "site", "date","q.cfs","code", "tz", "q.cms")
instFlow<-select(instFlow, c(-'agency', -site, -q.cfs, -code, -tz))

metab<-merge(metab, instFlow, by="date", all.x=TRUE)

# Depth calc based on depth by discharge relationship generated from field depth measurements
# log(depth) = -0.962 + 0.35*log(discharge)
metab$depth <- exp(-0.962 + 0.35 * log(metab$q.cms))

#Format date before model
dat$solar.time<- as.POSIXct(dat$solar.time,
                            tz= "UTC",
                            tryFormats = c("%Y-%m-%d %H:%M:%OS",
                                           "%Y/%m/%d %H:%M:%OS",
                                           "%Y-%m-%d %H:%M",
                                           "%Y/%m/%d %H:%M",
                                           "%Y-%m-%d",
                                           "%Y/%m/%d",
                                           "%m/%d/%y %H:%M"))

#Compile necessary parts for sM model
mod_dat <- select(metab,
                  solar.time, DO.obs, DO.sat,
                  depth, discharge = q.cms,
                  temp.water, light)

# remove leading and ending NAs
w <- which(!is.na(mod_dat$DO.obs))
mod_dat <- mod_dat[min(w):max(w),]

# Make the time steps 15 minutes
mod_dat <- data.frame(solar.time = seq(min(mod_dat$solar.time),
                                       max(mod_dat$solar.time),
                                       by = '15 min')) %>%
    left_join(mod_dat, by = 'solar.time')


# Model Metabolism using cleaned and prepped data: ####
# Visualize the data #####
head(mod_dat) ; tail(mod_dat)

mod_dat %>% unitted::v() %>%
    mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
    select(solar.time, starts_with('DO')) %>%
    gather(type, DO.value, starts_with('DO')) %>%
    mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
    ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() +
    facet_grid(units ~ ., scale='free_y') + theme_bw() +
    scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)',
            light='PAR\n(umol m^-2 s^-1)', discharge='Q\n(cms)')
mod_dat %>% unitted::v() %>%
    mutate(discharge = log(discharge)) %>%
    select(solar.time, depth, temp.water, light, discharge) %>%
    gather(type, value, depth, temp.water, light, discharge) %>%
    mutate(
        type=ordered(type, levels=c('depth','temp.water','light','discharge')),
        units=ordered(labels[type], unname(labels))) %>%
    ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() +
    facet_grid(units ~ ., scale='free_y') + theme_bw() +
    scale_color_discrete('variable')

#sM model specs--------------- ####
bayes_name <- mm_name(type='bayes', pool_K600='normal',
                      err_obs_iid=TRUE, err_proc_iid=TRUE)

#
bayes_specs <- specs(bayes_name,
                     burnin_steps = 1000,
                     saved_steps = 1000,
                     K600_daily_meanlog_meanlog = 2.484906649788,
                     K600_daily_meanlog_sdlog = 0.75,
                     K600_daily_sdlog_sigma = 0.001,
                     n_cores=4, verbose=T)


#Run model
fit <- metab(bayes_specs, data=mod_dat)

