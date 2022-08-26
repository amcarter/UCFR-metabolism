# prep raw O2 data from UCFR for running metabolism
# - add air pressure, temperature, and light data
# - download discharge from nearby USGS gages

#Reinstall unitted and streamMetabolizer if needed
# remotes::install_github('appling/unitted', force = TRUE)
# remotes::install_github("USGS-R/streamMetabolizer", force = TRUE)

#load all packages
library(tidyverse)
library(streamMetabolizer)
library(lubridate)
library(dataRetrieval)
library(zoo)

##set working directory
setwd("C:/Users/alice.carter/git/UCFR-metabolism")
site_dat <- read_csv('data/site_data.csv')
depth_Q_fits <- read_csv('data/depth_discharge_relationships_allsites.csv')
##Location for estimating light. Same for all sites.
lat <- 46.790812
long <- -113.748535

##Save USGS gage numbers for downloading
# usgs.GC<-'12324680' # Gold Creek USGS gage
# usgs.DL<-'12324200' # Deer Lodge USGS gage
# usgs.PL<-'12323800' # Perkins Ln. USGS gage
# usgs.GR<-'12324400' # Clark Fork ab Little Blackfoot R nr Garrison MT
# usgs.BM<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita
# usgs.BN<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita
#
# add usgs gage info to site metadata file:
# usgs_mdat <- dataRetrieval::whatNWISsites(stateCd = 'MT') %>%
#     filter(site_no %in% c(usgs.GC, usgs.DL, usgs.PL, usgs.GR, usgs.BM, usgs.BN))
# usgs_mdat <- bind_rows(usgs_mdat[1,], usgs_mdat)
# usgs_mdat <- usgs_mdat %>%
#     mutate(sitecode = c('BN', 'BM', 'GC', 'GR', 'PL', 'DL')) %>%
#     select(sitecode, nwis_code = site_no, usgs_sitename = station_nm,
#            lat = dec_lat_va, long = dec_long_va)
#
# site_dat <- left_join(site_dat, usgs_mdat[,1:3], by = 'sitecode')
#
# usgs <- sf::st_as_sf(usgs_mdat, coords = c('long', 'lat'), crs = 4326)
# sites <- site_dat %>%
#     filter(!is.na(sitecode)) %>%
#     sf::st_as_sf( coords = c('long', 'lat'), crs = 4326)
#
# write_csv(sites, 'data/site_data.csv')
# check USGS site locations relative to sites
# mapview(sites, color = 1) +
#     mapview(usgs, color = 2, cex = 5)

air_press_temp <- read_csv('data/air_pressure_temp.csv')
# prep raw data files for running metabolism -------------------------------####
f <- list.files('data/prepared_data/cleaned_data/')
site_dat$filename_2020 <- c(f[c(12, 9, 10, 11, 7, 8)])
site_dat$filename_2021 <- c(f[c(6, 3, 4, 5, 1, 2)])
site_dat <- depth_Q_fits %>%
    rename(sitecode = site) %>%
    left_join(site_dat)
for(i in 1:6){
    for(j in c(14,15)){

    dat <- read_csv(paste0('data/prepared_data/cleaned_data/', site_dat[i,j]))
    names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery",
                  "temp.c", "do.mgl", "do.sat","q")

    dat$date<-as.Date(dat$date.UTC,format="%m-%d-%Y")


    # convert to solar time
    time<-lubridate::as_datetime(dat$date.UTC)
    dat$solar.time <- convert_UTC_to_solartime(
      time,
      long,
      time.type = c("mean solar"))

    #Alternate conversion (both work)
    # convert to solar time
    posix.time.localtz <- as.POSIXct(dat$unix.time, origin='1970-01-01', tz='UTC')
    lubridate::tz(posix.time.localtz)

    dat$solar.time <- convert_UTC_to_solartime (
      posix.time.localtz,
      long,
      time.type="mean solar")

    # clean up data
    metab <- data.frame(dat$solar.time, dat$do.mgl,dat$temp.c)
    metab$date <- as.Date(as.POSIXct(dat$solar.time))
    colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

    # Generate light
    metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)

    #Add pressure and air temp data from DL to metab data frame
    #pressure data is Denver time zone; Convert to solar time first.
    air_press_temp$solar.time <- calc_solar_time(air_press_temp$datetime_MST,
                                                 longitude=long)
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

    # Discharge

    instFlow <- readNWISdata(sites = site_dat$nwis_code[i],
                             service = "dv",
                             parameterCd = "00060",
                             startDate = "2020-05-01",
                             endDate = "2021-11-30")

    instFlow$dateTime <- as.Date(instFlow$dateTime)
    instFlow$q.m3s<-instFlow$X_00060_00003/35.31
    names(instFlow)<-c("agency", "site", "date","q.cfs","code", "tz", "q.cms")
    instFlow<-select(instFlow, c(-'agency', -site, -q.cfs, -code, -tz))

    metab<-merge(metab, instFlow, by="date", all.x=TRUE)

    # Depth calc based on Deer Lodge discharge
    # metab$depth<-0.02*metab$q.cms+0.44
    metab$depth<-exp(site_dat$intercept[i] + site_dat$slope[i]*log(metab$q.cms))

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
    mod_dat<-data.frame(solar.time=metab$solar.time,
                       DO.obs=metab$DO.obs, DO.sat=metab$DO.sat,
                       depth=metab$depth,discharge = metab$q.cms,
                       temp.water=metab$temp.water,
                       light=metab$light)

    # remove leading and ending NAs
    w <- which(!is.na(mod_dat$DO.obs))
    mod_dat <- mod_dat[min(w):max(w),]


    if(j == 14){
    write_csv(mod_dat, paste0('data/prepared_data/', site_dat$sitecode[i], '_2020.csv'))
    } else{
        write_csv(mod_dat, paste0('data/prepared_data/', site_dat$sitecode[i], '_2021.csv'))
    }

    }
}


