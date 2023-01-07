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
site_dat <- read_csv('data/site_data.csv') %>%
    filter(!is.na(sitecode)) %>%
    mutate(filecode = c('PERKINS', 'DL', 'GC', 'Gar', 'BEAR.*', 'BONITA'))
depth_Q_fits <- read_csv('data/depth_discharge_relationships_allsites.csv')
bad_days <- read_csv('data/days_with_poor_DO_fits_initial_run.csv')
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
site_dat <- depth_Q_fits %>%
    rename(sitecode = site) %>%
    left_join(site_dat)

remove_bad_days = TRUE
# remove_bad_days = FALSE
enddate <- ymd_hms('2021-10-14 12:00:00')

for(i in 1:6){
    yy = data.frame()
    for(y in c(2020, 2021)){

    filename <- grep(paste0(site_dat$filecode[i], '_', y), f, value = TRUE)
    dat <- read_csv(paste0('data/prepared_data/cleaned_data/', filename))

    # convert to solar time
    time<-lubridate::as_datetime(dat$UTC)
    dat$solar.time <- convert_UTC_to_solartime(
      time,
      site_dat$Longitude[i],
      time.type = c("mean solar"))

    # clean up data
    metab <- select(dat, solar.time, DO.obs, temp.water = water.temp) %>%
        mutate(date = as.Date(solar.time))

    # Generate light
    metab$light<-calc_light(metab$solar.time,
                            latitude = site_dat$Latitude[i],
                            longitude= site_dat$Longitude[i])

    #Add pressure and air temp data from DL to metab data frame
    #pressure data is Denver time zone; Convert to solar time first.
    air_press_temp$solar.time <- convert_UTC_to_solartime(
        air_press_temp$datetime_UTC,
        longitude = site_dat$Longitude[i])
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
    metab$depth <- exp(site_dat$intercept[i] +
                           site_dat$slope[i]*log(metab$q.cms))

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

    yy <- bind_rows(yy, mod_dat)

    }

    # remove bad days
    if(remove_bad_days){
        bds <- filter(bad_days, site == site_dat$sitecode[i])

        yy <- yy %>%
            mutate(DO.obs = ifelse(substr(solar.time, 1, 13) %in%
                                       paste0(bds$bad_days, ' 12'), NA, DO.obs))
    }

    yy <- filter(yy, solar.time <= enddate)

    write_csv(yy, paste0('data/prepared_data/', site_dat$sitecode[i],
                              '2020_2021.csv'))

}


