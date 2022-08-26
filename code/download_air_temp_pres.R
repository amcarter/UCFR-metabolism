# Download air temperature and pressure for UCFR sites

#load all packages
library(tidyverse)
library(httr)
library(jsonlite)

##set working directory
setwd("C:/Users/alice.carter/git/UCFR-metabolism")

##Download pressure and air temp data from DeerLodge
# for 2020:
res<-GET("https://api.synopticdata.com/v2/stations/timeseries?stid=k38s&start=202001030000&end=202112200000&timeformat=%s&token=bc312c7c369043269a005c049228c49b")

#reorganize pressure data
data = fromJSON(rawToChar(res$content))
press.dl<-as.data.frame(data[["STATION"]][["OBSERVATIONS"]])
press.dl.unlist<-unlist(press.dl$date_time)
press.dl.unlist<-as.data.frame(press.dl.unlist)
press.dl.unlist$pressure<-unlist(press.dl$pressure_set_1d)
press.dl.unlist$press.dl.unlist<-as.POSIXct(as.numeric(press.dl.unlist$press.dl.unlist),origin='1970-01-01')
press.dl.unlist$date<-as.Date(press.dl.unlist$press.dl.unlist, tz = 'MST')

pressure_day <- press.dl.unlist %>%
  group_by(date) %>%
  summarise(mean.press=mean(pressure/100))

names(pressure_day)<-c("date", "press.mb")

#reorganize air temp data
data = fromJSON(rawToChar(res$content))
temp.dl<-as.data.frame(data[["STATION"]][["OBSERVATIONS"]])
temp.dl.unlist<-unlist(temp.dl$date_time)
temp.dl.unlist<-as.data.frame(temp.dl.unlist)
temp.dl.unlist$air_temp<-unlist(temp.dl$air_temp_set_1)
temp.dl.unlist$temp.dl.unlist<-as.POSIXct(as.numeric(temp.dl.unlist$temp.dl.unlist),origin='1970-01-01')
temp.dl.unlist$date<-as.Date(temp.dl.unlist$temp.dl.unlist, tz = 'MST')

air_temp_day<-temp.dl.unlist %>%
  group_by(date) %>%
  summarise(mean.temp=mean(air_temp))

names(air_temp_day)<-c("date", "mean.temp")
air_pres_temp_daily <- full_join(pressure_day, air_temp_day) %>%
    rename(temp_C = mean.temp)

air_pres_temp <- temp.dl.unlist %>%
    rename(press.dl.unlist = temp.dl.unlist) %>%
    full_join(press.dl.unlist, by = c('press.dl.unlist', 'date')) %>%
    mutate(press.mb = pressure/100) %>%
    select(datetime_MST = press.dl.unlist, date, press.mb, air_temp)

write_csv(air_pres_temp, 'data/air_pressure_temp.csv')
