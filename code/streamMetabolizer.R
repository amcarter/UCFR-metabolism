
#load all packages
library(rstan)
options(mc.cores = parallel::detectCores())
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(streamMetabolizer)
library(lubridate)
library(rnoaa)
library(rgdal)
library(sp)
library(Rcpp)
library(GSODR)
library(dataRetrieval)
options(noaakey = "RocTXcnTfTiprFMtvhljUQnloOuZtpeO")

lat<-46.790812
long<--113.748535
usgs.GC<-'12324680' # Gold Creek USGS gage
usgs.DL<-'12324200' # Deer Lodge USGS gage
usgs.PL<-'12323800' # Perkins Ln. USGS gage
usgs.BM<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita
usgs.BN<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita


# Generate temperature and O2
dat <- read.csv("312751_080421_MTT_forR.csv",header = TRUE)
names(dat)<-c("unix.time", "date", "time", "time.MST", "battery", "temp.c", "do.mgl", "do.sat","q")
#lubridate::as_datetime(dat$unix.time)
dat$date<-as.Date(dat$date,format="%m-%d-%Y")


# convert to solar time
time<-as.POSIXct(paste(dat$date, dat$time),format="%Y-%m-%d %H:%M:%S",tz="UTC")
dat$solar.time <- convert_UTC_to_solartime(
  time,
  -112.724167,
  time.type = c("apparent solar", "mean solar")
)

# Interpolate missing data
starttime <- round_date(dat$solar.time[2],"15 minutes")
endtime <- round_date(dat$solar.time[nrow(dat)-1],"15 minutes")
number <- as.numeric(difftime(endtime,starttime,15,units = c("mins")))/15+1
timeframe <- seq(starttime,endtime,length.out = number)


#interpolate data to time frame
do.inter <- approx(
  x = dat$solar.time,
  y = dat$do.mgl,
  xout = timeframe
);
temp.water <- approx(
  x = dat$solar.time,
  y = dat$temp.c,
  xout = timeframe
);


metab <- data.frame(timeframe,do.inter$y,temp.water$y)
metab$date <- as.Date(as.POSIXct(metab$timeframe))

colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)

# Get station (not sea level) pressure data from Missoula airport
station<-'727730-24153' #MSO airport
GSOD<-get_GSOD(years=2021, station=station)
pressure<-data.frame(GSOD$YEARMODA, GSOD$STP)
names(pressure)<-c("date", "press.mb")

metab<-left_join(metab,pressure)

# Generate sat
metab$DO.sat <- calc_DO_sat(
  temp=metab$temp.water,
  press=metab$press.mb,
  salinity.water = 0,
)

# Discharge

instFlow <- readNWISdata(sites = usgs.GC,
                         service = "dv", 
                         parameterCd = "00060",
                         startDate = "2021-01-01",
                         endDate = "2021-10-31") 

instFlow$dateTime <- as.Date(instFlow$dateTime)
instFlow$q.m3s<-instFlow$X_00060_00003/35.31
names(instFlow)<-c("agency", "site", "date","q.cfs","code", "tz", "q.cms")
instFlow<-select(instFlow, c(-'agency', -site, -q.cfs, -code, -tz))

metab<-merge(metab, instFlow, by="date", all.x=TRUE)

# Depth























###
meteo_nearby_stations()

ncdc_stations()

ncdc_locs(locationcategoryid='MISSOULA', sortfield='name', sortorder='desc')
out <- ncdc(datasetid='NORMAL_DLY', stationid='KMSO:USW00024153
', startdate = '2010-05-01', enddate = '2010-05-10')

out <- ncdc(datasetid='GSOM', stationid='GHCND:US1MTMS0022',  startdate = '2018-05-01', enddate = '2018-10-31', limit=500)



ncdc_plot(out, breaks="1 month", dateformat="%d/%m")