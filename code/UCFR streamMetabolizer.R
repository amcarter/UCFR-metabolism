
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
library(shinystan)
library(tidyverse)
options(noaakey = "RocTXcnTfTiprFMtvhljUQnloOuZtpeO")

##Location for estimating light. Same for all sites. 
lat<-46.790812
long<--113.748535

##Save USGS gage numbers for downloading
usgs.GC<-'12324680' # Gold Creek USGS gage
usgs.DL<-'12324200' # Deer Lodge USGS gage
usgs.PL<-'12323800' # Perkins Ln. USGS gage
usgs.GR<-'12324400' # Clark Fork ab Little Blackfoot R nr Garrison MT
usgs.BM<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita
usgs.BN<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita

##Save MSO airport ID for downloading pressure
station<-'727730-24153' #MSO airport

##set working directory
setwd("~/GitHub/UCFR-metabolism/data")


# Deer Lodge --------------------------------------------------------------

dat <- read_csv("JAP_DO_MINIDOT_DL_2020.csv")
names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery", "temp.c", "do.mgl", "do.sat","q")

dat$date<-as.Date(dat$date.UTC,format="%m-%d-%Y")


# convert to solar time
time<-lubridate::as_datetime(dat$unix.time)
dat$solar.time <- convert_UTC_to_solartime(
  time,
  long,
  time.type = c("apparent solar", "mean solar")
)
#Alt: as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")

# clean up
metab <- data.frame(dat$solar.time,dat$do.mgl,dat$temp.c)
metab$date <- as.Date(as.POSIXct(dat$solar.time))
colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)


# Get station (not sea level) pressure data from Missoula airport
GSOD<-get_GSOD(years=2020, station=station)
pressure<-data.frame(GSOD$YEARMODA, GSOD$STP)
names(pressure)<-c("date", "press.mb")

metab<-left_join(metab,pressure)

# Calculate DO concentration at saturation
metab$DO.sat <- calc_DO_sat(
  temp=metab$temp.water,
  press=metab$press.mb,
  salinity.water = 0,
)

# Discharge

instFlow <- readNWISdata(sites = usgs.DL,
                         service = "dv", 
                         parameterCd = "00060",
                         startDate = "2020-05-01",
                         endDate = "2020-10-31") 

instFlow$dateTime <- as.Date(instFlow$dateTime)
instFlow$q.m3s<-instFlow$X_00060_00003/35.31
names(instFlow)<-c("agency", "site", "date","q.cfs","code", "tz", "q.cms")
instFlow<-select(instFlow, c(-'agency', -site, -q.cfs, -code, -tz))

metab<-merge(metab, instFlow, by="date", all.x=TRUE)

# Depth- assign depth =1m so we can look at how different depths affect rates post modeling.
metab$depth<-rep(1, times=length(metab$date))


DL_dat<-tibble(solar.time=metab$solar.time,DO.obs=metab$DO.obs,DO.sat=metab$DO.sat,depth=metab$depth,temp.water=metab$temp.water,light=metab$light)

#Run model
nb_DL <- mm_name(type='bayes', pool_K600='normal', err_obs_iid=TRUE, err_proc_iid=TRUE)

DL_specs <- specs(nb_DL,K600_daily_meanlog_meanlog=1.986474396, K600_daily_meanlog_sdlog=0.75, 
                  K600_daily_sdlog_sigma=0.1, burnin_steps=500, saved_steps=500,verbose=T)

DL_fit <- metab(DL_specs, data=DL_dat)

#Check model output
DL_params<-get_params(DL_fit , uncertainty='ci')
DL_mcmc<-get_mcmc(DL_fit)
print(DL_fit)

launch_shinystan(DL_mcmc)
pairs(DL_mcmc)
# Garrison ----------------------------------------------------------------
rm(list = ls())
dat <- read_csv("JAP_DO_MINIDOT_GAR_2020.csv")
names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery", "temp.c", "do.mgl", "do.sat","q")
#lubridate::as_datetime(dat$unix.time)
dat$date<-as.Date(dat$date.UTC,format="%m-%d-%Y")

dat<-subset(dat,dat$date != as.Date("2020-08-02") & dat$date != as.Date("2020-08-03"))

# convert to solar time
time<-as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")
dat$solar.time <- convert_UTC_to_solartime(
  time,
  long,
  time.type = c("apparent solar", "mean solar")
)



metab<-dat[c(10,6,5,9)]

colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)


# Get station (not sea level) pressure data from Missoula airport
GSOD<-get_GSOD(years=2020, station=station)
pressure<-data.frame(GSOD$YEARMODA, GSOD$STP)
names(pressure)<-c("date", "press.mb")

metab<-left_join(metab,pressure)

# Calculate DO concentration at saturation
metab$DO.sat <- calc_DO_sat(
  temp=metab$temp.water,
  press=metab$press.mb,
  salinity.water = 0,
)

# Discharge

instFlow <- readNWISdata(sites = usgs.GR,
                         service = "dv", 
                         parameterCd = "00060",
                         startDate = "2020-05-01",
                         endDate = "2020-10-31") 

instFlow$dateTime <- as.Date(instFlow$dateTime)
instFlow$q.m3s<-instFlow$X_00060_00003/35.31
names(instFlow)<-c("agency", "site", "date","q.cfs","code", "tz", "q.cms")
instFlow<-select(instFlow, c(-'agency', -site, -q.cfs, -code, -tz))

metab<-left_join(metab, instFlow)

# Depth- assign depth =1m so we can look at how different depths affect rates post modeling.
metab$depth<-rep(1, times=length(metab$date))


GAR_dat<-tibble(solar.time=metab$solar.time,DO.obs=metab$DO.obs,DO.sat=metab$DO.sat,depth=metab$depth,temp.water=metab$temp.water,light=metab$light)

#Run model
nb_GAR<- mm_name(type='bayes', pool_K600='normal', err_obs_iid=TRUE, err_proc_iid=TRUE)

GAR_specs <- specs(nb_GAR, K600_daily_meanlog_meanlog=1.986474396, K600_daily_meanlog_sdlog=0.75, K600_daily_sdlog_sigma=0.1, burnin_steps=500, saved_steps=500,verbose=T)

GAR_fit <- metab(specs=GAR_specs, data=GAR_dat)

GAR_params<-get_params(GAR_fit , uncertainty='ci')
GAR_mcmc<-get_mcmc(GAR_fit)
launch_shinystan(GAR_mcmc)


GAR_fit

plot_DO_preds(GAR_fit)
get_fit(GAR_fit)$overall %>%
  select(ends_with('Rhat'))

rstan::traceplot(GAR_mcmc, pars='err_proc_iid_sigma', nrow=3)




###OLD

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





meteo_nearby_stations()

ncdc_stations()

ncdc_locs(locationcategoryid='MISSOULA', sortfield='name', sortorder='desc')
out <- ncdc(datasetid='NORMAL_DLY', stationid='KMSO:USW00024153
', startdate = '2010-05-01', enddate = '2010-05-10')

out <- ncdc(datasetid='GSOM', stationid='GHCND:US1MTMS0022',  startdate = '2018-05-01', enddate = '2018-10-31', limit=500)



ncdc_plot(out, breaks="1 month", dateformat="%d/%m")