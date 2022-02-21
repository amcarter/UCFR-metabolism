
# Beginning ---------------------------------------------------------------


#load all packages
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
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
library(dygraphs)
library(httr)
library(jsonlite)
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

##Save weather station ID's for downloading pressure
mso.station<-'727730-24153' #MSO airport
dl.station<-'GSOD:USC00242275' #Deer lodge


##Elevation of sites (masl)
ele.BN<-1103
ele.BM<-1134
ele.GC<-1275
ele.GR<-1340
ele.DL<-1378
ele.PL<-1451


##Download air temp data from DeerLodge (used for all sites for now)
res<-GET("https://api.synopticdata.com/v2/stations/timeseries?stid=drlm&start=202001030000&end=202012200000&timeformat=%s&token=bc312c7c369043269a005c049228c49b")

#reorganize
data = fromJSON(rawToChar(res$content))
temp.dl<-as.data.frame(data[["STATION"]][["OBSERVATIONS"]])
temp.dl.unlist<-unlist(temp.dl$date_time)
temp.dl.unlist<-as.data.frame(temp.dl.unlist)
temp.dl.unlist$air_temp<-unlist(temp.dl$air_temp_set_1)
temp.dl.unlist$temp.dl.unlist<-as.POSIXct(as.numeric(temp.dl.unlist$temp.dl.unlist),origin='1970-01-01')
temp.dl.unlist$date<-as.Date(temp.dl.unlist$temp.dl.unlist)

air_temp_day<-temp.dl.unlist %>%
  group_by(date) %>% 
  summarise(mean.temp=mean(air_temp))


##set working directory
setwd("~/GitHub/UCFR-metabolism/data")


# Deer Lodge --------------------------------------------------------------

dat <- read_csv("JAP_DO_MINIDOT_DL_2020.csv")
names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery", "temp.c", "do.mgl", "do.sat","q")

dat$date<-as.Date(dat$date.UTC,format="%m-%d-%Y")
#<-as_datetime(dat$unix.time)

dat<-subset(dat,dat$date != as.Date("2020-08-24") & dat$date != as.Date("2020-08-25"))



# convert to solar time
posix.time.localtz <- as.POSIXct(dat$unix.time, origin='1970-01-01', tz='UTC')
lubridate::tz(posix.time.localtz)

dat$solar.time <- convert_UTC_to_solartime (
  posix.time.localtz,
  long, 
  time.type="mean solar"
)

#Alt: as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")
#lubridate::as_datetime(dat$unix.time)


# clean up
metab <- data.frame(dat$solar.time,dat$do.mgl,dat$temp.c)
metab$date <- as.Date(as.POSIXct(dat$solar.time))
colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)


# Get station (not sea level) pressure data from Missoula airport
GSOD<-get_GSOD(years=2020, station=mso.station)
pressure<-data.frame(GSOD$YEARMODA, GSOD$STP)
names(pressure)<-c("date", "press.mb")

pressure<-left_join(pressure,air_temp_day)


# Correct BP for the elveation at the site
## Function to correct barometric pressure for altitude. From Colt (2012).
## This function gives bp in mmHg for altitude given nearby measurement of standardized barometric pressure. 
## temp is degC 
## alt is m
## bpst is in inches of Hg and the sort of data you will find from U.S. weather stations.  Delete the *25.4 if you in the rest of the world where pressure is given in metric units

bpcalc<- function(bpst, alt, degC) {
  bpst*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+degC)))
}

pressure$press.st.mb<-bpcalc((pressure$press.mb), ele.DL,pressure$mean.temp )

metab<-left_join(metab,pressure)

# Calculate DO concentration at saturation
metab$DO.sat <- calc_DO_sat(
  temp=metab$temp.water,
  press=metab$press.st.mb,
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
metab$depth<-0.02*metab$q.cms+0.44
  #rep(1, times=length(metab$date))

dat$solar.time<- as.POSIXct(dat$solar.time, 
                            tz= "UTC", 
                            tryFormats = c("%Y-%m-%d %H:%M:%OS",
                                           "%Y/%m/%d %H:%M:%OS",
                                           "%Y-%m-%d %H:%M",
                                           "%Y/%m/%d %H:%M",
                                           "%Y-%m-%d",
                                           "%Y/%m/%d",
                                           "%m/%d/%y %H:%M"))

DL_dat<-data.frame(solar.time=metab$solar.time,
               DO.obs=metab$DO.obs,DO.sat=metab$DO.sat,depth=metab$depth,temp.water=metab$temp.water,
               light=metab$light)#, discharge=metab$q.cms)
#write.csv(DL_dat,"example.csv")



#Run model
bayes_name <- mm_name(type='bayes', pool_K600='normal',  err_obs_iid=FALSE, err_proc_iid=TRUE)#GPP_fun="satlight",

night<-metab_night(data=DL_dat)
DL_params<-get_params(night)


bayes_specs <- specs(bayes_name, burnin_steps=1000, saved_steps=200, n_cores=4, verbose=T,
                K600_daily_meanlog_meanlog=1.986474396, K600_daily_meanlog_sdlog=.75, K600_daily_sdlog_sigma=0.5)
                 #alpha_meanlog=0.3, alpha_sdlog=1, Pmax_mu=10, Pmax_sigma=1)

DL_fit_noobs<- metab(bayes_specs, data=DL_dat)

#Check model output
DL_params<-get_params(DL_fit_noobs , uncertainty='ci')
DL_mcmc<-get_mcmc(DL_fit_noobs)
print(DL_fit_normalk)

plot_DO_preds(DL_fit_noobs, y_var=c( "pctsat"), style='dygraphs',y_lim = list(conc = c(NA, NA), pctsat = c(NA, NA), ddodt = c(-50, 50)))#ddodt,conc,pctsat

#write_rds(p, "C:/Users/matt/Documents/GitHub/UCFR-metabolism/data/check_DO_sat.RDS")
#write_rds(p2, "C:/Users/matt/Documents/GitHub/UCFR-metabolism/data/check_dDOdt.RDS")
#write_rds(p3, "C:/Users/matt/Documents/GitHub/UCFR-metabolism/data/check_DO_conc.RDS")



launch_shinystan(DL_mcmc)
pairs(DL_mcmc)
plot(DL_params$ER.daily, DL_params$K600.daily)
model<-lm(DL_params$ER.daily~ DL_params$K600.daily)
summary(model)
get_fit(DL_fit_normalk)$overall %>%
  select(ends_with('Rhat'))


# Garrison ----------------------------------------------------------------
#rm(list = ls())
dat <- read_csv("JAP_DO_MINIDOT_GAR_2020.csv")
names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery", "temp.c", "do.mgl", "do.sat","q")
dat$date<-as.Date(dat$date.UTC,format="%m-%d-%Y")

#dat<-subset(dat,dat$date != as.Date("2020-08-02") & dat$date != as.Date("2020-08-03"))

# convert to solar time
time<-lubridate::as_datetime(dat$unix.time)
  #as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")
dat$solar.time <- convert_UTC_to_solartime(
  time,
  long,
  time.type = c( "apparent solar")
)

#Begin cleaning up
metab<-dat[c(10,6,5,9)]

colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = 47, longitude= 113)


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

# Depth- assign depth based on site specific equation so we can look at how different depths affect rates post modeling.
metab$depth<-0.04*metab$q.cms+0.299
  #rep(1, times=length(metab$date))


GAR_dat<-data.frame(solar.time=metab$solar.time,DO.obs=metab$DO.obs,DO.sat=metab$DO.sat,depth=metab$depth,temp.water=metab$temp.water,light=metab$light
, discharge=metab$q.cms)#Add discharge if binned or linear pooling.


ggplot(data=GAR_dat, aes(x=solar.time, y=light))+
  geom_point(color="black", size=3)+
  #geom_line(aes(x=solar.time, y=discharge),color="blue")+
  #geom_point(aes(y=y, x=time),color="red", size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))
  

#Run model
nb_GAR<- mm_name(type='bayes', pool_K600='linear', err_obs_iid=TRUE, err_proc_iid=TRUE)

GAR_specs <- specs(nb_GAR, burnin_steps=500, saved_steps=500,verbose=T)

GAR_fit <- metab(specs=GAR_specs, data=GAR_dat)

GAR_params<-get_params(GAR_fit , uncertainty='ci')

ggplot(data=GAR_params, aes(x=date, y=K600.daily))+
  geom_point(color="black", size=3)+
#geom_point(aes(y=y, x=time),color="red", size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))

GAR_mcmc<-get_mcmc(GAR_fit)
launch_shinystan(GAR_mcmc)


GAR_fit

plot_DO_preds(GAR_fit)
get_fit(DL_fit)$overall %>%
  select(ends_with('Rhat'))

rstan::traceplot(DL_mcmc, pars='elp___', nrow=3)

plot(GAR_params$ER.daily, GAR_params$K600.daily)
#saveRDS(GAR_fit, file = "GAR_fit_normalpooling.RDS")












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



date_table <- table(DL_dat$solar.time)
num_dates <- length(date_table)
num_daily_obs <- unique(unname(date_table))
if(length(num_daily_obs) > 1) {
  warning(paste0(
    sapply(num_daily_obs, function(ndo) {
      tslabel <- paste(ndo, 'rows per day')
      tsdates <- names(date_table)[date_table == ndo]
      paste0(tslabel, ': ', paste0(tsdates, collapse=', '))
    }),
    collapse='\n')
  )
  stop("dates have differing numbers of rows; observations cannot be combined in matrix")
  



#Fixing K
  nb_DL <- mm_name(type='bayes', pool_K600='normal_sdzero', err_obs_iid=TRUE, err_proc_iid=TRUE)
  
  DL_specs <- specs(nb_DL,K600_daily_meanlog_meanlog=log(11), K600_daily_meanlog_sdlog=0.1,burnin_steps=500, saved_steps=500,verbose=T)
  









