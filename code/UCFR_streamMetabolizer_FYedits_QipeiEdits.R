#Reinstall unitted and streamMetabolizer if needed
remotes::install_github('appling/unitted', force = TRUE)
remotes::install_github("USGS-R/streamMetabolizer", force = TRUE)
#Install for interacting with MesoWest data
install.packages("remotes")
remotes::install_github("fickse/mesowest")

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
#To access mesowest data
library(mesowest)
requestToken(apikey = "KyW2GDQXm7iuxMZIttsuCK15e1zWMzvOTrUpX9k3HN")


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

##Elevation of sites (masl)
ele.BN<-1103
ele.BM<-1134
ele.GC<-1275
ele.GR<-1340
ele.DL<-1378
ele.PL<-1451

##Download pressure and air temp data from DeerLodge
res<-GET("https://api.synopticdata.com/v2/stations/timeseries?stid=k38s&start=202001030000&end=202012200000&timeformat=%s&token=bc312c7c369043269a005c049228c49b")

#reorganize pressure data
data = fromJSON(rawToChar(res$content))
press.dl<-as.data.frame(data[["STATION"]][["OBSERVATIONS"]])
press.dl.unlist<-unlist(press.dl$date_time)
press.dl.unlist<-as.data.frame(press.dl.unlist)
press.dl.unlist$pressure<-unlist(press.dl$pressure_set_1d)
press.dl.unlist$press.dl.unlist<-as.POSIXct(as.numeric(press.dl.unlist$press.dl.unlist),origin='1970-01-01')
press.dl.unlist$date<-as.Date(press.dl.unlist$press.dl.unlist)

pressure_day<-press.dl.unlist %>%
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
temp.dl.unlist$date<-as.Date(temp.dl.unlist$temp.dl.unlist)

air_temp_day<-temp.dl.unlist %>%
    group_by(date) %>%
    summarise(mean.temp=mean(air_temp))

names(air_temp_day)<-c("date", "mean.temp")

##set working directory
setwd("~/Documents/GitHub/UCFR-metabolism/data")


# Deer Lodge --------------------------------------------------------------

dat <- read_csv("JAP_DO_MINIDOT_DL_2020.csv")
names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery", "temp.c", "do.mgl", "do.sat","q")

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
metab <- data.frame(dat$solar.time,dat$do.mgl,dat$temp.c)
metab$date <- as.Date(as.POSIXct(dat$solar.time))
colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)

#Add pressure data from DL to metab data frame
#pressure data is Denver time zone; Convert to solar time first.
press.dl.unlist$press.dl.unlist<- calc_solar_time(press.dl.unlist$press.dl.unlist, longitude=long)
#interpolate to sensor time
press.dl.solar<-approx(x=press.dl.unlist$press.dl.unlist, y = press.dl.unlist$pressure,
                       xout=metab$solar.time)
press.dl.solar$y<-press.dl.solar$y/100
metab<-tibble(metab,press.dl.solar$y)
metab<-rename(metab,"press.mb"="press.dl.solar$y")

#Add air temp data from DL to metab data frame
metab<-left_join(metab,air_temp_day)

# Calculate DO concentration at saturation
metab$DO.sat <- calc_DO_sat(
    temp=metab$temp.water,
    press=metab$press.mb,
    salinity.water = 0,)

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

# Depth calc based on Deer Lodge discharge
metab$depth<-0.02*metab$q.cms+0.44

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
DL_dat<-data.frame(solar.time=metab$solar.time,
                   DO.obs=metab$DO.obs,DO.sat=metab$DO.sat,depth=metab$depth,temp.water=metab$temp.water,
                   light=metab$light)

#-----------------sM model specs---------------
nb_DL <- mm_name(type='bayes', pool_K600='normal', err_obs_iid=TRUE, err_proc_iid=TRUE)

DL_specs <- specs(nb_DL,K600_daily_meanlog_meanlog=1.986474396, K600_daily_meanlog_sdlog=0.75,
                  K600_daily_sdlog_sigma=0.5, burnin_steps=500, saved_steps=500, n_cores=4, verbose=T)

#Run model
DL_fit <- metab(DL_specs, data=DL_dat)

#Check model output
DL_params<-get_params(DL_fit , uncertainty='ci')
DL_mcmc<-get_mcmc(DL_fit)
print(DL_fit)

predict_metab(DL_fit)

plot_metab_preds(DL_fit)

get_params(DL_fit)

predict_DO(DL_fit)

plot_DO_preds(DL_fit, y_var=c( "pctsat"), style='dygraphs',y_lim = list(conc = c(NA, NA), pctsat = c(NA, NA), ddodt = c(-50, 50)))

mcmc <- get_mcmc(DL_fit)
rstan::traceplot(DL_mcmc, pars='K600_daily', nrow=3)

get_fit(DL_fit)$overall %>%
    select(ends_with('Rhat'))
get_fit(DL_fit) %>%
    lapply(names)
get_fit(DL_fit)$overall %>%
    select('err_proc_iid_sigma_mean')
#launch_shinystan(DL_mcmc)
#pairs(DL_mcmc)

###-----------remove observation errors-----------------
nb_DL_obs <- mm_name(type='bayes', pool_K600='normal', err_obs_iid=FALSE, err_proc_iid=TRUE)

DL_specs_obs <- specs(nb_DL_obs,K600_daily_meanlog_meanlog=1.986474396, K600_daily_meanlog_sdlog=0.75,
                      K600_daily_sdlog_sigma=0.5, burnin_steps=500, saved_steps=500, n_cores=4, verbose=T)

#Run model
DL_fit_obs <- metab(DL_specs_obs, data=DL_dat)

#Check model output
DL_params_obs<-get_params(DL_fit_obs , uncertainty='ci')
DL_mcmc_obs<-get_mcmc(DL_fit_obs)

###----------------------------------------------------------------
###following parts aim to run DL data with mle and night regression
###Used to troubleshooting models

###mle algorithms
DL_mle <- metab_mle(data=DL_dat)
predict_metab(DL_mle)
DL_mle_params<-get_params(DL_mle)
write.csv(DL_mle_params, "DL_mle_params.csv")

###night regression to get k600
DL_night_reg <- metab_night(data=DL_dat)
predict_metab(DL_night_reg)
DL_night_reg_params<-get_params(DL_night_reg)
write.csv(DL_night_reg_params, "DL_night_reg.csv")

###import k600 from night regression to mle models
DL_night_reg_params$date<-as.Date(DL_night_reg_params$date)
DL_night_reg_params<-select(DL_night_reg_params,date,K600.daily)
DL_mle_fixk <- metab_mle(data=DL_dat,data_daily=DL_night_reg_params)
predict_metab(DL_mle_fixk)
DL_mle_fixk_params<-get_params(DL_mle_fixk)
write.csv(DL_mle_fixk_params, "DL_mle_fixk_params.csv")
###save workspace
save.image(file = "DL_troubleshoot.RData")
