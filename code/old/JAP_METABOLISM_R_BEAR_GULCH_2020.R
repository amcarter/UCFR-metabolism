#Library all necessary packages
library(rstan)
options(mc.cores = parallel::detectCores())
library(streamMetabolizer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

#Choose and read data in CSV
dat<-read.csv(file.choose())

#convert miniDOT unix time to solar time 
dat$solar.time<-as.POSIXct(dat$solar.time, origin="1970-01-01")
dat$solar.time <- lubridate::force_tz(dat$solar.time, "UTC")
head(dat)
dat

# convert degrees to radians 
radi<-function(degrees){(degrees*pi/180)}

## estimate light for UCFR sampling site (adjust the lat and long for site)
dat$light<- calc_light(dat$solar.time, 46.69611, 113.1114) #solar.time,lat,long
head(dat)
#check data to ensure light is estimated correctly 
dat
#plot estimated light to double check 
l<-plot(dat$solar.time, dat$light)
l<-plot(dat$solar.time, dat$light)
l1<-l+scale_x_datetime(date_breaks = "1 week",labels = date_format("%W"))
#Define model and specs (Re-name files here for site )
#K600_Daily is an estimated prior  
#See slope calculations code to estimate prior meanlog for K600
nb_BG <- mm_name(type='bayes', pool_K600='normal', err_obs_iid = T, err_proc_iid =T)
BG_specs <- specs(nb_BG, K600_daily_meanlog_meanlog=2.52, K600_daily_meanlog_sdlog=0.75, K600_daily_sdlog_sigma=0.1, burnin_steps=500, 
                        saved_steps=500,verbose=T)

#data.  May include correction for DO.obs and DO.sat (/1.0)
#Also adjust average depth (depth=rep(0.7512,length(dat$water.temp)))
BG_sm<-tibble(DO.obs= dat$DO.obs/1.0, DO.sat=calc_DO_sat(dat$water.temp,855),temp.water= dat$water.temp, depth=rep(0.5,length(dat$water.temp)), light=dat$light, solar.time=dat$solar.time)
#check sm
BG_sm
#create fit
BG_fit <- metab(BG_specs, data=BG_sm, info=c(site='BEAR GULCH', source='J.A. PRATER'))

#plot fit 
plot_DO_preds(predict_DO(BG_fit))
plot_metab_preds(predict_metab(BG_fit))

#Check if model is working by comparing K values, and GPP vs ER 
#Should be no relationship between K values and GPP/ER - relationship between GPP/ER
BG_params<-get_params(BG_fit , uncertainty='ci')
BG_mcmc<-get_mcmc(BG_fit)
print(BG_fit)

mean(BG_params$K600.daily, na.rm=T)  #131
mean(BG_params$GPP.daily, na.rm=T) #0.95
mean(BG_params$ER.daily, na.rm=T) #-10.6

plot(BG_params$K600.daily,BG_params$ER.daily)
plot(BG_params$K600.daily,BG_params$GPP.daily)
plot(BG_params$GPP.daily,BG_params$ER.daily)

#Interpolate gaps in data using zoo
inter_metab <- data.frame(BG_params$date, BG_params$GPP.daily, BG_params$ER.daily, BG_params$GPP.daily.lower,
                          BG_params$GPP.daily.upper, BG_params$ER.daily.lower, BG_params$ER.daily.upper,
                          BG_params$K600.daily, BG_params$K600.daily.lower, BG_params$K600.daily.upper)

df <- read.zoo(inter_metab)

metab_final <- na.approx(df)

df2 <- data.frame(metab_final)

write.csv(BG_sm, file=file.choose(), row.names=TRUE, sep=',', 
          col.names=TRUE)
write.csv(BG_specs, file=file.choose(), row.names=TRUE, sep=',', 
          col.names=TRUE)
write.csv(df2, file=file.choose(), row.names=TRUE, sep=',', 
          col.names=TRUE)