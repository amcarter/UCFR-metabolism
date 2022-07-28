#load packages
library(rstan)
options(mc.cores = parallel::detectCores())
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(streamMetabolizer)
library(lubridate)
library(dataRetrieval)
library(reshape2)
library(tidyverse)
library(lme4)

##Load depth data
setwd("C:/Users/alice.carter/git/UCFR-metabolism")
site_dat <- read_csv('data/site_data.csv')
UCFR_depth<- read_csv("data/UCFR_depth_summary.csv")
UCFR_depth$date<-as.Date(UCFR_depth$date, format="%m-%d-%Y")
start.20<-as.Date("2020-07-13")
end.20<-as.Date("2020-10-20")
start.21<-as.Date("2020-06-14")
end.21<-as.Date("2020-11-01")

##OPTIONAL-Make BG and BN the same data since they are very close together
#BM.index<-which(UCFR_depth$site=="BM")
#UCFR_depth[BM.index,]$site<-"BN"

## Download average daily Discharge directly from USGS for each gage
dailyflow<-vector("list",6) # 6 = number of sites
for (i in 1:6){
    dailyflow[[i]] <- readNWISdata(sites = site_dat$nwis_code[i], #download
                         service = "dv",
                         parameterCd = "00060",
                         startDate = "2020-7-15",
                         endDate = "2021-10-31")

    dailyflow[[i]]$dateTime <- as.Date(dailyflow[[i]]$dateTime) #reformat date
    dailyflow[[i]]$q.m3s<-dailyflow[[i]]$X_00060_00003/35.31 #transform from cubic feet per second to cubic meters per second
    names(dailyflow[[i]])<-c("agency", "site", "date","q.cfs","code", "tz", "q.cms") # change column header names
    dailyflow[[i]]<-select(dailyflow[[i]],
                           c(-'agency', -'site', -'q.cfs', -'code', -'tz')) # remove unecessary data
    dailyflow[[i]]$site<-rep(site_dat$sitecode[[i]],
                             length(dailyflow[[i]]$date)) # add column with site name
}


## Turn list into data frame in long format
daily.q<-do.call(rbind.data.frame, dailyflow)

daily.q.sub<-subset(daily.q, date<end.20  & date> start.21)
write_csv(daily.q.sub, 'discharge_UCFRsites_2020.csv')

## Join discharge with depth and width data (by date)
data.sub<-left_join(daily.q.sub, UCFR_depth)
data<-left_join(daily.q, UCFR_depth)

## Make sites report in order from upstream to downstream
data$site<-factor(data$site, levels=c("PL", "DL", "GR", "GC", "BM", "BN"))
data.sub$site<-factor(data.sub$site, levels=c("PL", "DL", "GR", "GC", "BM", "BN"))
data.sum <- data.sub %>%
  group_by(site) %>%
  summarise(Min = min(q.cms,na.rm=TRUE), Max=max(q.cms,na.rm=TRUE))

## plot depth vs Q relationship by site
ggplot(data=data, aes(x=q.cms, y=depth.m, color=site))+
  geom_point(size=4)+
  theme_classic()+
  xlab("Discharge (cms)")+
  ylab("Depth (m)")+
  geom_hline(data = data.sum, aes(yintercept = Min)) +
  geom_hline(data = data.sum, aes(yintercept = Max)) +
  scale_x_continuous(limits=c(0,20))+
  scale_y_continuous(limits=c(0.3,0.83), breaks=c(0.3,0.5,0.7))+
  theme(axis.title.x=element_text(size=12,colour = "black"))+
  theme(axis.title.y=element_text(size=12,colour = "black"))+
  theme(axis.text.y=element_text(size=12,colour = "black"))+
  theme(axis.text.x=element_text(size=12,colour = "black"))+
  facet_wrap(~site, ncol=1)

ggplot(data=data, aes(x=log(q.cms), y=log(depth.m), color=site))+
  geom_abline(intercept=-1.2088,slope=0.3386, size=1.5, color='grey')+
  geom_abline(intercept=-1.003,slope=0.23633, size=1.5, color='black')+
  geom_point(size=4)+
  geom_smooth(method='lm',formula= y~x,aes(color=site),se = FALSE)+
  theme_classic()+
  ylab("log Discharge (cms)")+
  xlab("log Depth (m)")+
  #scale_x_continuous(limits=c(0,20))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))

ggplot(data=data, aes(x=date, y=q.cms))+
  geom_point(size=3)+
  theme_classic()+
  ylab("Discharge (cms)")+
  xlab("Date")+
  #scale_y_continuous(limits=c(0,25))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  facet_grid(~site, scales="free")

####Analysis
model<-lm(log(depth.m)~ log(q.cms)+ site, data=data)
model2<-lm(log(depth.m)~ log(q.cms), data=data)
summary(model)
summary(model2)

coef <- model$coefficients
fits <- tibble(site = factor(site_dat$sitecode,
                             levels = c('PL', "DL", "GR", "GC", "BM", "BN"))) %>%
    mutate(intercept = coef[1] + c(0, coef[3:7, drop = T]),
           slope = rep(coef[2], 6),
           min_q = log(sapply(!!dailyflow, function(x) min(x$q.cms, na.rm = T))),
           max_q = log(sapply(!!dailyflow, function(x) max(x$q.cms, na.rm = T))),
           min_d = intercept + slope * min_q,
           max_d = intercept + slope * max_q)

png('figures/depth_discharge_relationships.png')
    plot(log(data$q.cms), log(data$depth.m), col = data$site, pch = 19,
         xlim = c(min(fits$min_q), max(fits$max_q)),
         ylim = c(min(fits$min_d), max(fits$max_d)),
         ylab = 'log depth (m)', xlab = 'log discharge (cms)')
    for(i in 1:6){
        lines(fits[i,4:5], fits[i,6:7], col = i)
    }
    legend('topleft', legend = levels(fits$site),
           col = 1:6, pch = 19, bty = 'n')
dev.off()

fits$formula = 'log(d) = a + b * log(q)'
write_csv(fits, 'data/depth_discharge_relationships_allsites.csv')

####Individual analysis ####
data.GR<-subset(data, site=="GR")
model<-lm(log(depth.m)~ log(q.cms), data=data.GR)
summary(model)

####Individual analysis
data.GR<-subset(data, site=="DL")
model<-lm(log(depth.m)~ log(q.cms), data=data.GR)
summary(model)
