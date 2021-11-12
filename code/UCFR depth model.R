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

##Load depth data
setwd("~/GitHub/UCFR-metabolism/data")
UCFR_depth<- read_csv("UCFR_depth_summary.csv")
UCFR_depth$date<-as.Date(UCFR_depth$date, format="%m-%d-%Y")

##Save USGS gage numbers for downloading
usgs.GC<-'12324680' # Gold Creek USGS gage
usgs.DL<-'12324200' # Deer Lodge USGS gage
usgs.PL<-'12323800' # Perkins Ln. USGS gage
usgs.GR<-'12324400' # Clark Fork ab Little Blackfoot R nr Garrison MT
usgs.BM<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita
usgs.BN<-'12331800' # Near Drummond, but pretty close to Bear Mouth and Bonita

##Turb into data frame
gage.id<-c(usgs.GC,usgs.DL,usgs.PL,usgs.BM,usgs.BN, usgs.GR)
gage.name<-c("GC", "DL", "PL", "BM","BN", "GR")
USGS.gage<-data.frame(gage.id,gage.name)

## Download average daily Discharge directly from USGS for each gage
dailyflow<-vector("list",6) # 6 = number of sites
for (i in 1:6){
dailyflow[[i]] <- readNWISdata(sites = USGS.gage$gage.id[i], #download
                         service = "dv", 
                         parameterCd = "00060",
                         startDate = "2020-7-15",
                         endDate = "2021-10-31") 

dailyflow[[i]]$dateTime <- as.Date(dailyflow[[i]]$dateTime) #reformat date
dailyflow[[i]]$q.m3s<-dailyflow[[i]]$X_00060_00003/35.31 #transform from cubic feet per second to cubic meters per second
names(dailyflow[[i]])<-c("agency", "site", "date","q.cfs","code", "tz", "q.cms") # change column header names
dailyflow[[i]]<-select(dailyflow[[i]], c(-'agency', -'site', -'q.cfs', -'code', -'tz')) # remove unecessary data
dailyflow[[i]]$site<-rep(USGS.gage$gage.name[[i]], length(dailyflow[[i]]$date)) # add column with site name
}

## Turn list into data frame in long format
daily.q<-do.call(rbind.data.frame, dailyflow)

## Join discharge with depth and width data (by date)
data<-left_join(daily.q, UCFR_depth )

## Make sites report in order from upstream to downstream
data$site<-factor(data$site, levels=c("PL", "DL", "GR", "GC", "BM", "BN"))

dMin <- data %>%
  group_by(site) %>%
  summarise(Min = min(q.cms), Max=max(q.cms))

## plot depth vs Q relationship by site
ggplot(data=data, aes(x=depth.m, y=q.cms, color=site))+
  geom_point(size=3)+
  theme_classic()+
  ylab("Discharge (cms)")+
  xlab("Depth (m)")+
  scale_y_continuous(limits=c(0,20))+
  theme(axis.title.x=element_text(size=18,colour = "black"))+
  theme(axis.title.y=element_text(size=18,colour = "black"))+
  theme(axis.text.y=element_text(size=18,colour = "black"))+
  theme(axis.text.x=element_text(size=18,colour = "black"))+
  facet_grid(~site, scales="free")

ggplot(data=data, aes(x=depth.m, y=q.cms, color=site))+
  geom_point(size=3)+
  theme_classic()+
  ylab("Discharge (cms)")+
  xlab("Depth (m)")+
  scale_y_continuous(limits=c(0,20))+
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

