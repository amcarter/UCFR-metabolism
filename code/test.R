
# Test ---------------------------------------------------------------


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


setwd("~/GitHub/UCFR-metabolism/data")
DL_dat<-read_csv("DL_example.csv")
#Run model
bayes_name <- mm_name(type='bayes', pool_K600='normal', err_obs_iid=TRUE, err_proc_iid=TRUE)

bayes_specs <- specs(bayes_name, burnin_steps=100, saved_steps=200, n_cores=1, verbose=T,
                 K600_daily_meanlog_meanlog=1.986474396, K600_daily_meanlog_sdlog=0.75, K600_daily_sdlog_sigma=0.1)

DL_fit_normalk <- metab_bayes(bayes_specs, data=DL_dat)




#Check model output
DL_params<-get_params(DL_fit_medk , uncertainty='ci')
DL_mcmc<-get_mcmc(DL_fit_medk)
print(DL_fit_normalk)

launch_shinystan(DL_mcmc)
pairs(DL_mcmc)
plot(DL_params$ER.daily, DL_params$K600.daily)
model<-lm(DL_params$ER.daily~ DL_params$K600.daily)
summary(model)
