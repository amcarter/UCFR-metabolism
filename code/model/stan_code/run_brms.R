# Script to run time series models on UCFR data
#
# library(rstan)
# options(mc.cores = parallel::detectCores())
# library(loo)
library(tidyverse)
library(brms)
library(viridis)

# read in helper functions
source('code/UCFR_helpers.R')

# prep data ####
dat <- read_csv('data/model_fits/biomass_metab_model_data.csv')

dat <- mutate(dat,
              biomass = epil_chla_mgm2_fit + fila_chla_mgm2_fit,
              log_biomass = log(epil_chla_mgm2_fit + 0.1) +
                  log(fila_chla_mgm2_fit + 0.1))
q1 <- dat %>% select(starts_with(c('epil', 'fila')), -ends_with('_se')) %>%
    mutate(across(.fns = ~case_when(. == 0 ~ NA_real_,
                                    TRUE ~ .))) %>%
    summarize(across(.fns =  ~quantile(., 0.01, na.rm = T)))
dat <- dat %>%
    mutate(year = lubridate::year(date),
           month = lubridate::month(date),
           epil_gm2_fit = epil_gm2_fit + q1$epil_gm2_fit,
           epil_chla_mgm2_fit = epil_chla_mgm2_fit + q1$epil_chla_mgm2_fit,
           fila_gm2_fit = fila_gm2_fit + q1$fila_gm2_fit,
           fila_chla_mgm2_fit = fila_chla_mgm2_fit + q1$fila_chla_mgm2_fit,
           across(starts_with(c('epil', 'fila', 'light', 'GPP')), ~ log(.),
                  .names = 'log_{.col}'),
           GPP_sd = (GPP.upper - GPP.lower)/3.92,
           log_GPP_sd = (log_GPP.upper - log_GPP.lower)/3.92)
# light = light/max(light),
# log_light = log_light - min(log_light),
# log_light = log_light/max(log_light)) %>%
dd <- dat %>%
    mutate(across(contains(c('epil', 'fila', 'light')), ~ scale(.)[,1])) %>%
    group_by(site, year) %>%
    mutate(
        across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
        site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_chla_mgm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    ungroup()

sites <- unique(dd$site)

# run on real data: ####
master_X <- select(dd, log_light,
                   epil_afdm = log_epil_gm2_fit,
                   fila_afdm = log_fila_gm2_fit,
                   epil_chla = log_epil_chla_mgm2_fit,
                   fila_chla = log_fila_chla_mgm2_fit) %>%
    mutate(across(-log_light, .fn = ~log_light * .,
                  .names = 'log_light_{.col}'))

colnames(master_X)

model_combinations <- list(1, c(1,2), c(1,3), c(1,4), c(1,5),
                           c(1,2,3), c(1,4,5),
                           c(1,2,6), c(1,3,7),
                           c(1,4,8), c(1,5,9),
                           c(1,2,3,6,7), c(1,4,5,8,9),
                           6,7,8,9, c(6,7), c(8,9))

y = arima.sim(model = list(ar = 0.5), n = 100)

N = nrow(dd)
K = 2
S = 6
ss = as.numeric(as.factor(dd$site))
new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1
y_sd = abs(rnorm(N,0,0.2))

beta = rnorm(3)
tau = 0.5
beta0 = rnorm(S, beta[1], tau)
X <- as.matrix(select(dd,log_light, log_fila_chla_mgm2_fit), ncol = K)
X = matrix(rnorm(N*K), ncol = K, nrow = N)
theta = 0.6
sigma = 1

mu = c(X %*% beta[2:(K+1)])
y = vector(length = N)
for(i in 1:length(which(new_ts == 1))){
    start = which(new_ts == 1)[i]
    if(i == length(which(new_ts == 1))){
        end = N
    } else{
        end = which(new_ts == 1)[i + 1] - 1
    }
    y[start:end] = arima.sim(model = list(ma = theta), n = end-start+1) +
        mu[start:end] + beta0[ss[start]]
}

y = rnorm(N, y, y_sd)

ma1_mod <- stan_model('code/model/stan_code/MA1.stan',
                      model_name = 'ma_mod')
arma_mod <- stan_model('code/model/stan_code/ARMA1.stan',
                       model_name = 'arma_mod')
datlist <- list(
    N = N, K = K, S = S,
    ss = ss, X = X, new_ts = new_ts,
    y = y#, y_sd = y_sd
)

ma_fit <- sampling(
    ma1_mod,
    datlist,
    iter = 4000,
    chains = 4
)

print(ma_fit, pars = c('sigma', 'theta', 'gamma', 'beta', 'tau'))
pairs(ma_fit, pars = c('sigma', 'theta', 'gamma', 'tau'))

shinystan::launch_shinystan(ma_fit)

ggplot(dd, aes(date, log_fila_chla_mgm2_fit))+
    geom_line() + facet_grid(site~year, scales = 'free')
dd %>% mutate(log_GPP = P) %>%
    ggplot( aes(log_fila_chla_mgm2_fit, log_GPP, col = site))+
    geom_point()


for(i in 1:length(model_combinations)){
    i = 1
    paste(colnames(master_X)[model_combinations[i]], sep = ' + ')
    bform <- brms::bf(GPP | mi() ~ light + ar(p = 1, gr = site) + (1|site))
    mod <- brms::brm(bform, data = dd)
    bform1 <- brms::bf(GPP | mi() ~ light + ar(p = 1, gr = site) + (1|site))
    mod1 <- brms::brm(bform1, data = dd)
    bform2 <- brms::bf(GPP | mi() ~ light + fila_chla_mgm2_fit + ar(p = 1, gr = c(site, year))
                       + 1|site + 1|year)
    mod2 <- brms::brm(bform2, data = dd)
    mccres=parll::tCs#eepatesS-# numberosieN- 200#numbeofobservtons persiteK<- 3#number ofcovariatespi-0.5eta<-rnorm(K+1)#ppulationlelefetstau<-.1bea<-rnorm(S,beta[1],tau) # randomintercept centered aroundgammasigma_proc<- 1sigma_obs <-0.1new_year<- 101# index where a newyearstartsX<-array(rnrm(N*K*S), dim =c(S,N,K))#X<- matrix(rnorm(N*K),nrow =N)theta <- vector()theta[1] <-sigma_proc/sqrt(1-phi)for(i in 2:N){if(i ==new_year){   theta[i] <- sigma_proc/sqrt(1-phi)    }else{theta[i]= phi* theta[i-1] +rnorm(1,0,sigma_proc)  }}mu<- matrix(NA, nrow=N, ncol =S)for(s in 1:S){    mu[,s] <-X[s,,] %*% beta[2:(K+1)] +theta + beta0[s]}y<- matrix(NA,nrow = N, ncol= S)for(t in 1:N){y[t,]=mvtnorm::rmvnorm(n = 1,mu[t,],sigma_obs * diag(S))}plot(theta,ylim = c(min(y),max(y)),type ='l', col= 2)for(s in 1:S){  lines(y[,s])}ar_mod<- stan_model("code/model/stan_code/AR1_multiste.stan")datlit<-list(  S=S,N=N, K =K,X=X, y =y,new_year =new_year)mod_fit <-sampling(ar_mod,           datlist)print(mod_fit, pars=c('phi', 'beta','sdp', 'sdo','beta0','tau')
