# simulate AR1 data for hierarcical model
library(forecast)
library(xts)
library(nlme)
library(tidyverse)

# prep data ####
dd <- read_csv('data/biomass_metab_model_data.csv')
dd <- dd %>%
    mutate(across(starts_with('epil_chla'),
                  function(x) x/max(dd$epil_chla_mgm2_fit, na.rm = T)),
           across(starts_with('epil_gm'),
                  function(x) x/max(dd$epil_gm2_fit, na.rm = T)),
           across(starts_with('fila_chla'),
                  function(x) x/max(dd$fila_chla_mgm2_fit, na.rm = T)),
           across(starts_with('fila_gm'),
                  function(x) x/max(dd$fila_gm2_fit, na.rm = T)))

sites <- unique(dd$site)
# fill in observation data to use as covariates in simulation
dd <- dd %>%
    mutate(year = lubridate::year(date)) %>%
    group_by(site, year) %>%
    mutate(across(-date, ~zoo::na.approx(., na.rm = F))) %>%
    filter(!is.na(GPP), !is.na(light), !is.na(fila_gm2_fit))
sim_dat <- filter(dd, site == 'GC')

# simulate data
tt <- seq(as.Date('2015-01-01'), as.Date('2023-01-10'), by = 'day')[1:2000]
n <- 2000
beta <- c(1, 3, 1)
phi <- 0.5              # ar1 parameter
sigma <- 1

light <- rep(sim_dat$light, 3) %>% sort()
light <- rep(light, 5)[1:2000]
biomass <- rep(sim_dat$fila_chla_mgm2_fit, 20)[1:2000]

X = matrix(c(light, biomass), ncol = 2)
# simulate response vector
mu <- vector(mode = 'double', length = n)
mu[1] = (beta[1] + beta[2]*light[1] + beta[3] * biomass[1])/(1-phi)

for(i in 2:n){
    mu[i] <- beta[1] + beta[2] * light[i] + beta[3] * biomass[i] +
        phi*mu[i-1] + rnorm(1, 0, sigma)
}

fit <- Arima(xts(mu, order.by = tt), order = c(1,0,0), xreg = X)
plot(tt, mu)
lines(tt, fit$fitted)
fit

# Test autoARIMA ####
dat <- dd %>% filter(site == 'GC')
X <- matrix(c(dat$light, dat$fila_chla_mgm2_fit), ncol = 2)
fit <- Arima(xts(dat$GPP, order.by = dat$date),
             order = c(1,0,0), xreg = X)

forecast::checkresiduals(fit)
plot(dat$date, dat$GPP)
lines(dat$date, fit$fitted, col = 'blue', lwd = 2)
fit
dat$fit <- fit$fitted
ggplot(dat, aes(date, GPP)) +
    geom_point() +
    geom_line( aes(y = fit), col = 'blue')+
    facet_wrap(.~year, scales = 'free_x') +
    theme_bw()
# the overall structure that seems to fit decently across sites is ARIMA(1,0,1)
