#######################################################################
# Script to simulate data to test different PI curve models of
# GPP as a function of light and biomass
#######################################################################

library(rstan)
library(tidyverse)

# prep data ####
dd <- read_csv('data/model_fits/biomass_metab_model_data.csv')

dd <- dd %>%
    mutate(year = lubridate::year(date),
           # across(starts_with(c('epil', 'fila')), ~ log(.+ 0.1),
           #        .names = 'log_{.col}'),
           # across(starts_with(c('GPP')), ~ log(.), .names = 'log_{.col}'),
           # log_GPP_sd = (log_GPP.upper - log_GPP.lower)/3.92,
           GPP_sd = (GPP.upper - GPP.lower)/3.92) %>%
    # mutate(across(contains(c('epil','fila', 'light')), ~ scale(.)[,1])) %>%
    group_by(site, year) %>%
    mutate(
        across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
        site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    ungroup()

sites <- unique(dd$site)

# Simulate data for a basic PI curve model
N <- nrow(dd)
K <- 1
light <- dd$light
# light <- dd$mean_PAR_umolm2s
biomass <- matrix(dd$fila_chla_mgm2_fit + dd$epil_chla_mgm2_fit, ncol = 1)
S <- length(sites)
ss <- as.numeric(dd$site)

#set parameters

alpha <- rnorm(K, 0.2, 0.1)
Pmax <- rnorm(K, 1, 0.1)
beta <- 4
tau_b <- 1
beta_s <- rnorm(S, beta, tau_b)
sigma <- 1


rates <- matrix(rep(NA_real_, N*K), ncol = K)
for(k in 1:K){
    rates[,k] <-  Pmax[k] * tanh((alpha[k] * light)/Pmax[k])
}
plot(light, rates[,1])
mu<- vector()
for(i in 1:N){
    mu[i] = beta_s[ss[i]]+ rates[i,]*t(biomass[i,])
}

dd$P <- rnorm(N, mu, sigma)

ggplot(dd, aes(light, P, col = fila_chla_mgm2_fit))+
    geom_point()+
    geom_point(aes(y = GPP), col = 'brown')+
    facet_grid(site~year)

PI_curve <- stan_model("code/model/stan_code/PIcurve_model_3.stan",
                       model_name = 'PI_curve')

datlist <- list(
    N = N, biomass = biomass,
    light = light, S = S, ss = ss, P = dd$P)#, P_sd = rep(P_sd, N)

fit <- sampling(
    PI_curve,
    datlist
)
plot(fit)

p1 <- plot(fit, pars = c("beta", "Pmax", "alpha",
                         "tau_b", "tau_a", "tau_P", "sigma"))
           xlim(0,2)

new_ts <- rep(0, nrow(dd))
new_ts[rle2(paste0(dd$year, dd$site))$starts] <- 1

datlist <- list(
    N = N, K = 2,
    S = S, ss = ss,
    P = dd$GPP, P_sd = dd$GPP_sd,
    light = light,
    biomass = matrix(c(dd$fila_chla_mgm2_fit, dd$epil_chla_mgm2_fit),
                     ncol = 2),
    new_ts = new_ts
)

fit_datb <- sampling(
    PI_curve,
    datlist
)

plot(fit_datb, pars = c("Pmax", "alpha", "beta", "tau_b",
                               "sigma"))
print(fit_dat, pars = c( "Pmax", "alpha",
                              "sigma"))
pairs(fit_dat, pars = c("Pmax", "alpha",
                             "sigma"))

preds <- get_model_preds(fit_datb, dd)
1-sum((preds$P_mod - preds$GPP)^2)/sum((preds$GPP - mean(preds$GPP))^2)
plot_model_fit(preds, log_ests = F)

calculate_r2_adj(fit_dat, dd, log = F)

sixtyfiftyfortythirtytwenty1nnineeightsevensixfivefourthreetwoonem((preds$GPP - preds$P_mod)^2)/sum((preds$GPP - mean(preds$GPP))^2)
