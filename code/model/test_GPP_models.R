# Compare GPP to light - are the residuals related to biomass?

library(tidyverse)
library(zoo)
library(rstan)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')
dat <- read_csv('data/metabolism_biomass_compiled_all_sites.csv')

sites <- unique(dat$site)

# check for autocorrelation ####
png('figures/GPP_acf_all_sites.png', width = 800, height = 600)
    par(mfrow = c(3,2), mar = c(2,2,1,1), oma = c(1,1,1,1))

    acf(na.approx(dat$GPP[dat$site == 'PL'], na.rm = T))
    mtext('PL', line = -2, adj = .9)
    acf(na.approx(dat$GPP[dat$site == 'DL'], na.rm = T))
    mtext('DL', line = -2, adj = .9)
    acf(na.approx(dat$GPP[dat$site == 'GR'], na.rm = T))
    mtext('GR', line = -2, adj = .9)
    acf(na.approx(dat$GPP[dat$site == 'GC'], na.rm = T))
    mtext('GC', line = -2, adj = .9)
    acf(na.approx(dat$GPP[dat$site == 'BM'], na.rm = T))
    mtext('BM', line = -2, adj = .9)
    acf(na.approx(dat$GPP[dat$site == 'BN'], na.rm = T))
    mtext('BN', line = -2, adj = .9)
dev.off()

# run AR1 model on data ####

dd <- dat %>%
    filter(site == 'PL', year == 2020) %>%
    mutate(GPP = na.approx(GPP, na.rm = F)) %>%
    filter(!is.na(GPP))
sigma_obs = median(dd$GPP.upper - dd$GPP.lower, na.rm = T)/4/sd(dd$GPP)

dd <- dd %>%
    mutate(light = (SW - mean(SW))/sd(SW),
           GPP = GPP/sd(GPP))

mod_dat <- list(P_obs = dd$GPP, N = nrow(dd),
                light = dd$light, mu_obs = sigma_obs)
fit <- stan(file = 'code/model/ar1_model_lognormal.stan',
            data = mod_dat,
            warmup = 2000, iter = 3000,
            chains = 4, cores = 4)

traceplot(fit, pars= c('phi', 'b1', 'sigma_proc', 'sigma_obs'))
plot(fit, pars= c('phi', 'b1', 'sigma_proc', 'sigma_obs'))
pairs(fit, pars= c('phi', 'b1', 'sigma_proc', 'sigma_obs'))
