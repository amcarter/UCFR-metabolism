# Compare GPP to light - are the residuals related to biomass?

library(tidyverse)
library(zoo)
library(rstan)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')
dat <- read_csv('data/metabolism_biomass_compiled_all_sites.csv') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

sites <- unique(dat$site)

# check for autocorrelation ####
png('figures/GPP_acf_all_sites.png', width = 800, height = 600)
    par(mfrow = c(3,2), mar = c(2,2,1,1), oma = c(1,1,1,1))

    acf(dat$GPP[dat$site == 'PL' & dat$year == 2020], na.action = na.pass,
        xlim = c(0, 20), ylim = c(-0.2, 1), lwd = 2)
    par(new = T)
    acf(dat$GPP[dat$site == 'PL' & dat$year == 2021], na.action = na.pass,
        col = 2, xlim = c(-0.15, 19.85), ylim = c(-0.2, 1), xaxt = 'n', lwd = 2)
    mtext('PL', line = -2, adj = .9)
    acf(dat$GPP[dat$site == 'DL' & dat$year == 2020], na.action = na.pass,
        xlim = c(0, 20), ylim = c(-0.2, 1), lwd = 2)
    par(new = T)
    acf(dat$GPP[dat$site == 'DL' & dat$year == 2021], na.action = na.pass,
        col = 2, xlim = c(-0.15, 19.85), ylim = c(-0.2, 1), xaxt = 'n', lwd = 2)
    mtext('DL', line = -2, adj = .9)
    acf(dat$GPP[dat$site == 'GR' & dat$year == 2020], na.action = na.pass,
        xlim = c(0, 20), ylim = c(-0.2, 1), lwd = 2)
    par(new = T)
    acf(dat$GPP[dat$site == 'GR' & dat$year == 2021], na.action = na.pass,
        col = 2, xlim = c(-0.15, 19.85), ylim = c(-0.2, 1), xaxt = 'n', lwd = 2)
    mtext('GR', line = -2, adj = .9)
    acf(dat$GPP[dat$site == 'GC' & dat$year == 2020], na.action = na.pass,
        xlim = c(0, 20), ylim = c(-0.2, 1), lwd = 2)
    par(new = T)
    acf(dat$GPP[dat$site == 'GC' & dat$year == 2021], na.action = na.pass,
        col = 2, xlim = c(-0.15, 19.85), ylim = c(-0.2, 1), xaxt = 'n', lwd = 2)
    mtext('GC', line = -2, adj = .9)
    acf(dat$GPP[dat$site == 'BM' & dat$year == 2020], na.action = na.pass,
        xlim = c(0, 20), ylim = c(-0.2, 1), lwd = 2)
    par(new = T)
    acf(dat$GPP[dat$site == 'BM' & dat$year == 2021], na.action = na.pass,
        col = 2, xlim = c(-0.15, 19.85), ylim = c(-0.2, 1), xaxt = 'n', lwd = 2)
    mtext('BM', line = -2, adj = .9)
    acf(dat$GPP[dat$site == 'BN' & dat$year == 2020], na.action = na.pass,
        xlim = c(0, 20), ylim = c(-0.2, 1), lwd = 2)
    par(new = T)
    acf(dat$GPP[dat$site == 'BN' & dat$year == 2021], na.action = na.pass,
        col = 2, xlim = c(-0.15, 19.85), ylim = c(-0.2, 1), xaxt = 'n', lwd = 2)
    mtext('BN', line = -2, adj = .9)
dev.off()


png('figures/GPP_all_sites.png', width = 800, height = 600)
    ggplot(dat, aes(doy, GPP, col = factor(year))) +
        geom_line()+
        scale_color_manual(values = c(1,2)) +
        facet_wrap(.~site, ncol = 2) +
        theme_bw()
dev.off()
png('figures/GPP_light_byyear.png', width = 800, height = 600)
    ggplot(dat, aes(SW, GPP, col = factor(year))) +
        geom_point()+
        # scale_color_manual(values = c(1,2)) +
        facet_wrap(.~site, ncol = 2) +
        theme_bw()
dev.off()

png('figures/GPP_light_bybiomass.png', width = 800, height = 600)
BM <-     dat %>%
        group_by(site, year) %>%
        mutate(biomass_gm2 = na.approx(biomass_gm2, na.rm = F)) %>%
        filter(doy < 280,
               site == 'BM') %>%
    ggplot(aes(SW, GPP, col = biomass_gm2)) +
        geom_point()+
        scale_color_continuous(type = 'viridis') +
        # facet_wrap(.~site, ncol = 2) +
        theme_bw()
ggpubr::ggarrange(PL, DL, GR, GC, BM, BN, ncol = 2, nrow = 3, common.legend = T)
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
fit <- stan(file = 'code/model/stan_code/ar1_model_lognormal.stan',
            data = mod_dat,
            warmup = 2000, iter = 3000,
            chains = 4, cores = 4)

traceplot(fit, pars= c('b0', 'phi', 'b1', 'sigma_proc', 'sigma_obs'), ncol = 1)
plot(fit, pars= c('b0', 'phi', 'b1', 'sigma_proc', 'sigma_obs'))
pairs(fit, pars= c('b0', 'phi', 'b1', 'sigma_proc', 'sigma_obs'))

