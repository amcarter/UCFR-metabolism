# Calculate slope of 90% quantile regression for each of the UCFR sites:

library(tidyverse)
library(lubridate)
# install.packages('lqmm')
library(lqmm)
# install.packages('qrLMM')
library(qrLMM)

source('code/metabolism/functions_examine_SM_output.R')

met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK.csv') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

met <- met%>%
    mutate(siteyear = paste(site, year, sep = '_')) %>%
    group_by(siteyear) %>%
    mutate(across(starts_with(c('GPP', 'ER')), zoo::na.approx, na.rm = F)) %>%
    filter(!is.na(GPP))

ggplot(met, aes(GPP, ER, col = factor(year))) +
    geom_point() +
    facet_wrap(.~site, scales = 'free') +
    geom_abline(slope = -1, intercept = 0) +
    ylim(-15,0) + xlim(0,20)

# fit quantile regression model:
y = met$ER
x = cbind(1, met$GPP)
z = cbind(1, met$GPP)
# qm90sy <- qrLMM::QRLMM(y,x,z, groups = met$siteyear, p = 0.9,
#              MaxIter = 500, M = 20)
#
# saveRDS(qm90sy, 'data/model_fits/quantile_PR_fits.rds')

qm90sy <- readRDS('data/model_fits/quantile_PR_fits.rds')
# qm90 <- lqmm(ER ~ GPP, random = ~GPP, group =site, data = met, tau = 0.9,
#              control = lqmmControl(LP_max_iter = 2000))
# qm90 <- lqmm(ER ~ GPP, random = ~GPP, group =site, data = met, tau = 0.9,
#              control = lqmmControl(LP_max_iter = 2000))
# qm90sy <- lqmm(ER ~ GPP, random = ~GPP, group =siteyear, data = met, tau = 0.9)

# summary(qm90)
# ranef.lqmm(qm90)

str(qm90sy)
beta <- qm90sy$res$beta
q90 <- data.frame(qm90sy$res$weights)
q90 <- q90 %>%
    mutate(siteyear = unique(met$siteyear),
           site = substr(siteyear, 1,2),
           year = as.numeric(substr(siteyear, 4,7)),
           intercept = X1 + beta[1,1],
           q90 = X2 + beta[2,1],
           ARf = -q90) %>%
    select(-siteyear, -X1, -X2)

plot(density(q90$ARf))

q90 <- mutate(q90, year = factor(year),
              site = factor(site,
                            levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))
met <- met %>%
    left_join(select(q90, site, year, ARf), by = c('site', 'year')) %>%
    mutate(NPP = GPP * (1-ARf))

png('figures/quantile_regression.png', width = 5, height = 3.5, units = 'in',
    res = 300)
    met %>%
        mutate(year = factor(year),
               site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
        ggplot(aes(GPP, ER, col = year)) +
        geom_point() +
        geom_abline(intercept = 0, slope = -1, lty = 2, col = 'grey40')+
        geom_abline(data = q90, aes(intercept = intercept, slope = q90, col = year))+
        facet_wrap(site~.)+
        ylim(-25,0)+
        ylab(expression(paste('ER (g ', O[2], m^-2, d^-1, ')'))) +
        xlab(expression(paste('GPP (g ', O[2], m^-2, d^-1, ')'))) +
        theme_classic()+
        theme(panel.border = element_rect(fill = NA))
dev.off()

NPP <- met %>%
    group_by(site, year) %>%
    summarize(NPP = mean(NPP), GPP = mean(GPP), ER = mean(ER)) %>%
    left_join(q90, by = c('site', 'year'))

ggplot(NPP, aes(ER, ARf, col = factor(year))) + geom_point(size = 2) +theme_minimal()

min(met$date)

met %>% filter(doy >=221) %>%
    group_by(siteyear) %>%
    mutate(year = factor(year),
           doy = doy - 220,
           cumNPP = cumsum(NPP) * 14/32) %>%
ggplot(aes(doy, cumNPP, col = year))+
    geom_line(linewidth = 1.2) +
    facet_wrap(.~site, scales = 'free') +
    ylab('Cumulative biomass production (g C)')+
    xlab('Days starting Aug 8')+
    theme_bw()

mm <- met %>% #filter(doy >=221) %>%
    group_by(siteyear) %>%
    mutate(year = factor(year),
           doy = doy - 167,
           cumNPP = cumsum(NPP) * 14/32) %>%
    select(site, date, year, doy, cumNPP) %>%
    ungroup() %>% group_by(site)

mm <- lapply(group_split(mm),
       function(x){
           # print(x$siteyear[1])
           day0 <- min(x$doy[x$year == 2020])
           val <- x$cumNPP[x$year == 2021 & x$doy == day0]
           x$cumNPP[x$year == 2020] <- x$cumNPP[x$year == 2020] + val

           return(x)
       }) %>% Reduce(bind_rows, .)

ggplot(mm, aes(doy, cumNPP, col = year)) +
    geom_line(linewidth = 1.2) +
    facet_wrap(.~site) +
    xlab('Days starting July 15') +
    ylab('Cumulative biomass production (g C)')+
    theme_bw()

met %>%
    ggplot(aes(date, ER + GPP * ARf, col = factor(year)))+
    geom_point() +
    geom_line(aes(y = ER)) +facet_grid(site~year, scales = 'free_x')
# compare to biomass Growth ####

biogams <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv')
bio <- read_csv('data/biomass_data/biomass_working_data_summary.csv')
NPP <- bio %>%
    mutate(year = lubridate::year(date)) %>%
    group_by(site, year) %>%
    summarize(fila_chla_meas = mean(filamentous_chla_mgm2),
              epil_chla_meas = mean(epilithon_chla_mgm2),
              fila_meas = mean(filamentous_gm2),
              epil_meas = mean(epilithon_gm2)) %>%
    left_join(NPP, by = c('site', 'year'))
NPP <- biogams %>%
    group_by(site, year) %>%
    summarize(fila_chla = mean(fila_chla_mgm2_fit),
              epil_chla = mean(epil_chla_mgm2_fit),
              fila_gm = mean(fila_gm2_fit),
              epil_gm = mean(epil_gm2_fit)) %>%
    left_join(NPP, by = c('site', 'year'))

ggplot(NPP, aes(epil_meas, epil_gm, col = site)) +
    geom_point(size = 2)
ggplot(NPP, aes(fila_meas, fila_gm, col = site)) +
    geom_point(size = 2)
ggplot(NPP, aes(epil_gm, ARf, col = site)) +
    geom_point(size = 2)
ggplot(NPP, aes(epil_chla/epil_gm, ARf, col = site)) +
    geom_point(size = 2)
ggplot(NPP, aes(fila_gm, ARf, col = site)) +
    geom_point(size = 2)

npp <- NPP %>%
    select(-ends_with('_meas')) %>%
    mutate(biomass_C = (epil_gm + fila_gm)/2,
           NPP_C = NPP * 14/32,
           days_epil = epil_gm/2/NPP_C,
           days_fila = fila_gm/2/NPP_C,
           days_bio = biomass_C/NPP_C)
ggplot(npp, aes(fila_gm/2/(biomass_C), days_bio, col = epil_gm))+
    geom_point(size = 2)+
    xlab('Filamentous fraction of total biomass')+
    ylab('Turnover time of biomass standing stock (d)')+
    theme_classic()
ggplot(npp, aes(NPP_C, fila_gm, col = factor(year)))+
    geom_point(size = 2)
summary(npp)

NPP$NPP * 30
mean(npp$days_bio)
plot(density(npp$days_bio))

# growth rate calculations ####
light <- read_csv('data/site_data/daily_modeled_light_all_sites.csv')
bm_met <- select(biogams, site, date, epil_gm2_fit, fila_gm2_fit,
       epil_chla_mgm2_fit, fila_chla_mgm2_fit) %>%
    rename_with(~gsub('_fit', '', .x)) %>%
    left_join(select(ungroup(met), site, date, year, GPP, ER, ARf, NPP),
              by = c('site', 'date')) %>%
    left_join(select(light, site, date, PAR_surface)) %>%
    mutate(light = PAR_surface/max(PAR_surface)) %>%
    filter(!is.na(GPP))

mod <- lm(GPP ~ 0 + epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod2 <- lm(NPP ~ 0 + epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod2 <- lm(NPP/light ~ 0 + epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod2 <- lm(NPP/light ~ 0 + epil_gm2 + fila_gm2, bm_met)
summary(mod2)

ggplot(bm_met, aes(NPP, 0.1719 * epil_chla_mgm2*light + 0.0154*fila_chla_mgm2*light,
                   col = site)) + geom_point() + geom_abline(intercept = 0, slope = 1)

ggplot(bm_met, aes(NPP, 0.435 * epil_gm2*light + 0.02394*fila_gm2*light,
                   col = site)) + geom_point() + geom_abline(intercept = 0, slope = 1)


bm_met <- bm_met %>%
    mutate(fila_prod_gCd = (0.01546 * fila_chla_mgm2*light)*14/32 ,
           fila_turnover = case_when(fila_prod_gCd > 0.04 ~ fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epil_prod_gCd = (0.1719 * epil_chla_mgm2*light)*14/32 ,
           epil_turnover = epil_gm2/2/epil_prod_gCd)
bm_met <- bm_met %>%
    mutate(fila_prod_gCd = (0.02394 * fila_gm2*light)*14/32 ,
           fila_turnover = case_when(fila_prod_gCd > 0.04 ~ fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epil_prod_gCd = (0.435 * epil_gm2*light)*14/32 ,
           epil_turnover = epil_gm2/2/epil_prod_gCd)

bm_fila <- filter(bm_met, fila_prod_gCd > 0.04)

plot(density(bm_fila$fila_turnover), ylim = c(0, 0.06), xlim = c(0, 80),
     main = 'Algal Biomass Turnover Time',
     xlab = 'Days', ylab = 'Density')
par(new = T)
plot(density(bm_met$epil_turnover), xlim = c(0, 80),
     xlab = '', ylab = '', yaxt = 'n', main = '')
mtext('Biofilm', line = -5, adj = 0.12)
mtext('Filamentous', line = -14, adj = 0.56)

ggplot(bm_met, aes(date, fila_prod_gCd))+
    geom_line(col = 'forestgreen') +
    geom_line(aes(y = epil_prod_gCd), col = 'sienna')+
    facet_grid(site~year, scales = 'free_x')

bm_met %>%
    mutate(doy = as.numeric(format(date, '%j')),
           site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
ggplot(aes(fila_prod_gCd, fila_gm2, col = doy))+
    geom_point() +
    facet_grid(site~year)+
    theme_classic()+
    theme(panel.border = element_rect(fill = NA))
bm_met %>%
    mutate(doy = as.numeric(format(date, '%j')),
           site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
ggplot(aes(epil_chla_mgm2, GPP, col = doy))+
    geom_point() +
    facet_grid(site~year)+
    theme_classic()+
    theme(panel.border = element_rect(fill = NA))


# Comparison to individual quantile regressions: ####

library(quantreg)
q90 <- data.frame(siteyear = unique(met$siteyear),
           slope = NA_real_,
           yintercept = NA_real_)

for(s in unique(met$siteyear)){
    d <- filter(met, siteyear == s)
    rqfit <- rq(ER ~ GPP, tau = 0.9, data = d)
    q90$slope[q90$siteyear == s] <- summary(rqfit)$coefficients[2,1]
    q90$yintercept[q90$siteyear == s] <- summary(rqfit)$coefficients[1,1]
}

q90 <- q90 %>%
    mutate(site = substr(siteyear, 1,2),
           year = as.numeric(substr(siteyear, 4,7)))

ggplot(met, aes(GPP, ER)) +
    geom_point() +
    geom_abline(data = q90, aes(intercept = yintercept,
                                slope = slope))
    facet_grid(site~year)
