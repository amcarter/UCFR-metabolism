# Calculate slope of 90% quantile regression for each of the UCFR sites:

library(tidyverse)
library(lubridate)
# install.packages('lqmm')
library(lqmm)
# install.packages('qrLMM')
library(qrLMM)
library(lme4)
library(viridis)

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
abline(v = -beta[2,1])

q90 <- mutate(q90, year = factor(year),
              site = factor(site,
                            levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

write_csv(q90, 'data/quantile_PR_fits_summary.csv')

met <- mutate(met, year = factor(year)) %>%
    left_join(select(q90, site, year, ARf), by = c('site', 'year')) %>%
    mutate(NPP = GPP * (1-ARf),
           NPP_globalARf = GPP * (1 + beta[2,1]),
           AR = GPP*(ARf))

png('figures/quantile_regression.png', width = 5, height = 3.5, units = 'in',
    res = 300)
    met %>%
        mutate(year = factor(year),
               site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
        ggplot(aes(GPP, ER, col = year)) +
        geom_point(size = 0.6) +
        geom_abline(intercept = beta[1,1], slope = beta[2,1], lty = 2, col = 'grey40')+
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
    summarize(NPP = mean(NPP), NPP2 = mean(NPP2), GPP = mean(GPP), ER = mean(ER)) %>%
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
# Maybe someday try making this into aplot that also shows net change in biomass
# through time. The difference should be the turnover?
# You could even add ER on here and then see what is left.
met %>% #filter(doy >=221) %>%
    group_by(siteyear) %>%
    mutate(year = factor(year),
           date = as.Date(paste0('2020-', doy), "%Y-%j"),
           cumNPP = cumsum(NPP) * 14/32) %>%
ggplot(aes(date, cumNPP))+
    geom_line(linewidth = 1.2) +
    facet_grid(site~year, scales = 'free') +
    ylab('Cumulative biomass production (g C)')+
    xlab('Date')+
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
    ggplot(aes(date, GPP*(1-ARf)))+
    geom_line(col = 'aquamarine3', size = 1.5) +
    geom_line(aes(y = GPP), size = 1.2) +
    facet_grid(site~year, scales = 'free_x')+
    theme_classic()+
    ylab('Respiration') +
    xlab('Date')+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'in'))
met %>%
    ggplot(aes(date, ER + GPP * ARf))+
    geom_line(col = 'aquamarine3', size = 1.5) +
    geom_line(aes(y = ER), size = 1.2) +
    facet_grid(site~year, scales = 'free_x')+
    theme_classic()+
    ylim(-22,0)+
    ylab('Respiration') +
    xlab('Date')+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'in'))

png('figures/poster/GPP_breakdown.png',
    width = 6.5, height = 1.8,
    res = 300, units = 'in')
met %>%
    filter(site == 'BN', year == 2021) %>%
    mutate(bottom = 0) %>%
    ggplot(aes(x = date)) +
    geom_ribbon(aes(ymin = AR, ymax = GPP),
                col = 'black', fill = alpha('aquamarine3', 0.5))+
    geom_ribbon(aes(ymin = bottom, ymax = AR), col = 'grey20', fill = 'grey70') +
    theme_classic()+
    ylab(expression(paste('GPP (g ', O[2], m^-2, d^-1,')'))) +
    xlab('Date')
dev.off()
png('figures/poster/GPP.png',
    width = 6, height = 2,
    res = 300, units = 'in')
met %>%
    filter(site == 'BN', year == 2021) %>%
    ggplot(aes(x = date, y = GPP)) +
    geom_line(col = 'forestgreen')+
    geom_point(col = 'forestgreen', size = 0.5)+
    geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), col = 'forestgreen') +
    theme_classic()+
    ylab('GPP') +
    xlab('Date')
dev.off()

png('figures/respiration_breakdown.png',
    width = 6.5, height = 6.5,
    res = 300, units = 'in')
# why is this broken for BM 2021?
met %>% #ungroup() %>%
    mutate(ER = AR,
           AR = 0,
           group = 1) %>%
    arrange(date, .by_group = TRUE)%>%
    bind_rows(arrange(met, desc(date), .by_group = TRUE)) %>%
    # arrange(site, year, group) %>%
    ggplot(aes(date, AR))+
    geom_polygon(fill = 'aquamarine3') +
    geom_polygon(aes(y = ER), fill = 'grey50') +
    facet_grid(site~year, scales = 'free_x')+
    theme_classic()+
    ylim(-22,0)+
    ylab('Respiration') +
    xlab('Date')+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'in'))
met$AR <- -met$AR
met %>% #ungroup() %>%
    mutate(GPP = AR,
           AR = 0,
           group = 1) %>%
    arrange(date, .by_group = TRUE)%>%
    bind_rows(arrange(met, desc(date), .by_group = TRUE)) %>%
    # arrange(site, year, group) %>%
    ggplot(aes(date, AR))+
    geom_polygon(fill = 'grey50') +
    geom_polygon(aes(y = GPP), fill = 'aquamarine3') +
    facet_grid(site~year, scales = 'free_x')+
    theme_classic()+
    # ylim(-22,0)+
    ylab('Productivity') +
    xlab('Date')+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'in'))

dev.off()

