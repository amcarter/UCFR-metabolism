# Calculate slope of 90% quantile regression for each of the UCFR sites:

library(tidyverse)
library(lubridate)
# install.packages('lqmm')
library(lqmm)
# install.packages('qrLMM')
library(qrLMM)
library(lme4)
library(viridis)
library(brms)
source('code/metabolism/functions_examine_SM_output.R')

met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK_correctedSE.csv') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

met <- met%>%
    mutate(siteyear = paste(site, year, sep = '_')) %>%
    group_by(siteyear) %>%
    mutate(GPP.se = (GPP.upper - GPP.lower)/(2*1.96),
           ER.se = (ER.upper - ER.lower)/(2*1.96))%>%
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

# Try running a bayesian quantile regression model:
library(brms)

qform <- bf(ER ~ GPP + (1 + GPP| siteyear), quantile = 0.9)
get_prior(qform, met, family = asym_laplace())


priors <- c(
    prior(normal(-0.5, 0.1), class = "b"),
    # prior(normal(0.5, 0.1), class = "meanme"),
    prior(normal(0, 0.02), class = "sd"),
    prior(normal(0, 1), class = "sd", group = "siteyear", coef = "Intercept")
)


fit_bqr <-
    brm(qform,
    # stancode(qform,
    family = asym_laplace(),
    prior = priors,
    data = met,
    iter = 6000,
    chains = 4, cores = 4,
    control = list(max_treedepth = 13))

saveRDS(fit_bqr, "data/model_fits/brms_quantile_regression_ARf.rds")
summary(fit_bqr)
ranef(fit_bqr)
plot(fit_bqr)
ests <- fitted(fit_bqr, dpar = "mu")
met$fit <- ests[,1]
met$fit.low <- ests[,3]
met$fit.high <- ests[,4]

brms::pp_check(fit_bqr, type = "error_scatter_avg_grouped", group = "siteyear")
par(mfrow = c(1,2))
pp1 <- brms::pp_check(fit_bqr) +
    ggtitle("Posterior Predictions of Data Distribution")
q90 <- function(y) quantile(y, 0.9)
pp2 <- brms::pp_check(fit_bqr, type = "stat", group = 'siteyear', stat = "q90") +
    ggtitle("Posterior Predictions of 90th Quantile")

png("figures/SI/ppred_checks_brms_quantile_reg.png", width = 8, height = 4.5,
    units = "in", res = 300)
    ggpubr::ggarrange(pp1, pp2)
dev.off()

ggplot(met, aes(GPP, ER, col = factor(year))) +
    geom_point() +
    geom_line(aes(y = fit)) +
    geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = factor(year)),
                col = NA, alpha = 0.3) +
    facet_wrap(.~site)

coef_df <- coef(fit_bqr)$siteyear
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
q90b <- q90 %>%
    mutate(siteyear = unique(met$siteyear),
           site = substr(siteyear, 1,2),
           year = as.numeric(substr(siteyear, 4,7)),
           intercept = X1 + beta[1,1],
           q90 = X2 + beta[2,1],
           ARf = -q90) %>%
    select(-X1, -X2)

q90 <- data.frame(coef_df[,,2]) %>%
    mutate(siteyear = rownames(coef_df[,,2])) %>%
    rename(slope = Estimate, slope.se = Est.Error, slope.lower = Q2.5, slope.upper = Q97.5) %>%
    left_join(q90b, by = "siteyear") %>%
    select(siteyear, site, year, q90, starts_with('slope')) %>%
    mutate(across(any_of(starts_with('slope')), \(x) abs(x)))

plot(density(q90$slope))
abline(v = -beta[2,1])

q90 <- mutate(q90, year = factor(year),
              site = factor(site,
                            levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

write_csv(q90, 'data/quantile_PR_fits_summary_brms.csv')

met <- mutate(met, year = factor(year)) %>%
    left_join(select(q90, site, year, starts_with('slope')), by = c('site', 'year')) %>%
    mutate(NPP = GPP * (1-slope),
           NPP.se = sqrt(GPP.se^2 + slope.se^2),
           NPP_globalARf = GPP * (1 + beta[2,1]),
           AR = GPP*(slope))

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

png('figures/quantile_regression_brms.png', width = 5, height = 3.5, units = 'in',
    res = 300)
    met %>%
        mutate(Year = factor(year),
               site = case_when(site == 'BM' ~ 'BG',
                                       TRUE ~ site),
               site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN'))) %>%
    ggplot( aes(GPP, ER, col = Year)) +
        geom_point(size = 0.3) +
        geom_abline(intercept = -3.52, slope = -0.51, lty = 2, col = 'grey40')+
        geom_point(size = 0.6) +
        geom_line(aes(y = fit)) +
        geom_ribbon(aes(ymin = fit.low, ymax = fit.high, fill = Year),
                    col = NA, alpha = 0.5) +
        facet_wrap(site~.)+
        ylim(-25,0)+
        ylab(expression(paste('ER (g ', O[2], m^-2, d^-1, ')'))) +
        xlab(expression(paste('GPP (g ', O[2], m^-2, d^-1, ')'))) +
        theme_classic()+
        theme(panel.border = element_rect(fill = NA))
dev.off()

NPP <- met %>%
    group_by(site, year) %>%
    summarize(NPP = mean(NPP),
              GPP = mean(GPP), ER = mean(ER)) %>%
    left_join(q90, by = c('site', 'year'))

ggplot(NPP, aes(ER, slope, col = factor(year))) + geom_point(size = 2) +theme_minimal()

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
    ggplot(aes(date, GPP*(1-slope)))+
    geom_line(col = 'aquamarine3', size = 1.5) +
    geom_line(aes(y = GPP), size = 1.2) +
    facet_grid(site~year, scales = 'free_x')+
    theme_classic()+
    ylab('Respiration') +
    xlab('Date')+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, units = 'in'))
met %>%
    ggplot(aes(date, ER + GPP * slope))+
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

