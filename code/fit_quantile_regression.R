# Calculate slope of 90% quantile regression for each of the UCFR sites:

library(tidyverse)
library(lubridate)
# install.packages('lqmm')
library(lqmm)
# install.packages('qrLMM')
library(qrLMM)
library(lme4)

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
met <- mutate(met, year = factor(year)) %>%
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
met$AR = met$ER + met$GPP*met$ARf

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

dev.off()

# compare to biomass Growth ####

biogams <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv')
bio <- read_csv('data/biomass_data/biomass_working_data_summary.csv')
NPP <- bio %>%
    mutate(year = factor(lubridate::year(date))) %>%
    group_by(site, year) %>%
    summarize(fila_chla_meas = mean(filamentous_chla_mgm2),
              epil_chla_meas = mean(epilithon_chla_mgm2),
              fila_meas = mean(filamentous_gm2),
              epil_meas = mean(epilithon_gm2)) %>%
    left_join(NPP, by = c('site', 'year'))
NPP <- mutate(biogams, year = factor(year)) %>%
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
mod1 <- lm(NPP ~ 0 + epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod2 <- lm(NPP/light ~ 0 + epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod3 <- lm(NPP/light ~ 0 + epil_gm2 + fila_gm2, bm_met)

modmix <- lme4::lmer(NPP/light ~ 0 + epil_chla_mgm2 + fila_chla_mgm2 +
                       (0+epil_chla_mgm2|site) + (0+fila_chla_mgm2|site),
                     data = bm_met,
                     REML = FALSE,
                     control = lmerControl(optimizer ="Nelder_Mead"))
summary(mod2)
summary(modmix)
ranef(modmix)

bm_met %>%
    mutate(preds = predict(mod2, newdata = bm_met)) %>%
ggplot(aes(NPP, preds*light, col = site)) +
    geom_point() + geom_abline(intercept = 0, slope = 1)
bm_met %>%
    mutate(preds = predict(modmix, newdata = bm_met)) %>%
ggplot(aes(NPP, preds*light, col = site)) +
    geom_point() + geom_abline(intercept = 0, slope = 1)

ggplot(bm_met, aes(NPP, 0.435 * epil_gm2*light + 0.02394*fila_gm2*light,
                   col = site)) + geom_point() + geom_abline(intercept = 0, slope = 1)

coefs <- fixef(modmix)

bm_met <- bm_met %>%
    mutate(fila_prod_gCd = (coefs[2] * fila_chla_mgm2*light) * 14/32 ,
           fila_turnover = case_when(fila_prod_gCd > 0.04 ~ fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epil_prod_gCd = (coefs[1] * epil_chla_mgm2*light) * 14/32 ,
           epil_turnover = epil_gm2/2/epil_prod_gCd)


# bm_met <- bm_met %>%
#     mutate(fila_prod_gCd = (0.02394 * fila_gm2*light)*14/32 ,
#            fila_turnover = case_when(fila_prod_gCd > 0.04 ~ fila_gm2/2/fila_prod_gCd,
#                                      TRUE ~ NA_real_),
#            epil_prod_gCd = (0.435 * epil_gm2*light)*14/32 ,
#            epil_turnover = epil_gm2/2/epil_prod_gCd)

bm_fila <- filter(bm_met, fila_prod_gCd > 0.04)
cols <- c()

p1 <- bm_fila %>%
    pivot_longer(cols = ends_with('turnover'),
                 values_to = 'turnover', names_to = 'Biomass',
                 names_pattern = '(epil|fila)_.*') %>%
    ggplot(aes(x=turnover, group = Biomass, fill = Biomass)) +
    geom_density(adjust=1.5, alpha=.4) +
    scale_fill_manual(values = c('#1B9EC9', '#97BB43'))+
    xlim(0, 40)+
    ylab('Density')+
    xlab('Turnover time of biomass fraction (days)')+
    theme_classic()+
    theme(legend.position = 'none')
p1 <- p1 + annotate(geom = 'text', x = 10, y = 0.25,
              label="Epilithon", col = '#1B9EC9')
p1 <- p1 + annotate(geom = 'text', x = 27, y = 0.05,
              label="Filamentous", col = '#97BB43')
p2 <- ggplot(npp, aes(fila_gm/2/(biomass_C)*100, days_bio, col = epil_gm))+
        geom_point(size = 1.5)+
        xlab('Filamentous fraction of total biomass (%)')+
        ylab('Turnover time of standing stock (days)')+
        scale_color_gradient('Epilithon \nmass (g/m2)',
                             low = 'grey97', high = '#1B9EC9')+
        geom_point(size = 1.6, pch = 1, col = 'black')+
        theme_classic()+
        # guide_colorbar(frame.colour = 'black')+
        theme(legend.position = c(0.2, 0.7))
png('figures/turnover_time_by_composition.png',
    width = 6.5, height = 3.5, units = 'in', res = 300)
    ggpubr::ggarrange(p1, p2)
dev.off()


plot(density(bm_fila$fila_turnover), ylim = c(0, 0.25), xlim = c(0, 50),
     main = 'Algal Biomass Turnover Time',
     xlab = 'Days', ylab = 'Density')
par(new = T)
plot(density(bm_met$epil_turnover), xlim = c(0, 50),
     xlab = '', ylab = '', yaxt = 'n', main = '')
mtext('Biofilm', line = -5, adj = 0.12)
mtext('Filamentous', line = -20, adj = 0.56)

ggplot(bm_met, aes(date, fila_prod_gCd))+
    geom_line(col = 'forestgreen') +
    geom_line(aes(y = epil_prod_gCd), col = 'sienna')+
    facet_grid(site~year, scales = 'free_x')
png('figures/NPP_vs_standing_crop_all_sites.png',
    width = 6.5, height = 6.5,
    res = 300, units = 'in')
    bm_met %>%
        mutate(doy = as.numeric(format(date, '%j')),
               Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
               site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
        select(site, Date, year, fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
        pivot_longer(cols = starts_with(c('fila', 'epil')),
                     values_to = 'value', names_to = c('biomass', 'measure'),
                     names_pattern = '(fila|epil)_(.*)') %>%
        pivot_wider(names_from = measure, values_from = value) %>%
        mutate(gCm2 = gm2/2)%>%
    ggplot(aes(prod_gCd, gCm2, col = Date, pch = biomass))+
        geom_point() +
        facet_grid(site~year)+
        theme_classic()+
        scale_shape_manual('Biomass', values = c(1, 19))+
        ylab(expression(paste('Biomass standing crop (gC ', m^-2, ')')))+
        xlab(expression(paste('Biomass production (gC ', m^-2, d^-1, ')')))+
        theme(panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, units = 'in'))
dev.off()


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
