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
met <- mutate(met, year = factor(year)) %>%
    left_join(select(q90, site, year, ARf), by = c('site', 'year')) %>%
    mutate(NPP = GPP * (1-ARf),
           NPP2 = GPP * (1 + beta[2,1]))

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
met$AR = met$GPP*(met$ARf)

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

# compare to biomass Growth ####

biogams <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv')
bio <- read_csv('data/biomass_data/biomass_working_data_summary.csv')
min_epil_gm2 <- min(bio$epilithon_gm2)
min_fila_gm2 <- filter(bio, filamentous_gm2>0) %>% select(filamentous_gm2) %>% min()

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
ggplot(NPP, aes(fila_chla, ARf, col = year)) +
    geom_point(size = 2)
ggplot(NPP, aes((fila_chla + epil_chla)/(fila_gm + epil_gm), ARf, col = year)) +
    geom_point(size = 2)
ggplot(NPP, aes((fila_chla)/(fila_gm), ARf, col = site)) +
    geom_point(size = 2)
ggplot(NPP, aes(fila_chla/(fila_chla + epil_chla), ARf, col = year)) +
    geom_point(size = 2)
ggplot(NPP, aes(epil_chla/(fila_chla + epil_chla), ARf, col = year)) +
    geom_point(size = 2)

npp <- NPP %>%
    select(-ends_with('_meas')) %>%
    mutate(biomass_C = (epil_gm + fila_gm)/2,
           NPP_C = NPP * 14/32,
           NPP_C2 = NPP2 * 14/32,
           days_epil = epil_gm/2/NPP_C,
           days_fila = fila_gm/2/NPP_C,
           days_bio = biomass_C/NPP_C) %>%
    ungroup() %>%
    mutate(fila_bloom = c(0,1,0,1,0,0,1,1,1,0,0,0),
           fila_bloom = case_when(fila_bloom == 0 ~ 'No Bloom',
                                  fila_bloom == 1 ~ 'Bloom'))

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
           fila_turnover = case_when(fila_gm2 > min_fila_gm2 ~ fila_gm2/2/fila_prod_gCd,
                                     TRUE ~ NA_real_),
           epil_prod_gCd = (coefs[1] * epil_chla_mgm2*light) * 14/32 ,
           epil_turnover = epil_gm2/2/epil_prod_gCd)
bm_met %>%
    filter(year == 2021) %>%
    select(-year, -ARf, -PAR_surface) %>%
    write_csv('data/biomass_data/2021_turnovers_for_Rafa.csv')

# bm_met <- bm_met %>%
#     mutate(fila_prod_gCd = (0.02394 * fila_gm2*light)*14/32 ,
#            fila_turnover = case_when(fila_prod_gCd > 0.04 ~ fila_gm2/2/fila_prod_gCd,
#                                      TRUE ~ NA_real_),
#            epil_prod_gCd = (0.435 * epil_gm2*light)*14/32 ,
#            epil_turnover = epil_gm2/2/epil_prod_gCd)

bm_fila <- filter(bm_met, fila_gm2 > min_fila_gm2)

p1 <- bm_fila %>%
    pivot_longer(cols = ends_with('turnover'),
                 values_to = 'turnover', names_to = 'Biomass',
                 names_pattern = '(epil|fila)_.*') %>%
    ggplot(aes(x=turnover, group = Biomass, fill = Biomass)) +
    geom_density(adjust=1.5, alpha=.4) +
    scale_fill_manual(values = c('#1B9EC9', '#97BB43'))+
    xlim(0, 40)+
    ylab('Density')+
    xlab('Turnover time (d)')+
    theme_classic()+
    theme(legend.position = 'none',
          panel.border = element_rect(fill = NA))
p1 <- p1 + annotate(geom = 'text', x = 12, y = 0.25,
              label="Epilithic", col = '#1B9EC9')
p1 <- p1 + annotate(geom = 'text', x = 32, y = 0.06,
              label="Filamentous", col = '#97BB43')
p2 <- ggplot(npp, aes(fila_gm/2/(biomass_C)*100, days_bio))+
        xlab('Filamentous fraction of \ntotal biomass (%)')+
        ylab('Turnover time (d)')+
        geom_point(aes(pch = factor(fila_bloom)), size = 1.6, col = 'black')+
        theme_classic()+
        theme(legend.title = element_blank(),
              legend.position = c(0.22, 0.85),
              panel.border = element_rect(fill = NA))
png('figures/turnover_time_by_composition.png',
    width = 6.5, height = 3.5, units = 'in', res = 300)
    ggpubr::ggarrange(p1, p2)
dev.off()

ggplot(bm_met, aes(date, fila_prod_gCd))+
    geom_line(col = 'forestgreen') +
    geom_line(aes(y = epil_prod_gCd), col = 'sienna')+
    facet_grid(site~year, scales = 'free_x')
# png('figures/NPP_vs_standing_crop_all_sites.png',
#     width = 6.5, height = 6.5,
#     res = 300, units = 'in')
#     bm_met %>%
#         mutate(doy = as.numeric(format(date, '%j')),
#                Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
#                site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
#         select(site, Date, year, fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
#         pivot_longer(cols = starts_with(c('fila', 'epil')),
#                      values_to = 'value', names_to = c('biomass', 'measure'),
#                      names_pattern = '(fila|epil)_(.*)') %>%
#         pivot_wider(names_from = measure, values_from = value) %>%
#         mutate(gCm2 = gm2/2)%>%
#     ggplot(aes(prod_gCd, gCm2, col = Date, pch = biomass))+
#         geom_point() +
#         facet_grid(site~year)+
#         theme_classic()+
#         scale_shape_manual('Biomass', values = c(1, 19))+
#         ylab(expression(paste('Biomass standing crop (gC ', m^-2, ')')))+
#         xlab(expression(paste('Biomass production (gC ', m^-2, d^-1, ')')))+
#         theme(panel.border = element_rect(fill = NA),
#               panel.spacing = unit(0, units = 'in'))
# dev.off()
# png('figures/NPP_vs_standing_crop_all_sites.png',
#     width = 6.5, height = 6.5,
#     res = 300, units = 'in')
    coeff = 125
    p <- bm_met %>%
        mutate(doy = as.numeric(format(date, '%j')),
               Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
               site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
        mutate(fila_gm2 = case_when(fila_gm2 < min_fila_gm2 ~ NA_real_,
                                    TRUE ~ fila_gm2)) %>%
        select(site, Date, year, light,
               fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
        pivot_longer(cols = starts_with(c('fila', 'epil')),
                     values_to = 'value', names_to = c('biomass', 'measure'),
                     names_pattern = '(fila|epil)_(.*)') %>%
        pivot_wider(names_from = measure, values_from = value) %>%
        mutate(Light = factor(rep(" ", 2*nrow(bm_met)), levels = c("1"," ")),
               gCm2 = gm2/2,
               biomass = case_when(biomass == 'fila' ~ 'Filamentous',
                                   biomass == 'epil'~'Epilithic'),
               biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%

    ggplot(aes(Date, prod_gCd/gCm2, col = biomass))+
        geom_area(aes(Date, gCm2/coeff, fill = biomass), color = NA, alpha = 0.4)+
        geom_line(size = 1.2)+
        geom_line(aes(y = light/2.5, lty = Light), col = 'grey20') +
        facet_grid(site~year, scales = 'free_x', ) +
        scale_y_continuous(
            name = expression(paste('Production rate (', d^-1, ')')),
            sec.axis = sec_axis(~.*coeff,
                                name = expression(paste('Biomass (g C', m^-2, ')')))
        )+
        scale_fill_manual('Biomass standing stock',
                          values = c('#97BB43', '#1B9EC9')) +
        scale_color_manual('Biomass production rate',
                           values = c('#97BB43', '#1B9EC9')) +
        scale_linetype_manual('Light', values = c(2)) +
        theme_classic()+
        xlab('Date')+
        theme(strip.text.y = element_blank(),
              panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, units = 'in'),
              # legend.spacing.x = unit(1, 'cm'),
              # legend.justification = 'top')
              # legend.box.margin = margin(0, 5, 0, 5, "cm"),
              legend.position = 'top')+
        guides(color = guide_legend(title.position = "top", title.hjust = 0,
                                      order = 1),
               fill = guide_legend(title.position = "top", title.hjust = 0,
                                   order = 2),
               linetype = guide_legend(title.position = 'top', title.hjust = 0,
                                  order = 3))

    ann_text <- data.frame(Date = rep(as.Date('2020-07-16'), 6),
                           prod_gCd = rep(0.48, 6),
                           gCm2 = rep(1, 6),
                           year = rep(2020,6),
                           biomass = factor(rep('Epilithic', 6),
                                            levels = c('Filamentous','Epilithic')),
                           site = factor(c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'),
                                        levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

    p3 <- p + geom_text(data = ann_text, aes(label = site), col = 'black')

   p4 <- bm_met %>%
        mutate(site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
        group_by(site, year) %>%
        # filter(date > as.Date('2020-08-07'),
        #        date < as.Date('2020-10-19') | date > as.Date('2021-07-16'),
        #        date < as.Date('2021-09-27'))%>%
       # summarize(start = min(date),
       #           end = max(date),
       #           n = n())
        summarize(n = n(),
                  fila_Biomass = max(fila_gm2)/2,
                  epil_Biomass = max(epil_gm2)/2,
                  fila_cumprod = sum(fila_prod_gCd)/n*100,
                  epil_cumprod = sum(epil_prod_gCd)/n*100) %>%
       mutate(bloom = factor(case_when(fila_Biomass > 30 ~ 'Bloom',
                                       TRUE ~ 'No Bloom'))) %>%
        pivot_longer(cols = starts_with(c('epil', 'fila')),
                     values_to = 'value',
                     names_to = c('biomass', 'measure'),
                     names_pattern = '(fila|epil)_(Biomass|cumprod)') %>%
        mutate(measure = case_when(measure == 'cumprod' ~ 'Cumulative \nProduction',
                                   measure == 'Biomass' ~ 'Maximum \nBiomass'),
               measure = factor(measure, levels = c('Maximum \nBiomass',
                                                    'Cumulative \nProduction')),
               biomass = case_when(biomass == 'epil' ~ 'Epilithic',
                                   biomass == 'fila' ~ 'Filamentous'),
               biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
        ggplot(aes(x = measure, y = value, fill = biomass)) +
        geom_boxplot(outlier.shape=NA, alpha = 0.4)+
        geom_point(aes(pch = bloom, group = biomass),
                   position = position_jitterdodge(jitter.width = 0.25))+
        # geom_jitter(color="black", alpha=0.9) +
        scale_fill_manual('',
                          values = c('#97BB43', '#1B9EC9')) +
        ylab(expression(paste('Mass (g C ', m^-2, ')')))+
        xlab('')+
        theme_classic()+
        theme(legend.title = element_blank(),
              panel.border = element_rect(fill = NA),
              legend.spacing.y = unit(-0.1, "cm"),
              legend.position = c(0.24, 0.79),
              legend.background = element_rect(fill = NA))
png('figures/biomass_cumulative_and_turnover.png', width = 9, height = 3.5,
    units = 'in', res = 300)
    ggpubr::ggarrange(p4, p1,p2, ncol = 3,
                      labels = c('A', 'B', 'C'),
                      align = 'h')

dev.off()

    p5 <- ggplot(data.frame(a = 1, b = 1), aes(a,b)) +
        geom_point(col = 'white') +
        xlab('') + ylab('')+
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank())

png('figures/biomass_prod_and_turnover_3panel.png', width = 8, height = 8,
    units = 'in', res = 300)
    ggpubr::ggarrange(p3,
                      ggpubr::ggarrange(p2, p1, p4, p5, ncol = 1,
                                        labels = c('B', 'C', 'D', ''),
                                        heights = c(5.1, 5, 5.2, 1.3),
                                        align = 'v'),
                      ncol = 2, labels = 'A', #label.y = 0.917,
                      widths = c(2,1))
dev.off()
png('figures/biomass_prod_and_turnover.png', width = 6, height =8,
    units = 'in', res = 300)
    p3
dev.off()
png('figures/biomass_turnover_frac.png', width = 3.5, height = 3,
    units = 'in', res = 300)
p2
dev.off()


p2 <- bm_met %>%
    filter(year == 2021, site == 'GC') %>%

    mutate(doy = as.numeric(format(date, '%j')),
           Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
           site = factor(site, levels = c('PL', 'DL', 'GR','GC','BM','BN'))) %>%
    mutate(fila_gm2 = case_when(fila_gm2 < min_fila_gm2 ~ NA_real_,
                                TRUE ~ fila_gm2)) %>%
    select(site, Date, year, light,
           fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
    pivot_longer(cols = starts_with(c('fila', 'epil')),
                 values_to = 'value', names_to = c('biomass', 'measure'),
                 names_pattern = '(fila|epil)_(.*)') %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(gCm2 = gm2/2,
           biomass = case_when(biomass == 'fila' ~ 'Filamentous',
                               biomass == 'epil'~'Epilithic'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
    filter(biomass == 'Epilithic') %>%
    ggplot(aes(Date, gCm2, fill = biomass))+
    geom_area(color = NA, alpha = 0.7)+
    ylab(expression(paste('Biomass (g C', m^-2, ')')))+
    scale_fill_manual('Biomass standing stock',
                      values = c( '#1B9EC9')) +
    theme_classic()+
    xlab('Date')+
    ylim(0,47)+
    theme(legend.position = 'none')
# p1 <-
    bm_met %>%
    filter(year == 2021, site == 'GC') %>%
    mutate(doy = as.numeric(format(date, '%j')),
           Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
           fila_cum = cumsum(fila_prod_gCd),
           epil_cum = cumsum(epil_prod_gCd)) %>%
    ggplot(aes(Date, fila_cum))+
    geom_line(color = '#97BB43', size = 1.5)+
    geom_line(aes(y = epil_cum), color = '#1B9EC9', size = 1.5)+
    ylab(expression(paste('Biomass production (g C', m^-2, d^-1, ')')))+
    scale_color_manual('Biomass production rate',
                       values = c('#97BB43', '#1B9EC9')) +
    theme_classic()+
    xlab('Date')




bm_met %>%
    group_by(site, year) %>%
    mutate(fila_rate = fila_prod_gCd/(fila_gm2/2),
           epil_rate = epil_prod_gCd/(epil_gm2/2)) %>%
    filter(fila_gm2 >= min_fila_gm2) %>%
               summary()
bm_met %>%
    group_by(site, year) %>%
    filter((year == 2020 & site %in% c('GR', 'GC'))|
               (year == 2021 & site %in% c('GC', 'BN','BM')))%>%
    summarize(n = n(),
              fila_Biomass = max(fila_gm2)/2,
              epil_Biomass = max(epil_gm2)/2,
              fila_prod = mean(fila_prod_gCd),
              epil_prod = mean(epil_prod_gCd),
              fila_cumprod = sum(fila_prod_gCd),
              epil_cumprod = sum(epil_prod_gCd)) %>%
    mutate(n_epil = epil_cumprod/epil_Biomass,
           n_fila = fila_cumprod/fila_Biomass)%>%
    summary()

bm_met %>%
    filter(!(site %in% c('DL', 'PL') ))%>%
    ggplot(aes(date, fila_turnover)) +
    geom_line() +
    facet_grid(site~year, scale = 'free_x') + ylim(0,100)

    ggplot(aes(Date, NPP/(fila_gm2 + epil_gm2),
               col = fila_gm2/(epil_gm2 + fila_gm2)))+
        geom_point() +
        facet_grid(site~year, scales = 'free')+
        geom_line(aes(y = light))+
        theme_classic()+
        ylab(expression(paste('Biomass standing crop (gC ', m^-2, ')')))+
        xlab(expression(paste('Biomass production (gC ', m^-2, d^-1, ')')))+
        theme(panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, units = 'in'))
dev.off()


bm_met %>%
        mutate(fila_gm2 = case_when(fila_gm2 < min_fila_gm2 ~ NA_real_,
                                    TRUE ~ fila_gm2)) %>%
    mutate(fila_mgchlagC = fila_chla_mgm2/(fila_gm2/2),
           epil_mgchlagC = epil_chla_mgm2/(epil_gm2/2)) %>%
    ggplot(aes(date, fila_mgchlagC))+
    geom_point() +
    geom_point(aes(y = epil_mgchlagC), col = 2) +
    facet_grid(site~year, scales = 'free_x')
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
