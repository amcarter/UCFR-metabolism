# Plot metabolism and model fits

library(tidyverse)
library(lubridate)
source('code/metabolism/functions_examine_SM_output.R')

met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK_correctedSE.csv') %>%
# met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK.csv') %>%
# metabolism_compiled_all_sites_2000iter_bdr_kss05.csv') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')))

good_met <- met %>%
    select(-errors) %>%
    mutate(across(starts_with(c('GPP', 'ER', 'K600')),
                  ~case_when((!is.na(DO_fit) & DO_fit == 'bad') ~ NA_real_,
                             TRUE ~ .))) %>%
    select(-DO_fit)

# compile the input data:
# dat <- data.frame()
# for(site in unique(met$site)){
#     d <- read_csv(paste0('data/prepared_data/', site, '2020_2021.csv'))
#     d$site = site
#     dat <- bind_rows(dat, d)
# }
#
# dat <- dat %>%
#     mutate(date = as.Date(solar.time)) %>%
#     group_by(site, date) %>%
#     summarize(across(-solar.time, mean, na.rm = T))
#
# write_csv(dat, 'data/prepared_data/compiled_prepared_data.csv')
dat <- read_csv('data/prepared_data/compiled_prepared_data.csv')%>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')))
bmet <- read_csv('data/model_fits/biomass_metab_model_data.csv') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')))

    bmet %>%
        mutate(NEP = ER + GPP,
               doy = as.numeric(format(date, '%j')),
               bmass = fila_chla_mgm2_fit+ epil_chla_mgm2_fit) %>%
        ggplot(aes(GPP, ER, col = site)) +
        geom_point() +
        # facet_grid(site~year) +
        # facet_grid(site~year, scales = 'free') +
        theme_bw() +
        # scale_color_viridis(option = 'D') +
        geom_abline(slope = -1, intercept = 0, col = 'grey50') +
        ylab(expression(paste('ER (g ', O[2], m^-2, d^-1, ')')))
    bmet %>%
        mutate(NEP = ER + GPP,
               doy = as.numeric(format(date, '%j')),
               bmass = fila_chla_mgm2_fit+ epil_chla_mgm2_fit,
               site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
        ggplot(aes(log(light), log(GPP), col = bmass)) +
        geom_point() +
        # facet_grid(site~year) +
        facet_grid(site~year, scales = 'free') +
        theme_bw() +
        # scale_color_viridis(option = 'D') +
        # geom_abline(slope = -1, intercept = 0, col = 'grey50')
        ylab(expression(paste('ER (g ', O[2], m^-2, d^-1, ')')))


lme4::lmer(log(GPP) ~log(light) + (1|site), bmet) %>% summary()

met <- met %>%
    mutate(NEP = GPP + ER,
           GPP.se = (GPP.upper - GPP.lower)/(2*1.96),
           ER.se = (ER.upper - ER.lower)/(2*1.96),
           NEP.se = sqrt(GPP.se^2 + ER.se^2 + GPP.se*ER.se),
           PR = -GPP/ER,
           PR.se = PR * sqrt((GPP.se/GPP)^2 + (ER.se/ER)^2 -
                                 2*(GPP.se * ER.se)/(GPP * ER)),
           month = month(date))
M_sum <- met %>%
    group_by(site, year) %>%
    summarize(PR = mean(PR, na.rm = T),
              PR.sd = sd(GPP/ER, na.rm = T),
              GPP.sd = sd(GPP, na.rm = T),
              ER.sd = sd(ER, na.rm = T),
              NEP.sd = sd(NEP, na.rm = T),
              GPP = mean(GPP, na.rm = T),
              ER = mean(ER, na.rm = T),
              NEP = GPP + ER) %>%
    # filter(site == 'GC', month %in% c(7,8,9)) %>% select(site, date, GPP, ER, NEP)
    mutate(trophic = case_when(NEP > 0 ~ 'A',
                               TRUE ~ ' '),
           PR = as.character(round(PR, 2)),
           PR.sd = as.character(round(PR.sd, 2)),
           PR_lab = (paste("PR = ", PR , "\U00B1" , PR.sd)))%>%
    select(-NEP) %>% ungroup() %>%
    mutate(date = case_when(year == 2020 ~ as.Date('2020-08-08'),
                            year == 2021 ~ as.Date('2021-09-17')),
           GPP = rep(-21, 12))


met %>% group_by(site, year) %>%
    summarize(PR_cor = cor(GPP, ER, use = "complete.obs"),
              NEP = mean(GPP + ER, na.rm = T)) %>%
    arrange(PR_cor)

# write_csv(M_sum, 'data/metabolism/auto_sites_figure_labels.csv')
#
png('figures/metabolism_across_sites.png', width = 6.5, height = 6.5,
    res = 300, units = 'in')
    M <- met %>%
        mutate(NEP = ER + GPP) %>%
        filter(is.na(DO_fit) | DO_fit != 'bad') %>%
        ggplot(aes(date, GPP, group = c(site))) +
        geom_hline(yintercept = 0, size = 0.5, col = 'grey75')+
        geom_line(col = 'grey30') +
        geom_ribbon(aes(ymin = GPP.lower, ymax = GPP.upper),
                      fill = '#007828', alpha = 0.6) +
        geom_line(aes(y = ER), col = 'grey30') +
        geom_line(aes(y = NEP), col = 'grey30') +
        geom_ribbon(aes(ymin = ER.lower, ymax = ER.upper),
                      fill = '#A84F06', alpha = 0.6 ) +
        geom_ribbon(aes(ymin = NEP - 1.96*NEP.se, ymax = NEP + 1.96*NEP.se),
                      fill = 'grey', alpha = 0.6 ) +
        facet_grid(site~year, scales = 'free_x') +
        theme_classic() +
        theme(panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, 'line'))+
        ylim(-22.5,22)+
        ylab(expression(paste('Metabolism (g ', O[2], m^-2, d^-1, ')'))) +
        xlab('Date')

    # M + geom_text(data = M_sum, aes(label = trophic), col = 'black')

    M + geom_text(data = M_sum, aes(label = PR_lab),
                   col = 'black', size = 3.2)
dev.off()


png('figures/NEP_all_sites.png', width = 4, height = 4,
    res = 300, units = 'in')
    met %>%
        mutate(NEP = ER + GPP) %>%
        mutate(date = case_when(year == 2021 ~ date - 365,
                                TRUE ~ date)) %>%
        filter(is.na(DO_fit) | DO_fit != 'bad') %>%
        ggplot(aes(date, NEP, group = c(site))) +
        geom_hline(yintercept = 0, size = 0.5, col = 'grey75')+
        geom_line(aes(y = NEP), col = 'black') +
        facet_grid(year~., scales = 'free_x') +
        theme_classic() +
        theme(panel.border = element_rect(fill = NA),
              panel.spacing = unit(0, 'line'))+
        # ylim(-24,24)+
        ylab(expression(paste('NEP (g ', O[2], m^-2, d^-1, ')'))) +
        xlab('Date')
dev.off()

png('figures/K600_plots_across_sites_normal_pool.png', width = 5.5, height = 7,
    res = 300, units = 'in')
    met %>%
        left_join(dat, by = c('site', 'date')) %>%
        mutate(Q = log10(discharge),
               DO_fit = if_else(DO_fit == 'bad', 'poor DO fit', NA_character_))%>%
        pivot_longer(cols = c('GPP', 'ER', 'Q'),
                     names_to = 'variable', values_to = 'value') %>%
        filter(is.na(DO_fit))%>%
        ggplot(aes(value, K600)) +
        geom_point(size = 1) +
        facet_grid(site~variable, scales = 'free_x') +
        ylab(expression(paste('K600 (', d^-1, ')'))) +
        xlab(expression(paste('       ER (g ', O[2], m^-2, d^-1, ')         GPP (g ', O[2], m^-2, d^-1, ')      Discharge (', m^3, s^-1,', lo',g[10],')'))) +
        theme_bw()+
        theme(legend.position = 'none',
              # legend.title = element_blank(),
              strip.background.x = element_blank(),
              strip.text.x = element_blank())
dev.off()

png('figures/K600xGPP_across_sites.png', width = 6, height = 5,
    res = 300, units = 'in')
    met %>%
        # filter(is.na(DO_fit) | DO_fit != 'bad') %>%
        # filter(year == 2020) %>%
        ggplot(aes(K600, GPP)) +
        geom_point(size = .8) +
        facet_wrap(.~site) +
        ylab(expression(paste('GPP (g ', O[2], m^-2, d^-1, ')'))) +
        xlab(expression(paste('K600 (', d^-1, ')'))) +
        theme_bw()
dev.off()

PL <- good_met %>% filter(site == 'PL') %>% plot_KxQ(dat = dat[dat$site == 'PL',] )
DL <- good_met %>% filter(site == 'DL') %>% plot_KxQ(dat = dat[dat$site == 'DL',] )
GR <- good_met %>% filter(site == 'GR') %>% plot_KxQ(dat = dat[dat$site == 'GR',] )
GC <- good_met %>% filter(site == 'GC') %>% plot_KxQ(dat = dat[dat$site == 'GC',] )
BM <- good_met %>% filter(site == 'BM') %>% plot_KxQ(dat = dat[dat$site == 'BM',] )
BN <- good_met %>% filter(site == 'BN') %>% plot_KxQ(dat = dat[dat$site == 'BN',] )

ggpubr::ggarrange(PL, DL, GR, GC, BM, BN, ncol = 3,nrow = 2, align = 'hv')

PL <- filter(met, site == 'PL')
mod <- lm(ER ~ K600, data = PL)
summary(mod)
