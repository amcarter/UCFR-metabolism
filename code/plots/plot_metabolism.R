# Plot metabolism and model fits

library(tidyverse)
library(lubridate)
source('code/metabolism/functions_examine_SM_output.R')

met <- read_csv('data/metabolism_compiled_all_sites_2.csv') %>%
    mutate(year = factor(year(date)),
           doy = as.numeric(format(date, '%j')),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))
good_met <- met %>%
    filter(is.na(DO_fit) | DO_fit != 'bad') %>%
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

# write_csv(dat, 'data/prepared_data/compiled_prepared_data.csv')
dat <- read_csv('data/prepared_data/compiled_prepared_data.csv')

png('figures/metabolism_across_sites.png', width = 5, height = 6,
    res = 300, units = 'in')
    met %>%
        mutate(NEP = ER + GPP,
               date = as.Date(paste0('2021-', doy), format = '%Y-%j')) %>%
        filter(is.na(DO_fit) | DO_fit != 'bad') %>%
        ggplot(aes(date, GPP)) +
        geom_line(col = '#007828') +
        geom_point(col = '#007828', size = 0.7) +
        geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper),
                      col = '#007828') +
        # geom_point() +
        geom_line(aes(y = ER), col = '#A84F06') +
        geom_point(aes(y = ER), col = '#A84F06', size = 0.7) +
        geom_errorbar(aes(ymin = ER.lower, ymax = ER.upper),
                      col = '#A84F06') +
        # geom_line(aes(y = NEP), lty = 2) +
        geom_hline(yintercept = 0, size = 0.5, col = 'grey50')+
        # geom_point(aes(y = ER)) +
        facet_grid(site~year, scales = 'free_x') +
        theme_bw() +
        ylim(-24,24)+
        ylab(expression(paste('Metabolism (g ', O[2], m^-2, d^-1, ')'))) +
        xlab('Date')
dev.off()

png('figures/K600xER_across_sites_normal_pool.png', width = 6, height = 5,
    res = 300, units = 'in')
    met %>%
        # filter(is.na(DO_fit) | DO_fit != 'bad') %>%
        # filter(year == 2020) %>%
        ggplot(aes(doy, ER, col = year)) +
        geom_line(size = .8) +
        facet_wrap(.~site) +
        ylab(expression(paste('ER (g ', O[2], m^-2, d^-1, ')'))) +
        xlab(expression(paste('K600 (', d^-1, ')'))) +
        theme_bw()
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
mod <- lm(ER ~ doy, data = PL)
summary(mod)
