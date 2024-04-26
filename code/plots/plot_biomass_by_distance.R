
library(tidyverse)

# Plot biomass by distance downstream
mdat <- read_csv('data/site_data/site_data.csv') %>%
    select(site = sitecode, dist = 'Distance downstream (km)') %>%
    filter(!is.na(site))

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    filter(!is.na(site) & site != 'CR') %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           year = as.factor(lubridate::year(date)),
           site_year = paste(site, year, sep = '_'))
meas <- select(biomass, date, site, sample,
               epil_gm2_meas = epilithon_gm2,
               fila_gm2_meas = filamentous_gm2,
               epil_chla_mgm2_meas = epilithon_chla_mgm2,
               fila_chla_mgm2_meas = filamentous_chla_mgm2) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 names_to = c('biomass_type', 'units', 'stat'),
                 names_pattern = '([a-z]+)_([a-z0-9_]+)_([a-z]+)',
                 values_to = 'value') %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    mutate(year = lubridate::year(date))

qq <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv')
qq <- qq %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 names_to = c('biomass_type', 'units', 'stat'),
                 names_pattern = '([a-z]+)_([a-z0-9_]+)_([a-z]+)',
                 values_to = 'value') %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    group_by(site, biomass_type, units) %>%
    mutate(across(c('fit', 'se'), ~zoo::na.approx(., x = date, na.rm = F)))


png('figures/biomass_fits_by_distance.png', width = 10, height = 6, units = 'in',
    res = 300)
qq %>% left_join(mdat, by = 'site') %>%
    filter(units == 'gm2') %>%
    mutate(dist = round(dist))%>%
    ggplot(aes(factor(dist), fit, fill = biomass_type))+
    geom_boxplot(alpha = 0.6) +
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    # geom_point(data = filter(left_join(meas, mdat), units == 'gm2'),
    #            aes(y = meas, col = biomass_type),
    #            position = position_jitter(width = 0.2, height = 0))+
    labs(y = 'Biomass (g/m2)', x = 'Distance downstream (km)')+
    ylim(0,150)+
    facet_wrap(year~., ncol = 1, strip.position = 'right') +
    theme_classic()
dev.off()
png('figures/biomass_meas_by_distance.png', width = 10, height = 6, units = 'in',
    res = 300)
meas %>% left_join(mdat, by = 'site') %>%
    filter(units == 'gm2') %>%
    mutate(dist = round(dist))%>%
    ggplot(aes(factor(dist), meas, fill = biomass_type))+
    geom_boxplot(alpha = 0.6) +
    scale_fill_discrete(type = c('#1B9EC9', '#97BB43'))+
    labs(y = 'Biomass (g/m2)', x = 'Distance downstream (km)')+
    ylim(0,150)+
    facet_wrap(year~., ncol = 1, strip.position = 'right') +
    theme_classic()
dev.off()


png('figures/biomass_meas_by_distance.png', width = 10, height = 6, units = 'in',
    res = 300)

    meas %>% left_join(mdat, by = 'site') %>%
        mutate(site = case_when(site == 'BM' ~ 'BG',
                                site == 'GR' ~ 'GAR',
                                site == 'PL' ~ 'WS',
                                TRUE ~ site),
               site = factor(site, levels = c('WS', 'DL', 'GAR', 'GC', 'BG', 'BN')),
               biomass_type = case_when(biomass_type == 'epil' ~ 'Epilithic',
                                        biomass_type == 'fila' ~ 'Filamentous')) %>%
        filter(units == 'gm2') %>%
        mutate(dist = round(dist))%>%
        ggplot(aes(factor(site), meas, col = biomass_type))+
        geom_boxplot(size = 1) +
        scale_color_discrete("", type = c('#1B9EC9', 'forestgreen'))+
        labs(y = 'Biomass (g AFDM /m2)', x = 'Site')+
        ylim(0,150)+
        facet_wrap(year~., ncol = 1, strip.position = 'right') +
        theme_bw()+
        theme(legend.position = "top",
              legend.text=element_text(size=20),
              legend.title = element_text(size = 20),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.text.x = element_text(angle =0, vjust = 0.5, hjust=0.5),
              plot.title = element_text(hjust = 0.5, size = 20),
              axis.text.y = element_text(margin = margin(0,0,0,0)),
              strip.text.y = element_text(size = 20, color = "black", face = "bold"))

dev.off()
