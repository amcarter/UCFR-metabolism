# Compare across metabolism fits with varying K600 sigma sigma

library(tidyverse)
source('code/UCFR_helpers.R')
source('code/metabolism/functions_examine_SM_output.R')

files <- grep('bdr\\.rds$', list.files('data/metabolism/metab_fits'), value = TRUE)

COMP_MET <- data.frame()

for(site in c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')){
    comp_met <- data.frame()

    sitefiles <- grep(paste0('^', site), files, value = TRUE)

    for(f in sitefiles){
        fit <- readRDS(paste0('data/metabolism/metab_fits/', f))
        K600_sigsig <- paste0('0.', str_match(f, '_([0-9]+)_bdr\\.rds$')[,2])
        K_pool <- str_match(f, '^[A-Z]{2}_([a-z]+)_')[2]
        if(K600_sigsig == '0.0'){
                met <- fit@fit %>%
                    select(date, GPP = GPP.daily, ER = ER.daily,
                           K600 = K600.daily) %>%
                    mutate(site = site, DO_fit = 'good')
        } else{
            bad_days <- unique(as.Date(get_bad_Rhats(fit, threshold = 1.1)))
            met <- extract_metab(fit, sitecode = site, bad_days = bad_days)
        }
        met <- select(met, site, date, GPP, ER, K600, DO_fit) %>%
            mutate(K_pool = K_pool,
                   K600_sigsig = K600_sigsig)
        comp_met <- bind_rows(comp_met, met)
    }

    COMP_MET <- bind_rows(COMP_MET, comp_met)
}

K_meds <- read_csv('data/metabolism/median_K600s.csv') %>%
    rename(kmed = K600)
COMP_MET <- COMP_MET %>%
    left_join(K_meds, by = 'site') %>%
    mutate(K600 = case_when(K_pool == 'kfull' ~ kmed,
                            TRUE ~ K600)) %>% select(-kmed)

KCOR <- COMP_MET %>% group_by(site, K600_sigsig, K_pool) %>%
    summarize(r2_er = cor(ER, K600, use = 'complete.obs'),
              r2_gpp = cor(GPP, K600, use = 'complete.obs')) %>%
    mutate(across(starts_with('r2_'), ~case_when(K_pool == 'kfull'~ NA_real_,
                                                TRUE ~ .)))

png('figures/SI/met_model_KER_relationships_comparison.png',
    width = 7, height = 5.5, units = 'in', res = 300)
    dd %>%
        mutate(K_group = paste(K_pool, K600_sigsig, sep = '_'),
               site = factor(site, level = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
        filter(K_group != 'kfull_0.0')%>%
        ggplot(aes(K600, ER, col = DO_fit)) +
        geom_point(size = 0.5) +
        scale_color_manual("Model fit",
                           values = c('brown3', 'black'))+
        facet_grid(site~K_group)+
        theme_bw()
dev.off()

cc <- COMP_MET %>%
    mutate(doy = as.numeric(format(date, "%j")),
           year = lubridate::year(date),
           model = paste(K_pool, K600_sigsig, sep = "_"),
           model = case_when(model == 'kb_0.1' ~ 'Binned K600',
                             model == 'knorm_0.01' ~ 'Normally pooled K600',
                             model == 'kfull_0.0' ~ 'Fixed K600',
                             TRUE ~ NA_character_),
           model = factor(model, levels = c('Binned K600',
                                            'Normally pooled K600',
                                            'Fixed K600')),
           site = factor(site, levels = c('PL', 'DL', 'GR',
                                          'GC', 'BM', 'BN'))) %>%
    filter(!is.na(model))
library(RColorBrewer)
a <- ggplot(cc, aes(date, GPP, col = model))+
    geom_line(linewidth = 0.5) +
    facet_grid(site~year, scales = 'free_x') +
    theme_classic()+
    scale_color_brewer('Model', palette = 'Set1')+
    theme(strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          panel.border = element_rect(fill = NA))
b <- cc %>%
    filter(model != 'Fixed K600') %>%
    ggplot(aes(K600, ER, col = model))+
    geom_point(size = 0.5) +
    facet_grid(site~year, scales = 'free_x') +
    scale_color_brewer('Model', palette = 'Set1')+
    theme_classic()+
    theme(panel.border = element_rect(fill = NA))

c <- KCOR %>%
    mutate(kgroup = paste(K_pool, K600_sigsig, sep = '_')) %>%
    ggplot( aes(K600_sigsig, r2_gpp, fill = kgroup)) +
    geom_bar(stat = 'identity',position = 'dodge',  width = 0.5) +
    geom_bar(aes(y = r2_er), stat = 'identity', position = 'dodge',  width = 0.5) +
    geom_hline(yintercept = 0)+
    facet_grid(site~.)+
    ylab('Correlation of K600 and GPP (positive) or ER (negative)')+
    ggtitle('')+
    theme_bw()

png('figures/SI/metab_change_model_choice.png', width = 9, height = 6,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE,
                      widths = c(2,1))#, align = 'v')
dev.off()

dat <- read_csv('data/prepared_data/compiled_prepared_data.csv')
dd <- left_join(COMP_MET, dat, by = c('site', 'date'))

COMP_MET %>%
    filter(site == 'GR') %>%
    mutate(doy = as.numeric(format(date, "%j")),
           year = lubridate::year(date)) %>%
    ggplot(aes(date, GPP, col = K600_sigsig, lty = K_pool))+
    geom_line(linewidth = 0.5) +
    facet_wrap(year~., scales = 'free_x', ncol = 1) +
    theme_bw()+
    scale_color_discrete('Prior for K600 sigma sigma')

dd %>%
    mutate(K_group = paste(K_pool, K600_sigsig, sep = '_'),
           site = factor(site, level = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    ggplot(aes(log(discharge), K600, col = DO_fit)) +
    geom_point() +
    facet_grid(site~K_group)+
    theme_bw()

add %>%
    mutate(K_group = paste(K_pool, K600_sigsig, sep = '_')) %>%
    group_by(K_group, site) %>%
    summarize(K600 = median(K600, na.rm = T)) %>%
    ggplot(aes(site, K600, fill = K_group)) +
    geom_bar(stat = 'identity', position = 'dodge') + theme_minimal()

k_meds <- dd %>% group_by(site) %>%
    summarize(K600 = median(K600, na.rm = T))
write_csv(k_meds, 'data/metabolism/median_K600s.csv')

KCOR %>%
    mutate(K600_sigsig = as.numeric(K600_sigsig)) %>%
    ggplot(aes(K600_sigsig, r2_er, col = K_pool, group = site))+
    geom_point() +
    geom_point(aes(y = r2_gpp)) +
    facet_wrap(.~site) + theme_bw() +
    geom_hline(yintercept = c(-0.4, 0.4))

COMP_MET %>%
    mutate(year = lubridate::year(date)) %>%
    tibble() %>%
    group_by(site, year) %>%
    summarize(GPP = sum(GPP, na.rm = T),
              ER = sum(ER, na.rm = T)) %>%
    mutate(NEP = GPP + ER,
           PR = -GPP/ER)

# Generate finalized file for results section
compiled_metab <- read_csv('data/metabolism/metabolism_compiled_all_sites_2000iter_kb_kss005.csv') %>%
    mutate(K600_sigsig = 0.05,
           K600_pool = 'binned_0.05')
compiled_metab <- read_csv('data/metabolism/metabolism_compiled_all_sites_2000iter_kn_kss01.csv') %>%
    mutate(K600_sigsig = 0.1,
           K600_pool = 'normal_0.1') %>%
    bind_rows(compiled_metab)

compiled_metab %>%
    ggplot(aes(K600, ER, col = DO_fit)) +
    geom_point() + facet_grid(site~K600_pool) +
    theme_bw()

met <- compiled_metab %>%
    select(-errors) %>%
    mutate(across(starts_with(c('GPP', 'ER', 'K600')),
                  ~case_when((!is.na(DO_fit) & DO_fit == 'bad') ~ NA_real_,
                             TRUE ~ .))) %>%
    select(-DO_fit)

dat <- read_csv('data/prepared_data/compiled_prepared_data.csv')
dd <- left_join(dat, met, by = c('site', 'date')) %>%
    select(-msgs.fit, -warnings, ends_with('Rhat') )

write_csv(dd, 'data/metabolism/metab_for_results.csv')


png('figures/K600xER_across_sites_knorm_K600sig01.png')
dd %>% ggplot(aes(K600, ER)) + geom_point() + facet_wrap(.~site)
dev.off()





# Comparison to Qipei's metabolism fit for Garrison
GR_20 <- read_csv('GR_metab.csv')
GR_2021 <- read_csv('GR_params_2020_2021.csv')
GR2 <-  GR_20 %>% select(date, GPP = GPP.daily, ER = ER.daily, K600 = K600.daily) %>%
    mutate(met = 'sensor light') %>%
    bind_rows(filter(COMP_MET, site == 'GR', K600_sigsig %in% c('0.0', '0.01'),
                     K_pool != 'kb'))
GR2 <-  mutate(GR2, est = case_when(is.na(met)~'theoretical light',
                               TRUE ~ met)) %>% select(-met)
# GR2 <-  GR_2021 %>% select(date, GPP = GPP.daily, ER = ER.daily, K600 = K600.daily) %>%
#     mutate(met = 'qipei') %>%
#     bind_rows(GR2) %>%
#     mutate(year = year(date),
#                est = case_when(is.na(met)~'alice',
#                                TRUE ~ met)) %>% select(-met)

GR2 %>%
    filter(date < ymd('2020-10-08')&
               date > ymd('2020-08-20')) %>%
    pivot_longer(cols = c('GPP', 'ER', 'K600'),
                 names_to = 'met', values_to = 'value') %>%
    mutate(K600_sigsig = case_when(is.na(K600_sigsig) ~ "0.0",
                                   TRUE ~ K600_sigsig))%>%
    ggplot(aes(date, value, col = est, lty = K600_sigsig)) +
        geom_line() +
        facet_wrap(.~met, scales = 'free', ncol = 1, strip.position = 'right')+
        scale_color_manual(values = c('brown','black'))+
        theme_bw()

GR2 %>%
    ggplot(aes(K600, ER, col = est)) +
        geom_point() + theme_bw()

GR <-  GR_20 %>% select(date, GPP_CO2 = GPP.daily, ER_CO2 = ER.daily, K600_CO2 = K600.daily) %>%
     left_join(filter(COMP_MET, site == 'GR'), by = 'date')
GR %>%
    rename(GPP_O2 = GPP, ER_O2 = ER) %>%
    select(-starts_with('K')) %>%
    pivot_longer(starts_with(c('ER', 'GPP')), names_to = c('met', 'gas'),
                 values_to = 'gm2', names_sep = '_') %>%
    pivot_wider(names_from = met, values_from = gm2)
ggplot(GR, aes(date)) +
    geom_line(aes(y = GPP, lty = K600_sigsig)) +
    geom_line(aes(y = ER, lty = K600_sigsig)) +
    geom_line(aes(y = GPP_CO2), col = 'brown') +
    geom_line(aes(y = ER_CO2), col = 'brown') +
    geom_hline(yintercept = 0)+
    ylab('metabolism')+
    theme_classic()
er <- GR2 %>% mutate(year = factor(lubridate::year(date))) %>%
ggplot(aes(K600, ER, col = year)) +
    geom_point() +
    facet_wrap(.~met, scales = 'free')+
    theme_bw()
gpp <- GR2 %>% mutate(year = factor(lubridate::year(date))) %>%
ggplot(aes(K600, GPP, col = year)) +
    geom_point() +
    facet_wrap(.~met, scales = 'free')+
    theme_bw()

ggpubr::ggarrange(gpp, er, ncol = 1, common.legend = T)
hist(GR$K600_CO2)
