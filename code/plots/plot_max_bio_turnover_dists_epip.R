library(tidyverse)
library(rstan)
library(lme4)
library(brms)
# Growth rate calculations ####

# read in datasets:

coefs <- readRDS("data/full_biomass_NPP_partition_mod_chla_posterior_coefs_epip.rds")
bm_met <- read_csv("data/biomass_data/log_gamma_brms_gam_fits_biomass_epip.csv")

biogams <- read_csv('data/biomass_data/log_gamma_brms_gam_fits_biomass_epip.csv') %>%
    mutate(site = ifelse(site == "BG", "BM", site))
bio <- read_csv('data/biomass_data/biomass_working_data_summary_epip.csv') %>%
    mutate(epi_gm2 = epilithon_gm2 + epiphyton_gm2,
           epi_chla_mgm2 = epilithon_chla_mgm2 + epiphyton_chla_mgm2,
           filaepip_gm2 = epiphyton_gm2 + filamentous_gm2,
           filaepip_chla_mgm2 = epiphyton_chla_mgm2 + filamentous_chla_mgm2)

# Calculate a chla -> biomass conversion using the measured data:
bio_full <- bio %>%
    mutate(
        filamentous_chla_mgm2 = ifelse(filamentous_gm2 <0.75, NA, filamentous_chla_mgm2 ),
        filamentous_gm2 = ifelse(filamentous_gm2 <0.75, NA, filamentous_gm2 ),
        fila_chla_per_gm = filamentous_chla_mgm2/filamentous_gm2,
        epiphyton_chla_per_gm = epiphyton_chla_mgm2/epiphyton_gm2,
        epilithon_chla_per_gm = epilithon_chla_mgm2/epilithon_gm2,
        epiphyton_chla_per_gm = ifelse(is.infinite(epiphyton_chla_per_gm), NA, epiphyton_chla_per_gm)) %>%
    right_join(select(biogams, site, date, PAR_bc_Jm2,
                      starts_with(c('frac_fila', 'NPP')),
                      epil_chla_mgm2_fit = epilithon_chla_mgm2,
                      epil_chla_mgm2_lower = epilithon_chla_lower,
                      epil_chla_mgm2_upper = epilithon_chla_upper,
                      epip_chla_mgm2_fit = epiphyton_chla_mgm2,
                      epip_chla_mgm2_lower = epiphyton_chla_lower,
                      epip_chla_mgm2_upper = epiphyton_chla_upper,
                      fila_chla_mgm2_fit = filamentous_chla_mgm2,
                      fila_chla_mgm2_lower = filamentous_chla_lower,
                      fila_chla_mgm2_upper = filamentous_chla_upper
    )) %>%
    group_by(site) %>% arrange(date) %>%
    mutate(across(ends_with('chla_per_gm'), ~zoo::na.approx(., na.rm = FALSE, x = date, rule = 2)),
           epil_gm2_est = case_when(epil_chla_mgm2_fit == 0 ~ 0,
                                    TRUE ~ epil_chla_mgm2_fit/epilithon_chla_per_gm),
           epil_gm2_lower = epil_chla_mgm2_lower/epilithon_chla_per_gm,
           epil_gm2_upper = epil_chla_mgm2_upper/epilithon_chla_per_gm,
           epip_gm2_est = case_when(epip_chla_mgm2_fit == 0 ~ 0,
                                    TRUE ~ epip_chla_mgm2_fit/epiphyton_chla_per_gm),
           epip_gm2_lower = epip_chla_mgm2_lower/epiphyton_chla_per_gm,
           epip_gm2_upper = epip_chla_mgm2_upper/epiphyton_chla_per_gm,
           fila_gm2_est = case_when(fila_chla_mgm2_fit == 0 ~ 0,
                                    TRUE ~ fila_chla_mgm2_fit/fila_chla_per_gm),
           fila_gm2_lower = fila_chla_mgm2_lower/fila_chla_per_gm,
           fila_gm2_upper = fila_chla_mgm2_upper/fila_chla_per_gm
    ) %>%
    select(-starts_with(c('epi_', 'filaepip_', 'pct_')))

bio_full <- bio_full[,-c(3:24)] %>%
    mutate(doy = as.numeric(format(date, '%j')),
           year = year(date))
# do calculations with bayesian model output:
# coefs_afdm <- coefs$coefs_afdm %>%
#     mutate(across(starts_with("mu"),
#                   ~ .x * 2 *12/32), # convert to units of 1/d
#            sigma = sigma * 12/32)
coefs_chl <- coefs$coefs_chla %>%
    mutate(across(starts_with("mu"),
                  ~ .x * 12/32),   # convert to units of mg C/mg chl a/d
           sigma = sigma * 12/32)

coefs_chl
data.frame(t(coefs_chl)) %>% mutate(model = "chl")

mu_epil_c <- coefs$coefs_chla$mu_epil * 12/32
mu_epip_c <- coefs$coefs_chla$mu_epip * 12/32
mu_fila_c <- coefs$coefs_chla$mu_fila * 12/32
fila_chla_C <- bm_met$filamentous_chla_mgm2/bm_met$filamentous_gm2 * 2
epip_chla_C <- bm_met$epiphyton_chla_mgm2/bm_met$epiphyton_gm2 * 2
epil_chla_C <- bm_met$epilithon_chla_mgm2/bm_met$epilithon_gm2 * 2
fila_chla_C <- if_else(is.infinite(fila_chla_C), NA, fila_chla_C)
epip_chla_C <- if_else(is.infinite(epip_chla_C), NA, epip_chla_C)
epil_chla_C <- if_else(is.infinite(epil_chla_C), NA, epil_chla_C)

mu_epil_cC <- mu_epil_c * mean(epil_chla_C, na.rm = T)
mu_epip_cC <- mu_epip_c * mean(epip_chla_C, na.rm = T)
mu_fila_cC <- mu_fila_c * mean(fila_chla_C, na.rm = T)
mu_epil_cC/mu_fila_cC

bm_met_chl <- bio_full %>%
    mutate(light = PAR_bc_Jm2/max(PAR_bc_Jm2)) %>%
    select(site, date, doy, year, light,
           starts_with(c("NPP", "epil", "epip", "fila", "frac"))) %>%
    mutate(lscale = light/(light + 0.5),
           mu_epil_med = mu_epil_c[1]*lscale,
           mu_epil_low = mu_epil_c[2]*lscale,
           mu_epil_high = mu_epil_c[3]*lscale,
           mu_epip_med = mu_epip_c[1]*lscale,
           mu_epip_low = mu_epip_c[2]*lscale,
           mu_epip_high = mu_epip_c[3]*lscale,
           mu_fila_med = mu_fila_c[1]*lscale,
           mu_fila_low = mu_fila_c[2]*lscale,
           mu_fila_high = mu_fila_c[3]*lscale,
           fila_prodgCd_med = fila_chla_mgm2_fit * mu_fila_c[1] * lscale,
           fila_prodgCd_low = fila_chla_mgm2_lower * mu_fila_c[2] * lscale,
           fila_prodgCd_high = fila_chla_mgm2_upper * mu_fila_c[3] * lscale,
           epip_prodgCd_med = epip_chla_mgm2_fit * mu_epip_c[1] * lscale,
           epip_prodgCd_low = epip_chla_mgm2_upper * mu_epip_c[3] * lscale,
           epip_prodgCd_high = epip_chla_mgm2_lower * mu_epip_c[2] * lscale,
           epil_prodgCd_med = epil_chla_mgm2_fit * mu_epil_c[1] * lscale,
           epil_prodgCd_low = epil_chla_mgm2_upper * mu_epil_c[3] * lscale,
           epil_prodgCd_high = epil_chla_mgm2_lower * mu_epil_c[2] * lscale,
           fila_turn_med = fila_gm2_est/fila_prodgCd_med/2,
           fila_turn_high = fila_gm2_lower/fila_prodgCd_high/2,
           fila_turn_low = fila_gm2_upper/fila_prodgCd_low/2,
           epip_turn_med = epip_gm2_est/epip_prodgCd_med/2,
           epip_turn_high = epip_gm2_lower/epip_prodgCd_high/2,
           epip_turn_low = epip_gm2_upper/epip_prodgCd_low/2,
           epil_turn_med = epil_gm2_est/epil_prodgCd_med/2,
           epil_turn_high = epil_gm2_lower/epil_prodgCd_high/2,
           epil_turn_low = epil_gm2_upper/epil_prodgCd_low/2
           )

write_csv(bm_met_chl, 'data/biomass_turnover_prod_chla_full_error_propagation_epip.csv')

# bm_met <- bm_met %>%
#     mutate(fila_prod_gCd = (coefs[2] * fila_chla_mgm2*light) * 12/32 ,
#            fila_turnover = case_when(fila_gm2 > min_fila_gm2 ~ fila_gm2/2/fila_prod_gCd,
#                                      TRUE ~ NA_real_),
#            epil_prod_gCd = (coefs[1] * epil_chla_mgm2*light) * 12/32 ,
#            epil_turnover = epil_gm2/2/epil_prod_gCd)
# bm_met %>%
#     filter(year == 2021) %>%
#     select(-year, -ARf, -PAR_surface) %>%
#     write_csv('data/biomass_data/2021_turnovers_for_Rafa.csv')

# bm_met <- bm_met %>%
#     mutate(fila_prod_gCd = (0.02394 * fila_gm2*light)*14/32 ,
#            fila_turnover = case_when(fila_prod_gCd > 0.04 ~ fila_gm2/2/fila_prod_gCd,
#                                      TRUE ~ NA_real_),
#            epil_prod_gCd = (0.435 * epil_gm2*light)*14/32 ,
#            epil_turnover = epil_gm2/2/epil_prod_gCd)
bloom_df <- data.frame(site = rep(c("PL","DL","GR", "GC", "BG", "BN"), each = 2),
                       year = rep(c(2020, 2021), 6),
                       bloom = c(0,0,0,0,1,0,1,1,0,1,0,1)) %>%
    mutate(fila_bloom = if_else(bloom == 0, "No Bloom", "Bloom"))


min_fila_gm2 <- 0.36129
bm_fila_chl <- filter(bm_met_chl, fila_gm2_est > min_fila_gm2) %>%
    pivot_longer(cols = ends_with(c('med', 'low', 'high')),
                 values_to = 'value', names_to = c('Biomass', 'measure', 'stat'),
                 names_pattern = '(epil|epip|fila|mu)_(.*)_(med|high|low)')  %>%
    filter(Biomass != 'mu') %>%
    left_join(bloom_df, by = c("site", "year")) %>%
    mutate(value = if_else(value > 10000, NA, value))

bm_fila_chl %>%
    group_by(Biomass, measure, stat) %>%
    summarize(value = mean(value, na.rm = T))

p1 <- bm_fila_chl%>%
    mutate(Rate = 'chl') %>%
    mutate(rate_biomass = paste0(Rate, Biomass)) %>%
    filter(measure == 'turn') %>%
    ggplot(aes(x=value, group = rate_biomass, fill = Biomass)) +
    geom_density(adjust=1.5, alpha=.4) +
    # scale_fill_manual(values = c('#1B9EC9', '#97BB43'))+
    scale_x_log10(limits = c(2, 1500))+
    ylab('Density')+
    xlab('Turnover time (d)')+
    theme_classic()+
    theme(legend.position = 'inside',
          legend.position.inside = c(0.9,0.6),
          panel.border = element_rect(fill = NA))
# p1 <- p1 + annotate(geom = 'text', x = 12, y = 0.25,
#                     label="Epilithic", col = '#1B9EC9')
# p1 <- p1 + annotate(geom = 'text', x = 32, y = 0.06,
#                     label="Filamentous", col = '#97BB43')
p2 <- bm_fila_chl %>% filter(!is.na(fila_bloom)) %>%
    ggplot(aes(frac_fila, value))+
    xlab('Filamentous fraction of \ntotal biomass (%)')+
    ylab('Turnover time (d)')+
    geom_boxplot(aes(pch = factor(fila_bloom), col = factor(fila_bloom)),
               size = 1.6)+
    theme_classic()+
    theme(legend.title = element_blank(),
          legend.position = c(0.22, 0.85),
          panel.border = element_rect(fill = NA))
png('figures/turnover_time_by_composition.png',
    width = 6.5, height = 3.5, units = 'in', res = 300)
ggpubr::ggarrange(p1, p2, common.legend = TRUE)
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
coeff = 150
ann_text2 <- read_csv('data/metabolism/auto_sites_figure_labels.csv') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')),
           year = factor(year),
           date = case_when(year == 2020 ~ as.Date('2020-10-20'),
                            TRUE ~ as.Date('2021-10-14')),
           prodgCd_med = rep(0.5, 12), gCm2 = rep(1, 12))
ann_text2.1 <- ann_text2 %>%
    mutate(date = case_when(year == 2020 ~ as.Date('2020-07-16'),
                            TRUE ~ date),
           sitename = case_when(year == 2020 ~ site,
                                TRUE ~ ""))

p <-
    bm_met_C %>%
        group_by(site, year) %>%
        mutate(across(c(filamentous_gm2, epilithon_gm2),
                      ~ zoo::rollmean(.x, k = 7, na.pad = TRUE))) %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR','GC','BG','BN'))) %>%
    mutate(filamentous_gm2 = case_when(filamentous_gm2 < min_fila_gm2 ~ NA_real_,
                                TRUE ~ filamentous_gm2)) %>%
        ungroup() %>%
    select(site, date, year, light,
           starts_with(c('fila_prodgCd', 'epil_prodgCd')),
                       fila_gm2_med = filamentous_gm2, epil_gm2_med = epilithon_gm2) %>%
    pivot_longer(cols = starts_with(c('fila', 'epil')),
                 values_to = 'value', names_to = c('biomass', 'measure'),
                 names_pattern = '(fila|epil).*_(prod.*|gm2)') %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(Light = factor(rep(" ", 2*nrow(bm_met)), levels = c("1"," ")),
           gCm2 = gm2/2,
           biomass = case_when(biomass == 'fila' ~ 'Filamentous',
                               biomass == 'epil'~'Epilithic'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
    mutate(prodgCd_med = if_else(prodgCd_med > 2, 2, prodgCd_med)) %>%
    ggplot(aes(date, prodgCd_med/gCm2, col = biomass))+
    geom_area(aes(date, gCm2/coeff, fill = biomass), color = NA, alpha = 0.4)+
    geom_line(size = 1.2)+
    geom_line(aes(y = light/2.5, lty = Light), col = 'grey20') +
    facet_grid(site~year, scales = 'free_x', ) +
    geom_text(data = ann_text2, aes(label = trophic), col = 'black') +
    geom_text(data = ann_text2.1, aes(label = sitename), col = 'black') +
    scale_y_continuous(
        expand = expand_scale(mult = c(0.1, 0.1)),
        name = expression(paste('Production rate (', d^-1, ')')),
        sec.axis = sec_axis(~.*coeff,
                            name = expression(paste('Biomass (g C', m^-2, ')')),
                            breaks = c(0, 25, 50))
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
png('figures/biomass_prod_and_turnover.png', width = 6, height = 8,
    units = 'in', res = 300)
p
dev.off()
p <-
    bm_met_C %>%
        group_by(site, year) %>%
        mutate(across(c(filamentous_gm2, epilithon_gm2),
                      ~ zoo::rollmean(.x, k = 7, na.pad = TRUE))) %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR','GC','BG','BN'))) %>%
    mutate(filamentous_gm2 = case_when(filamentous_gm2 < min_fila_gm2 ~ NA_real_,
                                TRUE ~ filamentous_gm2)) %>%
        ungroup() %>%
    select(site, date, year, light,
           starts_with(c('fila_prodgCd', 'epil_prodgCd')),
                       fila_gm2_med = filamentous_gm2, epil_gm2_med = epilithon_gm2) %>%
    pivot_longer(cols = starts_with(c('fila', 'epil')),
                 values_to = 'value', names_to = c('biomass', 'measure'),
                 names_pattern = '(fila|epil).*_(prod.*|gm2)') %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(Light = factor(rep(" ", 2*nrow(bm_met)), levels = c("1"," ")),
           gCm2 = gm2/2,
           biomass = case_when(biomass == 'fila' ~ 'Filamentous',
                               biomass == 'epil'~'Epilithic'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
    mutate(prodgCd_med = if_else(prodgCd_med > 2, 2, prodgCd_med)) %>%
    ggplot(aes(date, prodgCd_med/gCm2, col = biomass))+
    geom_area(aes(date, gCm2/coeff, fill = biomass), color = NA, alpha = 0.4)+
    geom_line(size = 1.2)+
    geom_line(aes(y = light/2.5, lty = Light), col = 'grey20') +
    facet_grid(site~year, scales = 'free_x', ) +
    geom_text(data = ann_text2, aes(label = trophic), col = 'black') +
    geom_text(data = ann_text2.1, aes(label = sitename), col = 'black') +
    scale_y_continuous(
        expand = expand_scale(mult = c(0.1, 0.1)),
        name = expression(paste('Production rate (', d^-1, ')')),
        sec.axis = sec_axis(~.*coeff,
                            name = expression(paste('Biomass (g C', m^-2, ')')),
                            breaks = c(0, 25, 50))
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
          axis.title = element_text(size = 14),
          # axis.text = element_text(size = 8),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
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
tiff('figures/biomass_prod_and_turnover.tiff', width = 18, height = 24,
    units = 'cm', res = 600)
p
dev.off()


png('figures/percent_standing_crop_and_prod.png', width = 6, height = 3,
    units = 'in', res = 300, type = 'cairo')
met_fig_dat <- bm_met %>%
    mutate(doy = as.numeric(format(date, '%j')),
           Date = as.Date(paste0('2020-', doy), format = '%Y-%j'),
           site = factor(site, levels = c('PL', 'DL', 'GR','GC','BG','BN'))) %>%
    mutate(fila_gm2 = case_when(fila_gm2 < min_fila_gm2 ~ NA_real_,
                                TRUE ~ fila_gm2)) %>%
    select(site, Date, year, light,
           fila_prod_gCd, fila_gm2, epil_prod_gCd, epil_gm2) %>%
    pivot_longer(cols = starts_with(c('fila', 'epil')),
                 values_to = 'value', names_to = c('biomass', 'measure'),
                 names_pattern = '(fila|epil)_(.*)') %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(Light = factor(rep(" ", 2*nrow(bm_met)), levels = c("1","")),
           gCm2 = gm2/2,
           biomass = case_when(biomass == 'fila' ~ 'Filamentous',
                               biomass == 'epil'~'Epilithic'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')),
           p_rate = prod_gCd/gCm2
    )%>%
    pivot_longer(cols = c('prod_gCd', 'gCm2'), values_to = 'Value',
                 names_to = 'measure') %>%
    mutate(measure = factor(measure, levels = c('gCm2', 'prod_gCd')))

per_means1 <- met_fig_dat %>%
    group_by(biomass, measure) %>%
    summarize(Value = mean(Value, na.rm = T))
per_means2 <- per_means1%>%
    pivot_wider(names_from = measure, values_from = Value) %>%
    mutate(gCm2 = 100*gCm2/19.04,
           prod_gCd = 100*prod_gCd/1.978) %>%
    pivot_longer(cols = c('gCm2', 'prod_gCd'),
                 names_to = 'measure', values_to = 'percent') %>%
    left_join(per_means1)

ggplot(met_fig_dat, aes(biomass, Value, fill = biomass))+
    geom_violin(alpha = 0.4, draw_quantiles = c(0.5)) +
    stat_summary(fun = mean, geom = "point", shape = 8, size = 2,
                 color = "black", show.legend = FALSE) +  # Mean point
    geom_text(
        aes(label = paste0(round(percent, 0), '%')),  # Label with rounded median values
        data = per_means2,
        vjust = -0.3,  # Adjust vertical position of text labels
        hjust = -0.35,   # Center text horizontally
        size = 2.75,      # Adjust text size as needed
        color = "black",
        fontface = "bold",
        position = position_dodge(width = 0.75)
    ) +
    scale_fill_manual('Biomass form',
                      values = c('#97BB43', '#1B9EC9')) +
    facet_wrap(measure~., scales = 'free',
               strip.position = 'left',
               labeller = as_labeller(c(gCm2 = 'Standing~Crop~(g~C~m^{-2})',
                                        prod_gCd = 'P[N]~(g~C~m^{-2}~d^{-1})'),
                                      default = label_parsed)) +
    theme_bw()+
    labs(y = NULL, x = NULL)+
    theme(legend.position = 'none',
          strip.background = element_blank(),
          strip.placement = 'outside')

dev.off()







ann_text <- data.frame(Date = rep(as.Date('2020-07-16'), 6),
                       prod_gCd = rep(0.48, 6),
                       gCm2 = rep(1, 6),
                       year = rep(2020,6),
                       biomass = factor(rep('Epilithic', 6),
                                        levels = c('Filamentous','Epilithic')),
                       site = factor(c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'),
                                     levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))


p3 <- p + geom_text(data = ann_text, aes(label = site), col = 'black')

bm_sum <-  bm_met_C %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR','GC','BG','BN'))) %>%
    group_by(site, year) %>%
    summarize(n = n(),
              fila_Biomass_med = quantile(filamentous_gm2, 0.975)/2,
              epil_Biomass_med = quantile(epilithon_gm2, 0.975)/2,
              fila_Biomass_high = quantile(filamentous_gm2_upper, 0.975)/2,
              epil_Biomass_high = quantile(epilithon_gm2_upper, 0.975)/2,
              fila_Biomass_low = quantile(filamentous_gm2_lower, 0.975)/2,
              epil_Biomass_low = quantile(epilithon_gm2_lower, 0.975)/2,
              fila_cumprod_med = sum(fila_prodgCd_med)/n*100,
              fila_cumprod_high = sum(fila_prodgCd_low)/n*100,
              fila_cumprod_low = sum(fila_prodgCd_high)/n*100,
              epil_cumprod_med = sum(epil_prodgCd_med)/n*100,
              epil_cumprod_high = sum(epil_prodgCd_low)/n*100,
              epil_cumprod_low = sum(epil_prodgCd_high)/n*100,
              # cumprod_gC = fila_cumprod_med + epil_cumprod_med,
              frac_fila = mean(frac_fila, na.rm =T)) %>%
    mutate(epil_cumprod_high = if_else(epil_cumprod_high>250,
                                       250, epil_cumprod_high)) %>%
    left_join(bloom_df, by = c('site', 'year'))
# p2 <- ggplot(bm_sum, aes(frac_fila, cumprod_gC))+
#     xlab('Filamentous fraction of \ntotal biomass (%)')+
#     ylab('Turnover time (d)')+
#     geom_point(aes(pch = factor(bloom)), size = 1.6, col = 'black')+
#     theme_classic()+
#     theme(legend.title = element_blank(),
#           legend.position = c(0.22, 0.85),
#           panel.border = element_rect(fill = NA))
jitterer <- position_jitter(width = 0.25, seed = 123)
p4 <- bm_sum %>%
    select(-bloom) %>% rename(bloom = fila_bloom) %>%
        mutate(bloom = factor(bloom, levels = c("No Bloom", "Bloom"))) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
                 values_to = 'value',
                 names_to = c('biomass', 'measure', 'stat'),
                 names_pattern = '(fila|epil)_(Biomass|cumprod)_(med|low|high)') %>%
    mutate(measure = case_when(measure == 'cumprod' ~ 'Cumulative \nProduction',
                               measure == 'Biomass' ~ 'Maximum \nBiomass'),
           measure = factor(measure, levels = c('Maximum \nBiomass',
                                                'Cumulative \nProduction')),
           biomass = case_when(biomass == 'epil' ~ 'Epilithic',
                               biomass == 'fila' ~ 'Filamentous'),
           biomass = factor(biomass, levels = c('Filamentous', 'Epilithic')))%>%
    pivot_wider(values_from = value, names_from = stat) %>%
    ggplot(aes(x = biomass, y = med, fill = biomass)) +
    geom_boxplot(outlier.shape=NA, alpha = 0.2)+
    geom_errorbar(aes(ymin = low, ymax = high, col = biomass), width = 0,
                  position = jitterer)+
    geom_point(aes(pch = bloom, group = biomass),
               position = jitterer)+
        facet_grid(.~measure)+
    # geom_jitter(color="black", alpha=0.9) +
    scale_color_manual('',
                      values = c('#97BB43', '#1B9EC9')) +
    scale_fill_manual('',
                      values = c('#97BB43', '#1B9EC9')) +
    ylab(expression(paste('Algal Biomass (g C ', m^-2, ')')))+
    xlab('')+
    theme_classic()+
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.spacing.y = unit(-0.1, "cm"),
          legend.position = c(0.24, 0.79),
          legend.background = element_rect(fill = NA))
p1 <- bm_fila_chl%>%
    mutate(Rate = 'chl') %>%
    bind_rows(mutate(bm_fila_C, Rate = "C")) %>%
    mutate(rate_biomass = paste0(Rate, Biomass),
           Biomass = case_when(Biomass == "epil" ~ "Epilithic",
                               Biomass == "fila" ~ "Filamentous")) %>%
    filter(measure == 'turn') %>%
    ggplot(aes(x=value, group = rate_biomass, fill = Biomass, lty = Rate)) +
    geom_density(adjust=1.5, alpha=.4) +
    scale_fill_manual(values = c('#1B9EC9', '#97BB43'))+
    scale_x_log10(limits = c(2, 1500))+
    ylab('Density')+
    xlab('Turnover time (d)')+
    theme_classic()+
    theme(legend.position = 'inside',
          legend.position.inside = c(0.7,0.6),
          panel.border = element_rect(fill = NA))
png('figures/biomass_cumulative_and_turnover.png', width = 6.5, height = 3.5,
    units = 'in', res = 300)
ggpubr::ggarrange(p4, p1, ncol = 2,
                  labels = c('(a)', '(b)'),
                  align = 'h', widths = c(1.5, 1))

dev.off()
tiff('figures/biomass_cumulative_and_turnover.tiff', width = 18, height = 9.7,
    units = 'cm', res = 600)
ggpubr::ggarrange(p4, p1, ncol = 2,
                  labels = c('(a)', '(b)'),
                  align = 'h', widths = c(1.5, 1))

dev.off()

bm_sum %>% group_by(bloom) %>%
    summarize(across(where(is.numeric), .fns = c(mean = ~mean(.), sd = ~sd(.))))

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


# Numbers for results:
dd <- bm_met %>% group_by(site, year) %>%
    summarize(NPP = mean(NPP, na.rm = T )*12/32,
              epil_prod_gCd = mean(epil_prod_gCd, na.rm = T),
              fila_prod_gCd = mean(fila_prod_gCd, na.rm = T))
summary(dd)
sd(dd$NPP)/sqrt(12)

mean(bm_met$epil_prod_gCd, na.rm = T)
calculate_ts_mean_se(bm_met$epil_prod_gCd)
mean(bm_met$fila_prod_gCd, na.rm = T)
calculate_ts_mean_se(bm_met$fila_prod_gCd)

bm_rate <- bm_met %>%
    group_by(site, year) %>%
    mutate(fila_rate = fila_prod_gCd/(fila_gm2/2),
           epil_rate = epil_prod_gCd/(epil_gm2/2)) %>%
    filter(fila_gm2 >= min_fila_gm2)
summary(bm_rate)

calculate_ts_mean_se(bm_rate$epil_rate)
calculate_ts_mean_se(bm_rate$fila_rate)

bm_prod <- bm_met %>%
    group_by(site, year) %>%
    summarize(n = n(),
              fila_Biomass = max(fila_gm2)/2,
              epil_Biomass = max(epil_gm2)/2,
              fila_prod = mean(fila_prod_gCd),
              epil_prod = mean(epil_prod_gCd),
              fila_cumprod = sum(fila_prod_gCd)*100/n,
              epil_cumprod = sum(epil_prod_gCd)*100/n) %>%
    mutate(n_epil = epil_cumprod/epil_Biomass,
           n_fila = fila_cumprod/fila_Biomass)

summary(bm_prod)
sd(bm_prod$epil_cumprod)/sqrt(12)

bm_bloom <-  filter(bm_prod, (year == 2020 & site %in% c('GR', 'GC'))|
                        (year == 2021 & site %in% c('GC', 'BN','BG')))
bm_notbloom <-  filter(bm_prod, (year == 2020 & !(site %in% c('GR', 'GC')))|
                           (year == 2021 & !(site %in% c('GC', 'BN','BG'))))

summary(bm_bloom)
sd(bm_bloom$fila_cumprod)/sqrt(5)
summary(bm_notbloom)
sd(bm_notbloom$fila_cumprod)/sqrt(7)

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

dd <- bm_met %>%
    select(site, year, date, epil_gm2, fila_gm2, GPP, ER, NPP, ARf,
           fila_prod_gCd, epil_prod_gCd) %>%
    mutate(GPP = GPP * 12/32,
           ER = ER *12/32,
           NPP = NPP *12/32,
           AR = GPP*ARf,
           HR = ER+AR,
           fila_gC = fila_gm2/2,
           epil_gC = epil_gm2/2)

dd$fila_diff = c(0, diff(dd$fila_gC) - dd$fila_prod_gCd[1:(nrow(dd)-1)])
dd$epil_diff = c(0, diff(dd$epil_gC) - dd$epil_prod_gCd[1:(nrow(dd)-1)])

ids <- dd %>%
    mutate(r_number = row_number()) %>%
    group_by(site, year) %>%
    summarise(index = min(r_number))
dd$epil_diff[ids$index] <- dd$fila_diff[ids$index] <- NA

dd %>%
    mutate(algal_loss = fila_diff + epil_diff - HR) %>%
    pivot_longer(cols = ends_with('diff'),
                 values_to = 'loss', names_to = 'algae') %>%
    mutate(loss = case_when(loss>0 ~ 0,
                            TRUE ~ loss))%>%
    ggplot(aes(date, loss, fill = algae))+
    geom_area()+
    geom_line(aes(y = HR))+
    facet_grid(site~year, scales = 'free_x') +
    theme_classic()

ggplot(data, aes(x=time, y=value, fill=group)) +
    geom_area()

filter(dd, site == 'DL', year == 2020)
fila_prod_gCd/(fila_prod_gCd + epil_prod_gCd),
epil_AR = AR * epil_prod_gCd/(fila_prod_gCd + epil_prod_gCd))



