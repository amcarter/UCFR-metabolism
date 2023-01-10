# Compare across metabolism fits with varying K600 sigma sigma

library(tidyverse)

files <- grep('\\.rds$', list.files('data/metab_fits'), value = TRUE)

COMP_MET <- data.frame()
KCOR <- data.frame()
for(site in c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')){
    comp_met <- data.frame()

    sitefiles <- grep(paste0('^', site), files, value = TRUE)

    for(f in sitefiles){
        fit <- readRDS(paste0('data/metab_fits/', f))
        K600_sigsig <- paste0('0.', str_match(f, '_([0-9]+)\\.rds$')[,2])
        if(K600_sigsig == '0.0'){
                met <- fit@fit %>%
                    select(date, GPP = GPP.daily, ER = ER.daily,
                           K600 = K600.daily) %>%
                    mutate(site = site)
        } else{
            met <- extract_metab(fit, sitecode = site)
        }
        met$K600_sigsig <- K600_sigsig
        met <- select(met, site, date, GPP, ER, K600, K600_sigsig)
        comp_met <- bind_rows(comp_met, met)
    }

    kcor <- data.frame()
    for(k in unique(comp_met$K600_sigsig)){
        kk <- filter(comp_met, K600_sigsig == k)
        if(k == '0.0'){ r2_er = r2_gpp = 0} else{
            r2_er <- cor(kk$ER, y = kk$K600, use = 'complete.obs')
            r2_gpp <- cor(kk$GPP, y = kk$K600, use = 'complete.obs')
        }
        kk <- data.frame(K600_sigsig = k,
                         r2_er = r2_er,
                         r2_gpp = r2_gpp)
        kcor <- bind_rows(kcor, kk)
    }
    kcor$site = site

    COMP_MET <- bind_rows(COMP_MET, comp_met)
    KCOR = bind_rows(KCOR, kcor)
}

a <- COMP_MET %>%
    mutate(doy = as.numeric(format(date, "%j")),
           year = year(date)) %>%
    ggplot(aes(date, GPP, col = K600_sigsig))+
    geom_line(linewidth = 0.5) +
    facet_grid(site~year, scales = 'free_x') +
    theme_bw()+
    theme(strip.background.y = element_blank(),
          strip.text.y = element_blank())+
    scale_color_discrete('Prior for K600 sigma sigma')
b <- ggplot(KCOR, aes(K600_sigsig, r2_gpp, fill = K600_sigsig)) +
    geom_bar(stat = 'identity', width = 0.5) +
    geom_bar(aes(y = r2_er), stat = 'identity', width = 0.5) +
    geom_hline(yintercept = 0)+
    facet_grid(site~.)+
    ylab('Correlation of K600 and GPP (positive) or ER (negative)')+
    ggtitle('')+
    theme_bw()

png('figures/SI/metab_change_with_K600_sigsig.png', width = 10, height = 7,
    units = 'in', res = 300)
    ggpubr::ggarrange(a,b, nrow = 1, common.legend = TRUE, widths = c(3, 1), align = 'v')
dev.off()
