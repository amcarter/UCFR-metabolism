# Calculate biomass numbers for results
library(tidyverse)
source('code/UCFR_helpers.R')

met <- read_csv('data/metabolism/metab_for_results.csv') %>%
    mutate(year = year(date),
           NEP = GPP + ER,
           PR = -GPP/ER,
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           site_year = paste(site, year, sep = '_'))%>%
    filter(!is.na(year))

head(met)

# results describing means/ranges:
min(met$GPP, na.rm = T)
max(met$GPP, na.rm = T)
mean(met$GPP, na.rm = T)
sd(met$GPP, na.rm = T)
calculate_ts_mean_se(met$GPP)

min(met$ER, na.rm = T)
max(met$ER, na.rm = T)
mean(met$ER, na.rm = T)
calculate_ts_mean_se(met$ER)

min(met$NEP, na.rm = T)
max(met$NEP, na.rm = T)
mean(met$NEP, na.rm = T)
calculate_ts_mean_se(met$NEP)

min(met$PR, na.rm = T)
max(met$PR, na.rm = T)
mean(met$PR, na.rm = T)
calculate_ts_mean_se(met$PR)

met %>% group_by(site_year) %>%
    summarize(across(any_of(c('GPP', 'ER', 'NEP', 'PR')), ~mean(., na.rm = T)))

year_sum <- data.frame(site_year = unique(met$site_year),
                       PRcor = NA_real_,
                       GPP = NA_real_,
                       GPP.se = NA_real_,
                       ER = NA_real_,
                       ER.se = NA_real_,
                       NEP = NA_real_,
                       NEP.se = NA_real_,
                       PR = NA_real_,
                       PR.se = NA_real_)

for(i in 1:nrow(year_sum)){
    mm <- filter(met, site_year == year_sum$site_year[i])
    year_sum$PRcor[i] <- cor(mm$GPP, -mm$ER, use = 'complete.obs')
    year_sum$GPP[i] <- mean(mm$GPP, na.rm = T)
    year_sum$GPP.se[i] <- calculate_ts_mean_se(mm$GPP)
    year_sum$ER[i] <- mean(mm$ER, na.rm = T)
    year_sum$ER.se[i] <- calculate_ts_mean_se(mm$ER)
    year_sum$NEP[i] <- mean(mm$NEP, na.rm = T)
    year_sum$NEP.se[i] <- calculate_ts_mean_se(mm$NEP)
    year_sum$PR[i] <- mean(mm$PR, na.rm = T)
    year_sum$PR.se[i] <- calculate_ts_mean_se(mm$PR)
}


summary(year_sum)
mean(year_sum$PR)
#PR_mean_se:
1/nrow(year_sum) * sqrt(sum(year_sum$PR.se ^2))

year_sum %>% mutate(year_sum,
                    PR_low = PR - PR.se - 1,
                    NEP_low = NEP - NEP.se) %>% tibble()

res_table <- year_sum %>%
    mutate(site = substr(site_year, 1,2),
           year = substr(site_year, 4,7)) %>%
    mutate(across(PRcor:PR.se,
                  ~round(., 2))) %>%
    mutate(GPP_mean = paste(GPP, GPP.se, sep = ' $\\pm$ '),
           ER_mean = paste(ER, ER.se, sep = ' $\\pm$ '),
           NEP_mean = paste(NEP, NEP.se, sep = ' $\\pm$ '),
           PR_mean = paste(PR, PR.se, sep =  ' $\\pm$ ')) %>%
    select(site, year, PRcor, ends_with('mean')) %>%
    arrange(site)


latex_table <- xtable::xtable(res_table, caption = "caption")
latex_table
