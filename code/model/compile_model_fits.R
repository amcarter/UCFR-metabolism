#############################################################
# Compare across model fits of GPP to biomass data
#############################################################

library(tidyverse)

mod_metrics <- read_csv('data/model_fits/linear_model_metrics.csv')
mod_metrics <- read_csv('data/model_fits/ar1_linear_model_metrics.csv') %>%
    bind_rows(mod_metrics)
mod_metrics <- read_csv('data/model_fits/PIcurve_model_metrics.csv') %>%
    bind_rows(mod_metrics)
mod_metrics <- read_csv('data/model_fits/ar_PIcurve_model_metrics.csv') %>%
    bind_rows(mod_metrics)

mod_metrics %>%
    arrange(rmse) %>% print(n = 50)
    arrange(loo) %>% filter(model == 'PI_curve')
on_dodge(0.7),alpha=0.5,,adjust=4,width=.5)+ill_viridis(discrete=T,option="G",name="")+ab("ParameterEstimate")+lab("")+theme_bw()

