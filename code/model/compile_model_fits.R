#############################################################
# Compare across model fits of GPP to biomass data
#############################################################

library(tidyverse)

mod_ests <- read_csv('data/model_fits/ar_PIcurve_parameter_ests.csv')
mod_ests <- read_csv('data/model_fits/ar1_linear_model_parameter_ests_ss.csv') %>%
    bind_rows(mod_ests)
chains <- read_csv('data/model_fits/ar_PIcurve_chains.csv')
chains <- read_csv('data/model_fits/ar1_linear_model_chains_ss.csv') %>%
    bind_rows(chains)
mod_metrics <- read_csv('data/model_fits/ar_PIcurve_model_metrics_unconditioned.csv')
mod_metrics <- read_csv('data/model_fits/ar1_linear_model_metrics_ss_uncond.csv') %>%
    bind_rows(mod_metrics)

mod_metrics <- mod_metrics %>%
    mutate(model = c('baseline', rep('linear model', 6),
                     rep('linear model with interaction', 6),
                     rep('linear model interaction only', 6),
                     rep('PI curve model', 6))) %>%
    arrange(model, rmse) %>% print(n = 30)


mod_ests <- mod_ests %>% select(-r2_adj, -rmse, -model) %>%
    rename(covariates = biomass_vars) %>%
    left_join(mod_metrics, by = 'covariates')

chains <- chains %>%
    select(-model) %>%
    rename(covariates = biomass_vars) %>%
    left_join(mod_metrics[,1:2], by = 'covariates')

best_mods_by_cat <- mod_metrics %>%
    group_by(model)%>%
    filter(r2adj > max(r2adj) - 0.012) %>%
    mutate(across(c('r2adj', 'rmse'), round, digits = 2),
           across(starts_with(c('loo', 'waic')), round, digits = 0))

write_csv(best_mods_by_cat, 'data/model_fits/best_models_by_category.csv')

best_mods <- best_mods_by_cat[c(1,2,5,8),]
best_chains <- chains %>%
    right_join(best_mods[,1:2])

best_chains <- data.frame(parameter = c(rep('gamma_fila_chla',2),
                               rep('gamma_epil_chla',2)),
                 chains = rep(0,4),
                 covariates = rep("light",4),
                 model = rep('baseline', 4)) %>%
    bind_rows(best_chains) %>%
    filter(!grepl('^beta', parameter))

interaction_chains
best_chains %>% tibble() %>%
    filter(model == 'linear model with interaction') %>%
    pivot_wider(names_from = parameter, values_from = chains)

best_chains %>%
    filter(grepl('^gamma', parameter)|parameter == 'phi',
           parameter != 'gamma_intercept',
           model != 'PI curve model',
           !is.na(parameter))%>%
    mutate(parameter = case_when(parameter == "gamma_light" ~ 'light',
                                 parameter == "gamma_epil_chla" ~"epil",
                                 parameter == "gamma_fila_chla" ~ "fila",
                                 TRUE ~ parameter),
           parameter = factor(parameter,
                              levels=c("phi", "light",
                                       "fila", "epil"))) %>%
    ggplot(aes(fill = model, y = chains, x = parameter)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ggtitle('GPP=f(bm)')+
    ylab("Parameter Estimate") +
    xlab("") +
    # ylim(-0.1,0.8)+
    # ylim(-0.5,5)+
    theme_bw()


mod_ests %>% filter(parameter == "phi") %>%
    select(mean, sd, covariates, model, r2adj) %>%
    data.frame()
best_chains <- left_join(best_mods[,1:2], chains, by = c('model', 'covariates'))


best_chains %>%
    ggplot(aes(parameter, chains, fill = model)) +
    geom_violin() + facet_wrap(.~parameter, scales = 'free')

a <- mod_ests %>% filter(parameter %in% c("phi",  "sigma")) %>%
    select(parameter, mean, sd, covariates, model) %>%
    right_join(best_mods_by_cat) %>%
    select(-starts_with(c('waic','loo')), -rmse, -r2adj) %>%
    pivot_wider(names_from = 'parameter', values_from = c('mean', 'sd')) %>%
    select(-covariates)%>%
    group_by(model) %>%
    summarize(across(everything(), mean, na.rm = T))
mod_ests %>% filter(grepl('^gamma*', parameter )) %>%
    filter(parameter != 'gamma_intercept')%>%
    select(parameter, mean, sd, covariates, model) %>%
    right_join(best_mods_by_cat) %>%
    select(-starts_with(c('waic','loo')), -rmse, -r2adj) %>%
    pivot_wider(names_from = 'parameter', values_from = c('mean', 'sd')) %>%
    select(-covariates)%>%
    group_by(model) %>%
    summarize(across(everything(), .fn = ~round(mean(., na.rm = T), 2))) %>%
    data.frame()

mod_ests %>% filter(grepl('^Pmax*', parameter )) %>%
    select(parameter, mean, sd, covariates, model) %>%
    right_join(best_mods_by_cat) %>%
    select(-starts_with(c('waic','loo')), -rmse, -r2adj) %>%
    pivot_wider(names_from = 'parameter', values_from = c('mean', 'sd')) %>%
    select(-covariates)%>%
    group_by(model) %>%
    summarize(across(everything(), .fn = ~round(mean(., na.rm = T), 2))) %>%
    data.frame()



write_csv(a, 'data/model_fits/gamma_estimates')

best_ests <- left_join(best_mods_by_cat[,1:2], mod_ests, by = c('model', 'covariates'))

models <- unique(best_ests$model)

# linear model ####
lin_ests <- best_ests %>% filter(model == models[2]) %>%
    slice(1:13) %>%
    data.frame()
lini_ests <- best_ests %>% filter(model == models[4]) %>%
    data.frame()%>%
    slice(1:15)
linio_ests <- best_ests %>% filter(model == models[3]) %>%
    data.frame()
pi_ests <- best_ests %>% filter(model == models[5]) %>%
    data.frame() %>% slice(13:26)

exp(mean(lin_ests$mean[1] * dd$log_GPP, na.rm = T))
exp(mean(lin_ests$mean[3] * dd$light, na.rm = T))
exp(mean(lin_ests$X_2.5[4] * dd$log_epil_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$X_97.5[4] * dd$log_epil_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$mean[4] * dd$log_epil_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$X_2.5[4] * dd$log_epil_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$X_97.5[4] * dd$log_epil_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$mean[5] * dd$log_fila_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$X_2.5[5] * dd$log_fila_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$X_97.5[5] * dd$log_fila_chla_mgm2_fit, na.rm = T))
exp(mean(lin_ests$mean[2])) * 0.00026
dd %>%
    mutate(frac_epil = exp(lin_ests$mean[4]*(1-lin_ests$mean[1]) * log_epil_chla_mgm2_fit),
           frac_fila = exp(lin_ests$mean[5]*(1-lin_ests$mean[1]) * log_fila_chla_mgm2_fit) + frac_epil,
           frac_light = exp(lin_ests$mean[3]*(1-lin_ests$mean[1]) * light) + frac_fila,
           frac_AR = exp(lin_ests$mean[1] * c(2, dd$log_GPP[1:(nrow(dd) - 1)])) + frac_light,
           frac_err = exp(rnorm(nrow(dd), 0, lin_ests$mean[13]*(1-lin_ests$mean[1])))+ frac_AR) %>%
    ggplot(aes(date, GPP)) +
    geom_line() +
    geom_line(aes(y = frac_epil), col = '#1B9EC9', size = 1) +
    geom_line(aes(y = frac_fila), col ='#97BB43', size = 1) +
    geom_line(aes(y = frac_light), col = 'gold', size = 1) +
    geom_line(aes(y = frac_AR), col = 'brown', size = 1) +
    geom_line(aes(y = frac_err), col = 'grey', size = 1) +
    ggtitle('linear model')+
    facet_grid(site~year, scales ='free_x')+
    theme_bw()


dd %>%
    mutate(frac_epil = exp(lini_ests$mean[4]*(1-lini_ests$mean[1]) * log_epil_chla_mgm2_fit +
                               lini_ests$mean[6]*(1-lini_ests$mean[1]) * log_epil_chla_mgm2_fit * light),
           frac_fila = exp(lini_ests$mean[5]*(1-lini_ests$mean[1]) * log_fila_chla_mgm2_fit+
                               lini_ests$mean[7]*(1-lini_ests$mean[1]) * log_fila_chla_mgm2_fit * light) + frac_epil,
           frac_light = exp(lini_ests$mean[3]*(1-lini_ests$mean[1]) * light) + frac_fila,
           frac_AR = exp(lini_ests$mean[1] * c(2, dd$log_GPP[1:(nrow(dd) - 1)])) + frac_light,
           frac_err = exp(rnorm(nrow(dd), 0, lini_ests$mean[15]*(1-lini_ests$mean[1])))+ frac_AR) %>%
    ggplot(aes(date, GPP)) +
    geom_line() +
    geom_line(aes(y = frac_epil), col = '#1B9EC9', size = 1) +
    geom_line(aes(y = frac_fila), col = '#97BB43', size = 1) +
    geom_line(aes(y = frac_light), col = 'gold', size = 1) +
    geom_line(aes(y = frac_AR), col = 'brown', size = 1) +
    geom_line(aes(y = frac_err), col = 'grey', size = 1) +
    ggtitle('linear model with interaction')+
    facet_grid(site~year, scales ='free_x')+
    theme_bw()


dd %>%
    mutate(frac_fila = exp(linio_ests$mean[3] *(1-linio_ests$mean[1])* log_fila_chla_mgm2_fit * light),
           frac_AR = exp(linio_ests$mean[1] * c(2, dd$log_GPP[1:(nrow(dd) - 1)])) + frac_fila,
           frac_err = exp(rnorm(nrow(dd), 0, linio_ests$mean[11]*(1-linio_ests$mean[1])))+ frac_AR) %>%
    ggplot(aes(date, GPP)) +
    geom_line() +
    geom_line(aes(y = frac_fila), col ='#97BB43', size = 1) +
    geom_line(aes(y = frac_AR), col = 'brown', size = 1) +
    geom_line(aes(y = frac_err), col = 'grey', size = 1) +
    ggtitle('linear model interaction only')+
    facet_grid(site~year, scales ='free_x')+
    ylim(0,22)+
    theme_bw()

dd %>%
    mutate(across(starts_with(c('epil', 'fila', 'light')), .fns  = ~. - min(.))) %>%
    mutate(frac_epil = epil_chla_mgm2_fit * pi_ests$mean[3] *(1-pi_ests$mean[5]) *
               tanh(pi_ests$mean[1] * light/pi_ests$mean[3]),
           frac_fila = fila_chla_mgm2_fit * pi_ests$mean[4] *(1-pi_ests$mean[5])*
               tanh(pi_ests$mean[2] * light/pi_ests$mean[4]) + frac_epil,
           frac_AR = pi_ests$mean[5] * c(exp(2), dd$GPP[1:(nrow(dd) - 1)]) + frac_fila,
           frac_err = rnorm(nrow(dd), 0, pi_ests$mean[14]*(1-pi_ests$mean[5]))+ frac_AR) %>%
    select(date, site, year, GPP, starts_with('frac'))  %>%
    # mutate(frac_err = frac_err - frac_AR,
    #        frac_AR = frac_AR - frac_fila,
    #        frac_fila = frac_fila - frac_epil) %>%
    # summarize(across(c(-date, -site, -GPP), .fns = c(mean, sd), na.rm = T))
    ggplot(aes(date, GPP)) +
    geom_line() +
    geom_line(aes(y = frac_epil), col = '#1B9EC9', size = 1 ) +
    geom_line(aes(y = frac_fila), col = '#97BB43', size = 1) +
    geom_line(aes(y = frac_AR), col = 'brown', size = 1) +
    geom_line(aes(y = frac_err), col = 'grey', size = 1) +
    ggtitle('PI_curve_model')+
    facet_grid(site~year, scales ='free_x')+
    ylim(0,22)+
    theme_bw()


mod_metrics %>% filter(model == 'linear model interaction only') %>%
    summarize(across(.cols = c(r2adj, rmse), .fns = c(mean, sd), na.rm = T))
mod_metrics %>% filter(grepl('epil', covariates)& grepl('fila', covariates)) %>%
    summarize(across(.cols = c(r2adj, rmse), .fns = c(mean, sd), na.rm = T))
mod_metrics %>% filter(grepl('chl', covariates)) %>%
    summarize(across(.cols = c(r2adj, rmse), .fns = c(mean, sd), na.rm = T))
mean(dd$GPP)
best_ests$model)
