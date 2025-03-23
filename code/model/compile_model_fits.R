#############################################################
# Compare across model fits of GPP to biomass data
#############################################################

library(tidyverse)
library(kableExtra)
options(scipen = 0)

dd <- read_csv('data/model_fits/biomass_metab_model_data.csv')

dd <- dd %>%
    mutate(year = lubridate::year(date),
           month = lubridate::month(date),
           across(starts_with(c('epil', 'fila')), ~ log(.+ 0.1),
                  .names = 'log_{.col}'),
           across(starts_with(c('GPP')), ~ log(.), .names = 'log_{.col}'),
           GPP_sd = (GPP.upper - GPP.lower)/3.92,
           log_GPP_sd = (log_GPP.upper - log_GPP.lower)/3.92) %>%
    mutate(across(contains(c('epil','fila', 'light')), ~ scale(.)[,1])) %>%
    group_by(site, year) %>%
    mutate(
        across(where(is.numeric), zoo::na.approx, na.rm = FALSE),
        site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN'))) %>%
    filter(!is.na(GPP), !is.na(epil_gm2_fit), !is.na(fila_gm2_fit),
           !is.na(light)) %>%
    ungroup()

dd <- dd %>% select(site, GPP, log_GPP, light, epil_chla_mgm2_fit, log_epil_chla_mgm2_fit,
              fila_chla_mgm2_fit, log_fila_chla_mgm2_fit)


mod_ests <- read_csv('data/model_fits/ar_err_log_model_parameter_ests_ss.csv')
chains <- read_csv('data/model_fits/ar_err_log_model_chains_ss.csv')
mod_metrics <- read_csv('data/model_fits/ar_err_log_model_metrics_ss_cond.csv')
mod_holdouts <- read_csv('data/model_fits/ar_err_log_model_holdout_RMSEs.csv')
mod_holdout_preds <- read_csv('data/model_fits/ar_err_log_model_holdout_predictions.csv')

mod_rmse <- mod_holdout_preds %>%
    group_by(biomass_vars) %>%
    mutate(estim = exp(estim),
           sqe = (GPP - estim)^2) %>%
    summarize(rmspe = sqrt(mean(sqe)))

mod_ests <- left_join(mod_ests, mod_rmse, by = 'biomass_vars')

mod_metrics <- mod_metrics %>%
    mutate(model = c('baseline', rep('linear model', 6),
                     rep('linear model with interaction', 6),
                     rep('linear model interaction only', 6),
                     rep('linear model', 2),
                     rep('linear model with interaction', 2),
                     rep('linear model interaction only', 2))) %>%
    arrange(model, rmse)


chains <- chains %>%
    select(-model) %>%
    rename(covariates = biomass_vars) %>%
    left_join(mod_metrics[,1:2], by = 'covariates')

best_mods <- mod_metrics %>%
    left_join(rename(mod_rmse, covariates = biomass_vars)) %>%
    group_by(model) %>% arrange(rmspe) %>%
    slice(1)
best_mods <- mod_metrics %>%
    left_join(rename(mod_rmse, covariates = biomass_vars)) %>%
    group_by(model) %>% arrange(rmse) %>%
    slice(1) %>% ungroup() %>% slice(-1) %>%
    bind_rows(best_mods) %>% arrange(model, rmspe) %>%
    select(-rmse,-starts_with(c('waic', 'loo')))

write_csv(best_mods, 'data/model_fits/best_models_by_category.csv')


best_mods %>%
    mutate(RMSE = round(rmspe, 2),
           r2 = round(r2, 2)) %>%
    select(-rmspe) %>%
    arrange(model, RMSE, desc = TRUE) %>%
    relocate(r2, .after = RMSE) %>%
    kbl(format = 'latex',
        linesep = '') %>%
    kable_classic(full_width = T, html_font = 'helvetica')

mod_ests %>%
    filter(grepl('chla', biomass_vars)| biomass_vars == 'log_light') %>%
    rename(covariates = biomass_vars,
           r2 = r2_adj, RMSE = rmspe) %>%
    select(-model, -rmse) %>%
    left_join(select(mod_metrics, covariates, model)) %>%
    mutate(value = case_when(mean >= 1 ~ paste0(signif(mean, 2), ' (', round(sd, 1), ')'),
                             mean < 1 ~ paste0(round(mean, 2), ' (', round(sd, 2), ')')),
           bio_frac = case_when(grepl('biomass', covariates) ~ 'sum',
                                grepl('epil.*fila', covariates) ~ 'both',
                                grepl('epil', covariates) ~ 'epil',
                                grepl('fila', covariates) ~ 'fila',
                                TRUE ~ 'none')) %>%
    select(model, covariates, bio_frac, parameter, value, r2, RMSE) %>%
    filter(!grepl('^beta', parameter)) %>%
    mutate(parameter = case_when(grepl('chla$', parameter) ~ substr(parameter, 1, nchar(parameter)-5),
                                 TRUE ~ parameter),
           parameter = case_when(grepl('^gamma_log', parameter) ~ substr(parameter, 11, nchar(parameter)),
                                 grepl('^gamma', parameter) ~ substr(parameter, 7, nchar(parameter)),
                                 TRUE ~ parameter),
           parameter = case_when(bio_frac == parameter ~ 'biomass',
                                 parameter == paste0('light_', bio_frac) ~ 'light_biomass',
                                 TRUE ~ parameter)) %>%
    pivot_wider(names_from = parameter, values_from = value) %>%
    mutate(epil = if_else(!is.na(epil), paste(epil, fila, sep = '\n'), NA_character_),
           biomass = if_else(is.na(biomass), epil, biomass),
           light_epil = if_else(!is.na(light_epil), paste(light_epil, light_fila, sep = '\n'),
                                NA_character_),
           light_biomass = if_else(is.na(light_biomass), light_epil, light_biomass)) %>%
    select(-epil, -light_epil, -fila, -light_fila, -intercept, -tau) %>%
    relocate(model, phi, sigma, light, biomass,
             light_biomass, r2, RMSE) %>%
    mutate(bio_frac = factor(bio_frac, levels = c('none', 'epil', 'fila',
                                                  'both', 'sum')),
           model = case_when(model == 'baseline' ~ 0,
                             model == 'linear model' ~ 1,
                             model == 'linear model interaction only' ~ 2,
                             model == 'linear model with interaction' ~ 3)) %>%
    arrange(bio_frac, model) %>%
    relocate(bio_frac)%>%
    select(-covariates) %>%
    mutate(r2 = round(r2, 3),
           RMSE = round(RMSE, 2)) %>%
    kbl(format = 'latex',
        linesep = '') %>%
    kable_classic( html_font = 'helvetica')

mod_ests %>%
    filter(grepl('afdm', biomass_vars)| biomass_vars == 'log_light') %>%
    rename(covariates = biomass_vars,
           r2 = r2_adj, RMSE = rmspe) %>%
    select(-model, -rmse) %>%
    left_join(select(mod_metrics, covariates, model)) %>%
    mutate(value = case_when(mean >= 1 ~ paste0(signif(mean, 2), ' (', round(sd, 1), ')'),
                             mean < 1 ~ paste0(round(mean, 2), ' (', round(sd, 2), ')')),
           bio_frac = case_when(grepl('biomass', covariates) ~ 'sum',
                                grepl('epil.*fila', covariates) ~ 'both',
                                grepl('epil', covariates) ~ 'epil',
                                grepl('fila', covariates) ~ 'fila',
                                TRUE ~ 'none')) %>%
    select(model, covariates, bio_frac, parameter, value, r2, RMSE) %>%
    filter(!grepl('^beta', parameter)) %>%
    mutate(parameter = case_when(grepl('afdm$', parameter) ~ substr(parameter, 1, nchar(parameter)-5),
                                 TRUE ~ parameter),
           parameter = case_when(grepl('^gamma_log', parameter) ~ substr(parameter, 11, nchar(parameter)),
                                 grepl('^gamma', parameter) ~ substr(parameter, 7, nchar(parameter)),
                                 TRUE ~ parameter),
           parameter = case_when(bio_frac == parameter ~ 'biomass',
                                 parameter == paste0('light_', bio_frac) ~ 'light_biomass',
                                 TRUE ~ parameter)) %>%
    pivot_wider(names_from = parameter, values_from = value) %>%
    mutate(epil = if_else(!is.na(epil), paste(epil, fila, sep = '\n'), NA_character_),
           biomass = if_else(is.na(biomass), epil, biomass),
           light_epil = if_else(!is.na(light_epil), paste(light_epil, light_fila, sep = '\n'),
                                NA_character_),
           light_biomass = if_else(is.na(light_biomass), light_epil, light_biomass)) %>%
    select(-epil, -light_epil, -fila, -light_fila, -intercept, -tau) %>%
    relocate(model, phi, sigma, light, biomass,
             light_biomass, r2, RMSE) %>%
    mutate(bio_frac = factor(bio_frac, levels = c('none', 'epil', 'fila',
                                                  'both', 'sum')),
           model = case_when(model == 'baseline' ~ 0,
                             model == 'linear model' ~ 1,
                             model == 'linear model interaction only' ~ 2,
                             model == 'linear model with interaction' ~ 3)) %>%
    arrange(bio_frac, model) %>%
    relocate(bio_frac)%>%
    select(-covariates) %>%
    mutate(r2 = round(r2, 3),
           RMSE = round(RMSE, 2)) %>%
    kbl(format = 'latex',
        linesep = '') %>%
    kable_classic( html_font = 'helvetica')

mutate(afdm, units = 'afdm') %>%
    bind_rows(mutate(chla, units = 'chla')) %>%
    pivot_wider(id_cols = c('bio_frac', 'model'),
                names_from = 'units', values_from = 'r2') %>%
    mutate(diff = (chla - afdm)/afdm) %>% slice(-13) %>%
    group_by(bio_frac) %>% summarize(diff = mean(diff))

# format model table for SI:
mod_ests_table <- mod_ests %>%
    rename(covariates = biomass_vars,
           r2 = r2_adj, RMSPE = rmspe) %>%
    select(-model) %>%
    left_join(select(mod_metrics, covariates, model) )%>%
    mutate(value = case_when(mean >= 1 ~ paste0(signif(mean, 2), ' (', round(sd, 1), ')'),
                             mean < 1 ~ paste0(round(mean, 2), ' (', round(sd, 2), ')')),
           units = case_when(grepl('chla$', covariates) ~ 'chla',
                             grepl('afdm$', covariates) ~ 'mass',
                             TRUE ~ NA_character_)) %>%
    select(model, covariates, units, parameter, value, r2, RMSPE) %>%
    filter(!grepl('^beta', parameter)) %>%
    mutate(parameter = case_when(grepl('chla$', parameter) ~ substr(parameter, 1, nchar(parameter)-5),
                                 grepl('afdm$', parameter) ~ substr(parameter, 1, nchar(parameter)-5),
                                 TRUE ~ parameter),
           parameter = case_when(grepl('^gamma_log', parameter) ~ substr(parameter, 11, nchar(parameter)),
                                 grepl('^gamma', parameter) ~ substr(parameter, 7, nchar(parameter)),
                                 TRUE ~ parameter)) %>%
    pivot_wider(names_from = parameter, values_from = value) %>%
    relocate(model, theta, sigma, intercept, tau, light, biomass, fila, epil,
             light_biomass, light_fila, light_epil, r2, RMSPE) %>%
    select(-covariates)

mod_ests_table %>%
    mutate(across(-any_of(c('model', 'r2', 'RMSPE')),
                  function(x) x = case_when(is.na(x) ~ '-',
                                            TRUE ~ x)),
           r2 = round(r2, 2), RMSPE = round(RMSPE, 2)) %>%
     kbl(format = 'latex',
         linesep = '') %>%
        kable_classic(full_width = F, html_font = 'helvetica')


# best_mods <- best_mods_by_cat[c(1,2,4,5,8),]
best_chains <- chains %>%
    right_join(best_mods[,1:2])

best_chains <- data.frame(parameter = c(rep('gamma_fila_chla',2),
                               rep('gamma_epil_chla',2)),
                 chains = rep(0,4),
                 covariates = rep("light",4),
                 model = rep('baseline', 4)) %>%
    bind_rows(best_chains) %>%
    filter(!grepl('^beta', parameter))

best_coeffs <- mod_ests %>% select(-r2_adj, -rmse, -model) %>%
    rename(covariates = biomass_vars) %>%
    left_join(best_mods, by = 'covariates')

mod_ests <- mod_ests %>% select(-r2_adj, -rmse, -model) %>%
    rename(covariates = biomass_vars) %>%
    left_join(mod_metrics, by = 'covariates')

baseline <- filter(best_chains, covariates == 'light' & model == 'baseline')
linear <- filter(best_chains, model == 'linear model')
linear_i <- filter(best_chains, model == 'linear model with interaction')

linear %>% group_by(parameter) %>% summarize(mean = mean(chains))

phi <- baseline$chains[baseline$parameter == 'phi']
g_light <- baseline$chains[baseline$parameter == 'gamma_light']
mean(phi)
mean(g_light * mean(dd$light))/mean(dd$log_GPP)

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

best_ests <- left_join(best_mods_by_cat[c(1,2,4,5,8),1:2],
                       mod_ests, by = c('model', 'covariates'))

models <- unique(best_ests$model)

# linear model ####
baseline_ests <- best_ests %>% filter(model == models[1]) %>%
    data.frame()
lin_ests <- best_ests %>% filter(model == models[2]) %>%
    data.frame()
lini_ests <- best_ests %>% filter(model == models[4]) %>%
    data.frame()
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

lin_ests

dd %>%
    mutate(frac_epil = exp(lin_ests$mean[4]*(1-lin_ests$mean[1]) * log_epil_chla_mgm2_fit),
           frac_fila = exp(lin_ests$mean[5]*(1-lin_ests$mean[1]) * log_fila_chla_mgm2_fit) + frac_epil,
           frac_light = exp(lin_ests$mean[3]*(1-lin_ests$mean[1]) * light) + frac_fila,
           frac_AR = exp(lin_ests$mean[1] * c(dd$log_GPP[1], dd$log_GPP[1:(nrow(dd) - 1)])) + frac_light,
           frac_err = rnorm(nrow(dd), 0, lin_ests$mean[13]*(1-lin_ests$mean[1]))+ frac_AR,
           GPP = GPP - exp(lin_ests$mean[2]*(1-lin_ests$mean[1]))) %>%
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

baseline_chains <- filter(best_chains, model == models[1])
lin_chains <- filter(best_chains, model == models[2])
linio_chains <- filter(best_chains, model == models[3])
lini_chains <- filter(best_chains, model == models[4])
pi_chains <- filter(best_chains, model == models[5])


hist(baseline_chains$chains[baseline_chains$parameter == 'phi'] +
         baseline_chains$chains[baseline_chains$parameter == 'gamma_light']
         )
baseline_chains  <- filter(baseline_chains, parameter == 'phi') %>%
    mutate(frac_AR = chains,
           frac_light = 1-chains) %>%
    dplyr::select(starts_with('frac'), model)

lin <- filter(lin_chains, parameter == 'phi') %>%
    select(model, frac_AR = chains)

light <- lin_chains$chains[lin_chains$parameter == 'gamma_light']
fila <- lin_chains$chains[lin_chains$parameter == 'gamma_fila_chla']
epil <- lin_chains$chains[lin_chains$parameter == 'gamma_epil_chla']
lin$frac_light = (1-lin$frac_AR) * light/(light + fila+ epil)
lin$frac_epil = (1-lin$frac_AR) * epil/(light + fila+ epil)
lin$frac_fila = (1-lin$frac_AR) * fila/(light + fila+ epil)


lini <- filter(lini_chains, parameter == 'phi') %>%
    select(model, frac_AR = chains)

light <- lini_chains$chains[lini_chains$parameter == 'gamma_light']
fila <- lini_chains$chains[lini_chains$parameter == 'gamma_fila_chla'] +
    lini_chains$chains[lini_chains$parameter == 'gamma_light_fila_chla']
epil <- lini_chains$chains[lini_chains$parameter == 'gamma_epil_chla']+
    lini_chains$chains[lini_chains$parameter == 'gamma_light_epil_chla']
lini$frac_light = (1-lini$frac_AR) * light/(light + fila+ epil)
lini$frac_epil = (1-lini$frac_AR) * epil/(light + fila+ epil)
lini$frac_fila = (1-lini$frac_AR) * fila/(light + fila+ epil)

linio <- filter(linio_chains, parameter == 'phi') %>%
    select(model, frac_AR = chains)

fila <- linio_chains$chains[linio_chains$parameter == 'gamma_light_fila_chla']
linio$frac_epil = (1-linio$frac_AR) * epil/(fila+ epil)
linio$frac_fila = (1-linio$frac_AR)

pic <- filter(pi_chains, parameter == 'phi') %>%
    select(model, frac_AR = chains)

epil = mean(dd$epil_chla_mgm2_fit - min(dd$epil_chla_mgm2_fit)) *
    pi_chains$chains[pi_chains$parameter == 'Pmax_epil_chla'] *
    (1-pi_chains$chains[pi_chains$parameter == 'phi']) *
    tanh(pi_chains$chains[pi_chains$parameter == 'alpha_epil_chla'] *
             mean(dd$light - min(dd$light))/
             pi_chains$chains[pi_chains$parameter == 'Pmax_epil_chla'])
fila = mean(dd$fila_chla_mgm2_fit - min(dd$fila_chla_mgm2_fit)) *
    pi_chains$chains[pi_chains$parameter == 'Pmax_fila_chla'] *
    (1-pi_chains$chains[pi_chains$parameter == 'phi']) *
    tanh(pi_chains$chains[pi_chains$parameter == 'alpha_fila_chla'] *
             mean(dd$light - min(dd$light))/
             pi_chains$chains[pi_chains$parameter == 'Pmax_fila_chla'])

pic$frac_epil = (1-pic$frac_AR) * epil/(fila+ epil)
pic$frac_fila = (1-pic$frac_AR) * fila/(fila+ epil)




bind_rows(baseline_chains, lin, lini, linio, pic) %>%
    pivot_longer(cols = starts_with('frac'), names_to = 'component',
                 names_pattern = 'frac_([A-Za-z]+)', values_to = 'fraction') %>%
    mutate(model = factor(model, levels = c('baseline', 'linear model', 'linear model with interaction', 'linear model interaction only', 'PI curve model')),
           component = factor(component, levels = c('AR', 'light', 'epil', 'fila'))) %>%
    ggplot(aes(fill = model, x = component, y = fraction)) +
    geom_violin(position=position_dodge(0.7), alpha=0.5,
                scale = 'width', adjust =4, width = .5) +
    # facet_grid(biomass_vars~model)+
    scale_fill_viridis(discrete=T, option = "G", name="") +
    ggtitle('')+
    ylab("Fraction of modeled GPP") +
    xlab("") +
    # ylim(-0.1,0.8)+
    # ylim(-0.5,5)+
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
    mutate(frac_err = frac_err - frac_AR,
           frac_AR = frac_AR - frac_fila,
           frac_fila = frac_fila - frac_epil) %>%
    summarize(across(c(-date, -site, -GPP), .fns = c(mean, sd), na.rm = T))
