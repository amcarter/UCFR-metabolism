# fit GAMS to biomass data
# 10/2022

library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(posterior)


# read in datasets:
met <- read_csv('data/metabolism/metabolism_compiled_all_sites_mle_fixedK_correctedSE.csv') %>%
  mutate(site = case_when(site == 'BM' ~ 'BG',
                          TRUE ~ site),
         site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')),
         year = factor(year),
         GPP.se = (GPP.upper - GPP.lower)/(2*1.96),
         ER.se = (ER.upper - ER.lower)/(2*1.96))

q90 <- read_csv('data/quantile_PR_fits_summary_brms.csv') %>%
  mutate(site = case_when(site == 'BM' ~ 'BG',
                          TRUE ~ site),
         year = factor(year))
met <- left_join(met, select(q90, site, year, ARf = slope, ARf.se = slope.se),
                 by = c('site', 'year')) %>%
  mutate(year = factor(year),
         NPP = GPP * (1-ARf),
         NPP.se = sqrt(GPP.se ^2 + ARf.se^2),
         # NPP_globalARf = GPP * (1 + beta[2,1]),
         AR = -GPP *(ARf),
         AR.se = sqrt(GPP.se^2 + ARf.se^2),
         HR = ER - AR,
         HR.se = sqrt(ER.se^2 + AR.se^2)) %>%
  select(-msgs.fit, -warnings, -errors, -K600, -DO_fit)

# need to get light into the correct units! ####
light <- read_csv('data/site_data/daily_modeled_light_all_sites.csv') %>%
  mutate(site = case_when(site == 'BM' ~ 'BG',
                          TRUE ~ site))

biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    filter(!is.na(site) & site != 'CR') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN')),
           year = as.factor(lubridate::year(date)),
           site_year = paste(site, year, sep = '_'))

glimpse(biomass)


# look at the data
# biomass %>%
#     select(date, doy, year, site, ends_with('_gm2')) %>%
#     pivot_longer(cols = -c('date', 'doy', 'year', 'site'),
#                  names_to = 'bm_category',
#                  values_to = 'gm2') %>%
#     ggplot(aes(doy, gm2, col = year)) +
#     geom_point() +
#     facet_grid(site~bm_category, scales = 'free_y')


# find 2%quantile ####
# add the 2% quantile to each of the biomass categories in order to
# fit log transformed GAMs
colnames(biomass)
min_vals <- biomass %>%
    select(epilithon_gm2, epilithon_chla_mgm2,
           filamentous_gm2, filamentous_chla_mgm2,
           fila_macro_gm2) %>%
    mutate(across(.fns = ~case_when(. == 0 ~ NA_real_,
                                    TRUE ~ .))) %>%
    pivot_longer(cols = everything(), names_to = 'biomass', values_to = 'value') %>%
    group_by(biomass) %>%
    summarize(min = quantile(value, 0.01, na.rm = T),
              q2 = quantile(value, 0.02, na.rm = T))

dates <- c(seq(min(biomass$date), max(biomass$date[biomass$year == 2020]), by = 'day'),
           seq(min(biomass$date[biomass$year == 2021]), max(biomass$date), by = 'day'))


# biomass <- biomass %>%
#     mutate(epilithon_gm2 = epilithon_gm2 + min_vals$min[min_vals$biomass == 'epilithon_gm2'],
#            epilithon_chla_mgm2 = epilithon_chla_mgm2 +
#                min_vals$min[min_vals$biomass == 'epilithon_chla_mgm2'],
#            filamentous_gm2 = filamentous_gm2 + min_vals$min[min_vals$biomass == 'filamentous_gm2'],
#            filamentous_chla_mgm2 = filamentous_chla_mgm2 +
#                min_vals$min[min_vals$biomass == 'filamentous_chla_mgm2'],
#            fila_macro_gm2 = fila_macro_gm2 + min_vals$min[min_vals$biomass == 'fila_macro_gm2'])
# fit gams
preds <- data.frame(date = dates,
                    doy = as.numeric(format(dates, '%j')),
                    year = lubridate::year(dates))%>%
    mutate(year = as.factor(year))
s_preds <- data.frame()
for(site in unique(biomass$site)){
    preds$site = site
    s_preds <- bind_rows(s_preds, preds)
}
biomass_sum <- select(biomass,
       date, year, site, filamentous_gm2, epilithon_gm2,
       filamentous_chla_mgm2, epilithon_chla_mgm2) %>%
    group_by(date, year, site) %>%
    summarize(across(filamentous_gm2:epilithon_chla_mgm2,
              ~mean(.x, na.rm = T)))

summary(biomass_sum)
min_vals

bm_met <- s_preds %>%
  tibble() %>%
  full_join(biomass,
            by = c("site", "date", "year", "doy")) %>%
  mutate(
    epilithon_gm2_min  = epilithon_gm2 + min_vals$min[min_vals$biomass == "epilithon_gm2"],
    filamentous_gm2_min = filamentous_gm2 + min_vals$min[min_vals$biomass == 'filamentous_gm2'],
    epilithon_chla_mgm2_min  = epilithon_chla_mgm2 + min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"],
    filamentous_chla_mgm2_min = filamentous_chla_mgm2 + min_vals$min[min_vals$biomass == 'filamentous_chla_mgm2'],
    site_year = factor(paste(site, year, sep = "_")),
    site = factor(site),
    year = factor(year)) %>%
  left_join(select(ungroup(met), site, date, GPP, ER, ARf, NPP,
                   GPP.se, ER.se, ARf.se, NPP.se),
            by = c('site', 'date')) %>%
  left_join(select(light, site, date, PAR_bc_Jm2)) %>%
  mutate(light = PAR_bc_Jm2/max(PAR_bc_Jm2),
         across(ends_with("_min"), \(x) if_else(is.na(x), Inf, x),
                .names = '{.col}2'),
         NPP2 = if_else(is.na(NPP), Inf, NPP))
  # filter(!is.na(GPP))
glimpse(bm_met)
sites <- unique(bm_met$site)
bm_sum2 <- bm_met %>%
  mutate(site = factor(site, levels = sites)) %>%
  group_by(site, date, doy) %>%
  summarize(across(epilithon_gm2_min:filamentous_chla_mgm2_min,
                   list(mean = ~mean(.x, na.rm = T),
                        se = ~sd(.x, na.rm = T))),
            NPP = mean(NPP, na.rm = T)) %>%
  ungroup() %>%
  mutate(NPP_chla = NPP,
         year = lubridate::year(date),
         site_year = paste(site, year, sep = "_"))

bform <- bf(epilithon_gm2_min | mi() ~ s(doy) +
              s(doy, site_year, bs = 'fs'))

# f2_bgamma <- brm(bform,
#       data = bm_met,
#       family = Gamma(link = 'log'),
#       chains = 4, cores = 4,
#       iter = 1000)
# saveRDS(f2_bgamma, "data/model_fits/brms_gam_for_spline_full.rds")
f2_bgamma <- readRDS("data/model_fits/brms_gam_for_spline_full.rds")

# stancode(f2_bgamma)
str(f2_bgamma)

DOY_spline <- standata(f2_bgamma)$Xs
bf2 <- bf(NPP|mi() ~ 0 + epilithon_gm2 + filamentous_gm2)
# stancode(bf2, data = bm_met)

# stan_mod <- stan_model("code/model/stan_code/full_biomass_NPP_mod.stan")
stan_mod_err <- stan_model("code/model/stan_code/full_biomass_NPPerr_mod.stan")

# fit single algal biomass gam:
# dat_list <- list(
#   N = nrow(bm_met),
#   epil_meas = bm_met$epilithon_gm2_min2,
#   Nmi_epil = sum(is.na(bm_met$epilithon_gm2_min)),
#   Jmi_epil = which(is.na(bm_met$epilithon_gm2_min)),
#   Ks = 1,
#   Xs = DOY_spline,
#   nb_1 = 1,
#   knots_1 = as.array(8),
#   Zs_1_1 = standata(f2_bgamma)$Zs_1_1,
#   nb_2 = 3,
#   knots_2 = as.array(c(96, 12, 12)),
#   Zs_2_1 = standata(f2_bgamma)$Zs_2_1,
#   Zs_2_2 = standata(f2_bgamma)$Zs_2_2,
#   Zs_2_3 = standata(f2_bgamma)$Zs_2_3
# )
#
# fite <- sampling(stan_mod,
#                  data = dat_list,
#                  iter = 1000,
#                  chains = 4,
#                  cores = 4)
# fit model on AFDM biomass: ####
dat_list <- list(
  N = nrow(bm_met),
  epil_meas = bm_met$epilithon_gm2_min2,
  fila_meas = bm_met$filamentous_gm2_min2,
  NPP_est = bm_met$NPP2,
  NPP_se = mean(bm_met$NPP.se, na.rm = T),
  # light = bm_met$PAR_bc_Jm2,
  light = bm_met$light,
  Nmi_epil = sum(is.na(bm_met$epilithon_gm2_min)),
  Nmi_fila = sum(is.na(bm_met$filamentous_gm2_min)),
  Nmi_NPP = sum(is.na(bm_met$NPP)),
  Jmi_epil = which(is.na(bm_met$epilithon_gm2_min)),
  Jmi_fila = which(is.na(bm_met$filamentous_gm2_min)),
  Jmi_NPP = which(is.na(bm_met$NPP)),
  K = 2,
  # data for fitting splines to DOY:
  Ks = 1,
  Xs = DOY_spline,
  nb_1 = 1,
  knots_1 = as.array(8),
  Zs_1_1 = standata(f2_bgamma)$Zs_1_1,
  # DOY splines grouped by siteyear
  nb_2 = 3,
  knots_2 = as.array(c(96, 12, 12)),
  Zs_2_1 = standata(f2_bgamma)$Zs_2_1,
  Zs_2_2 = standata(f2_bgamma)$Zs_2_2,
  Zs_2_3 = standata(f2_bgamma)$Zs_2_3
)

fitg <- sampling(stan_mod_err,
                 data = dat_list,
                 control = list(max_treedepth = 15,
                                adapt_delta = 0.9),
                 #                stepsize = 0.1),
                 iter = 4000,
                 chains = 4,
                 cores = 4)

saveRDS(fitg, "full_biomass_partition_NPPerr_mod_td15_ad9_ss1_iter4000.rds")


# fit model on chla biomass:####
dat_list <- list(
  N = nrow(bm_met),
  epil_meas = bm_met$epilithon_chla_mgm2_min2,
  fila_meas = bm_met$filamentous_chla_mgm2_min2,
  NPP_est = bm_met$NPP2,
  NPP_se = mean(bm_met$NPP.se, na.rm = T),
  # light = bm_met$PAR_bc_Jm2,
  light = bm_met$light,
  Nmi_epil = sum(is.na(bm_met$epilithon_chla_mgm2_min)),
  Nmi_fila = sum(is.na(bm_met$filamentous_chla_mgm2_min)),
  Nmi_NPP = sum(is.na(bm_met$NPP)),
  Jmi_epil = which(is.na(bm_met$epilithon_chla_mgm2_min)),
  Jmi_fila = which(is.na(bm_met$filamentous_chla_mgm2_min)),
  Jmi_NPP = which(is.na(bm_met$NPP)),
  K = 2,
  # data for fitting splines to DOY:
  Ks = 1,
  Xs = DOY_spline,
  nb_1 = 1,
  knots_1 = as.array(8),
  Zs_1_1 = standata(f2_bgamma)$Zs_1_1,
  # DOY splines grouped by siteyear
  nb_2 = 3,
  knots_2 = as.array(c(96, 12, 12)),
  Zs_2_1 = standata(f2_bgamma)$Zs_2_1,
  Zs_2_2 = standata(f2_bgamma)$Zs_2_2,
  Zs_2_3 = standata(f2_bgamma)$Zs_2_3
)

fitchla <- sampling(stan_mod_err,
                      data = dat_list,
                      control = list(max_treedepth = 15,
                                     adapt_delta = 0.95,
                                     stepsize = 0.1),
                      iter = 4000,
                      chains = 4,
                      cores = 4)

saveRDS(fitchla, "full_biomass_partition_chla_NPPerr_mod_td15_ad95_ss01_iter4000.rds")

# Evaluate model fits: ####
#AFDM model
fitg <- readRDS("full_biomass_partition_NPPerr_mod_td15_ad9_ss1_iter4000.rds")
traceplot(fitg, pars = c("b", "sigma", "lp__"))
pairs(fitg, pars = c("b", "sigma", "lp__"))

posterior_draws <- as_draws_df(fitg, .nchains = 4)
summary_stats <- summarize_draws(posterior_draws)
summary_stats$percent_ess <- summary_stats$ess_tail/20000
summary(summary_stats)
print(summary_stats[, c("variable", "rhat", "ess_bulk", "ess_tail")])

# Check for divergences
sampler_params <- get_sampler_params(fitg, inc_warmup = FALSE)
n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
cat("Number of divergent transitions:", n_divergent, "\n")


# Chla model
fitchla <- readRDS("full_biomass_partition_chla_NPPerr_mod_td15_ad95_ss01_iter4000.rds")
traceplot(fitchla, pars = c("b", "sigma", "lp__"))
pairs(fitchla, pars = c("b", "sigma", "lp__"))

posterior_draws <- as_draws_df(fitchla, .nchains = 4)
summary_stats <- summarize_draws(posterior_draws)
summary_stats$percent_ess <- summary_stats$ess_tail/4000
summary(summary_stats)
print(summary_stats[, c("variable", "rhat", "ess_bulk", "ess_tail")])

# Check for divergences
sampler_params <- get_sampler_params(fitchla, inc_warmup = FALSE)
n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
cat("Number of divergent transitions:", n_divergent, "\n")


# 2. Posterior predictive checks ----

# Simulated posterior predictive data
# Assume you generated `y_rep` inside Stan and it is accessible
y_rep <- posterior::as_draws_df(fitchla, variable = "y_rep") %>% select(starts_with("NPP["))  %>% as.matrix()
y_rep_g <- posterior::as_draws_df(fitg, variable = "y_rep") %>% select(starts_with("NPP["))  %>% as.matrix()
# Your observed data
bm_met$NPP[1] <- 5
y_obs <- bm_met$NPP # vector of observed responses
# Plot 1: Distributional overlap
p_chla <- ppc_dens_overlay(y = zoo::na.approx(y_obs), yrep = y_rep[1:50, ])+ # overlay first 100 posterior draws
  xlab("NPP (g C m-2 d-1)") +
  ggtitle('Model with chl a biomass')
p_g <- ppc_dens_overlay(y = zoo::na.approx(y_obs), yrep = y_rep[1:50, ])+ # overlay first 100 posterior draws
  xlab("NPP (g C m-2 d-1)") +
  ggtitle('Model with AFDM biomass')

png("figures/ppcheck_NPP_full_biomass_models.png", width = 3.6, height = 5, units = "in", res = 300)
  ggpubr::ggarrange(p_g, p_chla, ncol = 1)
dev.off()

# 90% credible intervals
y_pred_q <- apply(y_rep, 2, quantile, probs = c(0.05, 0.95))
y_pred_qg <- apply(y_rep_g, 2, quantile, probs = c(0.05, 0.95))
coverage <- (y_obs > y_pred_q[1, ]) & (y_obs < y_pred_q[2, ])
coverageg <- (y_obs > y_pred_qg[1, ]) & (y_obs < y_pred_qg[2, ])
coverage_rate_chla <- mean(coverage, na.rm = T)
coverage_rate_AFDM <- mean(coverageg, na.rm = T)

# Bayesian R^2
# Assume y_rep is matrix (iterations x observations)
var_fit <- apply(y_rep, 1, var) # variance of predictions
var_fitg <- apply(y_rep_g, 1, var) # variance of predictions
var_resid <- apply((y_rep - y_obs)^2, 1, function(x) mean(x, na.rm = T)) # variance of residuals
var_residg <- apply((y_rep_g - y_obs)^2, 1, function(x) mean(x, na.rm = T)) # variance of residuals
R2_bayes_chla <- mean(var_fit / (var_fit + var_resid))
R2_bayes_AFDM <- mean(var_fitg / (var_fitg + var_residg))

labs_dat <- data.frame(model = c("Chl a Model", "AFDM Model"),
           R2_bayes = c(R2_bayes_chla, R2_bayes_AFDM),
           coverage = c(coverage_rate_chla, coverage_rate_AFDM),
           obs = c(0.1,0.1), pred = c(12, 12)) %>%
  mutate(lab = paste0("Bayesian R² (posterior mean):", round(R2_bayes, 2), "\n",
                      "Empirical 90% coverage rate:", round(coverage, 3) * 100, "%"))
# Plot 2: Mean vs predicted mean
y_pred_mean <- colMeans(y_rep)
y_pred_g_mean <- colMeans(y_rep_g)
png("figures/mean_vs_predicted_fullbiomassmods.png", width = 7, height = 3.5, units = "in", res = 300)
  ggplot(data.frame(obs = c(y_obs, y_obs), pred = c(y_pred_mean, y_pred_g_mean),
                    model = rep(c("Chl a Model", "AFDM Model"), each = length(y_obs))),
         aes(x = obs, y = pred)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = "grey25") +
    labs(x = "NPP estimates", y = "Posterior Predicted Mean") +
    geom_text(data = labs_dat, aes(label = lab), adj = c(0,0))+
    facet_grid(.~model) + ylim(0, 12.5) +
    theme_bw()
dev.off()



# Extract fits: ####
# AFDM fits
draws_afdm <- as_draws_df(fitg)
coefs_afdm <- draws_afdm %>% select(starts_with(c("b[", "sigma")))
apply(coefs_afdm, 2, function(x) {
  c(median = median(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
})
# chla fits
draws_chla <- as_draws_df(fitchla)
coefs_chla <- draws_chla %>% select(starts_with(c("b[", "sigma")))
apply(coefs_chla, 2, function(x) {
  c(median = median(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
})

# extract missing data estimates
Ymi_epil <- draws_afdm %>% select(starts_with("Ymi_epil"))
Ymi_fila <- draws_afdm %>% select(starts_with("Ymi_fila"))
Ymi_NPP <- draws_afdm %>% select(starts_with("NPP["))
Ymi_epil_chla <- draws_chla %>% select(starts_with("Ymi_epil"))
Ymi_fila_chla <- draws_chla %>% select(starts_with("Ymi_fila"))
Ymi_NPP_chla <- draws_chla %>% select(starts_with("NPP"))
Ymi_sum <- data.frame(
  epilithon_gm2 = apply(Ymi_epil, 2, median),
  epilithon_gm2_lower = apply(Ymi_epil, 2, function(x) quantile(x, 0.025)),
  epilithon_gm2_upper = apply(Ymi_epil, 2, function(x) quantile(x, 0.975)),
  epilithon_gm2_lower1 = apply(Ymi_epil, 2, function(x) quantile(x, 0.25)),
  epilithon_gm2_upper1 = apply(Ymi_epil, 2, function(x) quantile(x, 0.75)),
  filamentous_gm2 = apply(Ymi_fila, 2, median),
  filamentous_gm2_lower = apply(Ymi_fila, 2, function(x) quantile(x, 0.025)),
  filamentous_gm2_upper = apply(Ymi_fila, 2, function(x) quantile(x, 0.965)),
  filamentous_gm2_lower1 = apply(Ymi_fila, 2, function(x) quantile(x, 0.25)),
  filamentous_gm2_upper1 = apply(Ymi_fila, 2, function(x) quantile(x, 0.75)),
  epilithon_chla_mgm2 = apply(Ymi_epil_chla, 2, median),
  epilithon_chla_lower = apply(Ymi_epil_chla, 2, function(x) quantile(x, 0.025)),
  epilithon_chla_upper = apply(Ymi_epil_chla, 2, function(x) quantile(x, 0.975)),
  epilithon_chla_lower1 = apply(Ymi_epil_chla, 2, function(x) quantile(x, 0.25)),
  epilithon_chla_upper1 = apply(Ymi_epil_chla, 2, function(x) quantile(x, 0.75)),
  filamentous_chla_mgm2 = apply(Ymi_fila_chla, 2, median),
  filamentous_chla_lower = apply(Ymi_fila_chla, 2, function(x) quantile(x, 0.025)),
  filamentous_chla_upper = apply(Ymi_fila_chla, 2, function(x) quantile(x, 0.975)),
  filamentous_chla_lower1 = apply(Ymi_fila_chla, 2, function(x) quantile(x, 0.25)),
  filamentous_chla_upper1 = apply(Ymi_fila_chla, 2, function(x) quantile(x, 0.75)),
  frac_fila = apply(Ymi_fila/(Ymi_epil + Ymi_fila), 2, median),
  frac_fila_lower = apply(Ymi_fila/(Ymi_epil + Ymi_fila), 2, function(x) quantile(x, 0.25)),
  frac_fila_upper = apply(Ymi_fila/(Ymi_epil + Ymi_fila), 2, function(x) quantile(x, 0.75))
  )
Ymi_NPP_sum <- data.frame(
  NPP = apply(Ymi_NPP, 2, median),
  NPP_chla = apply(Ymi_NPP_chla, 2, median),
  NPP_lower = apply(Ymi_NPP_chla, 2, function(x) quantile(x, 0.025)),
  NPP_upper = apply(Ymi_NPP_chla, 2, function(x) quantile(x, 0.975))
  )
NPP_mod <- bind_cols(bm_met, select(Ymi_NPP_sum, -NPP)) %>%
  group_by(date, site) %>%
  summarize(NPP_mod = mean(NPP_chla),
            NPP_upper = mean(NPP_upper),
            NPP_lower = mean(NPP_lower)) %>%
  ungroup()
png("figures/ppreds_NPP_full_biomass_mod.png", width = 10, height = 8, units = "in", res = 300)
left_join(select(bm_sum, date, doy, year, site, NPP_est = NPP),
          NPP_mod, by = c("site", "date")) %>%
    mutate(site = factor(site, levels = c("PL", "DL", "GR", "GC", "BG", "BN"))) %>%
    ggplot(aes(date, NPP_mod)) +
    geom_ribbon(aes(ymin = NPP_lower, ymax = NPP_upper),col = 'gray', alpha = 0.4)+
    geom_line() +
    geom_point(aes(y = NPP_est), pch = 3, size = 1, col = 'forestgreen') +
    facet_grid(site~year, scales = "free")+
    ylab("NPP (g C/m2/d)") + xlab("Date") +
    theme_bw()
dev.off()
# Add estimates into dataframe of missing:
bm_sum <- bm_sum2 %>%
  rename_with(~str_replace(.x, '_min_mean', '')) %>%
  mutate(epilithon_gm2_upper = epilithon_gm2 + 1.96 * epilithon_gm2_min_se,
         epilithon_gm2_lower = epilithon_gm2 - 1.96 * epilithon_gm2_min_se,
         filamentous_gm2_upper = filamentous_gm2 + 1.96 * filamentous_gm2_min_se,
         filamentous_gm2_lower = filamentous_gm2 - 1.96 * filamentous_gm2_min_se,
         epilithon_chla_upper = epilithon_chla_mgm2 + 1.96 * epilithon_chla_mgm2_min_se,
         epilithon_chla_lower = epilithon_chla_mgm2 - 1.96 * epilithon_chla_mgm2_min_se,
         filamentous_chla_upper = filamentous_chla_mgm2 + 1.96 * filamentous_chla_mgm2_min_se,
         filamentous_chla_lower = filamentous_chla_mgm2 - 1.96 * filamentous_chla_mgm2_min_se) %>%
  select(-ends_with("min_se"))
which_mi <- which(is.na(bm_sum$epilithon_gm2))

bm_sum$epilithon_gm2[is.na(bm_sum$epilithon_gm2)] <- Ymi_sum$epilithon_gm2
bm_sum$epilithon_gm2_upper[is.na(bm_sum$epilithon_gm2_upper)] <- Ymi_sum$epilithon_gm2_upper1
bm_sum$epilithon_gm2_lower[is.na(bm_sum$epilithon_gm2_lower)] <- Ymi_sum$epilithon_gm2_lower1
bm_sum$filamentous_gm2[is.na(bm_sum$filamentous_gm2)] <- Ymi_sum$filamentous_gm2
bm_sum$filamentous_gm2_upper[is.na(bm_sum$filamentous_gm2_upper)] <- Ymi_sum$filamentous_gm2_upper1
bm_sum$filamentous_gm2_lower[is.na(bm_sum$filamentous_gm2_lower)] <- Ymi_sum$filamentous_gm2_lower1


bm_sum$epilithon_chla_mgm2[is.na(bm_sum$epilithon_chla_mgm2)] <- Ymi_sum$epilithon_chla_mgm2
bm_sum$epilithon_chla_upper[is.na(bm_sum$epilithon_chla_upper)] <- Ymi_sum$epilithon_chla_upper1
bm_sum$epilithon_chla_lower[is.na(bm_sum$epilithon_chla_lower)] <- Ymi_sum$epilithon_chla_lower1
bm_sum$filamentous_chla_mgm2[is.na(bm_sum$filamentous_chla_mgm2)] <- Ymi_sum$filamentous_chla_mgm2
bm_sum$filamentous_chla_upper[is.na(bm_sum$filamentous_chla_upper)] <- Ymi_sum$filamentous_chla_upper1
bm_sum$filamentous_chla_lower[is.na(bm_sum$filamentous_chla_lower)] <- Ymi_sum$filamentous_chla_lower1


bm_sum <- bm_sum %>%
  mutate(across(starts_with("epilithon_gm2"),
                ~ .x - min_vals$min[min_vals$biomass == "epilithon_gm2"]),
         across(starts_with("filamentous_gm2"),
                ~ .x - min_vals$min[min_vals$biomass == "filamentous_gm2"]),
         across(starts_with("epilithon_chla"),
                ~ .x - min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"]),
         across(starts_with("filamentous_chla"),
                ~ .x - min_vals$min[min_vals$biomass == "filamentous_chla_mgm2"]),
         across(starts_with(c("epilithon", "filamentous")),
                ~ if_else(.x < 0, 0, .x)),
         frac_fila = NA,
         frac_fila_upper = NA,
         frac_fila_lower = NA)

bm_sum$frac_fila[which_mi] <- Ymi_sum$frac_fila
bm_sum$frac_fila_upper[which_mi] <- Ymi_sum$frac_fila_upper
bm_sum$frac_fila_lower[which_mi] <- Ymi_sum$frac_fila_lower
bm_sum <- left_join(bm_sum, NPP_mod, by = c("site", "date")) %>%
  left_join(select(light, site, date, PAR_bc_Jm2), by = c("site", "date"))
write_csv(bm_sum, 'data/biomass_data/log_gamma_brms_gam_fits_biomass.csv')


# plot GAMS  ####
qq <- read_csv('data/biomass_data/log_gamma_brms_gam_fits_biomass.csv')
qq <- qq %>%
    select(-NPP, -NPP_chla) %>%
    rename_with(.fn = ~paste0(.x, "_fit"), .cols = ends_with("gm2")) %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site),
           site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BG', 'BN'))) %>%
    pivot_longer(cols = starts_with(c('epil', 'fila')),
             names_to = c('biomass_type', 'units', 'stat'),
             names_pattern = '([a-z]+)_([a-z0-9_]+)_([a-z]+)',
             values_to = 'value') %>%
    mutate(units = if_else(units == "chla_mgm2", "chla", units)) %>%
    pivot_wider(names_from = 'stat', values_from = 'value') %>%
    group_by(site, biomass_type, units) %>%
    mutate(across(c('fit', 'upper', 'lower'), ~zoo::na.approx(., x = date, na.rm = F)))# %>%
    # mutate(fit = case_when(se > fit ~ NA_real_,
    #                        se > 500 ~ NA_real_,
    #                        TRUE ~ fit),
    #        se = case_when(is.na(fit)~NA_real_,
    #                       TRUE ~ se))
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
meas_chl <- filter(meas, units == 'chla_mgm2')
meas_mass <- filter(meas, units == 'gm2')


meas_mass2 <- mutate(meas_mass,
                     biomass_type = case_when(biomass_type == 'epil'~ 'Epilithon',
                                              biomass_type == 'fila' ~ 'Filamentous'),
                     meas = case_when(meas < min_vals$min[min_vals$biomass == "filamentous_gm2"] ~
                                        min_vals$min[min_vals$biomass == "filamentous_gm2"],
                                      TRUE ~ meas))
meas_chl2 <- mutate(meas_chl,
                    biomass_type = case_when(biomass_type == 'epil'~ 'Epilithon',
                                             biomass_type == 'fila' ~ 'Filamentous'),
                    meas = case_when(meas < min_vals$min[min_vals$biomass == "filamentous_chla_mgm2"] ~
                                       min_vals$min[min_vals$biomass == "filamentous_chla_mgm2"],
                                     TRUE ~ meas))

measmeas <- distinct(meas_mass2, site, date, biomass_type) %>% mutate(meas = TRUE)
measmeas <- meas %>% group_by(site, date, biomass_type, units) %>%
  summarize(meas = mean(meas, na.rm = T)) %>%
  mutate(biomass_type = case_when(biomass_type == 'epil'~ 'Epilithon',
                                  biomass_type == 'fila' ~ 'Filamentous'),
         units = case_when(units == "chla_mgm2" ~ "chla",
                           TRUE ~ units))

qq2 <- qq %>% mutate(biomass_type = case_when(biomass_type == 'epilithon'~ 'Epilithon',
                                       biomass_type == 'filamentous' ~ 'Filamentous')) %>%
  left_join(measmeas, by = c("site", "date", "biomass_type", "units"))
qq <- filter(qq2, is.na(meas))

mm <- qq %>% filter(units == 'gm2') %>%
    mutate(fit = case_when(fit < min_vals$min[min_vals$biomass == "filamentous_gm2"] ~
                             min_vals$min[min_vals$biomass == "filamentous_gm2"],
                           TRUE ~ fit),
           lower = case_when(lower < min_vals$min[min_vals$biomass == "filamentous_gm2"] ~
                             min_vals$min[min_vals$biomass == "filamentous_gm2"],
                           TRUE ~ lower),
           upper = case_when(upper < fit ~ fit,
                             TRUE ~ upper),
           year = year(date)) %>%
    ggplot(aes(date, fit, col = biomass_type)) +
    geom_line()+
    geom_ribbon(aes(ymax = upper, ymin = lower,
                    fill = biomass_type), alpha = 0.4, color = NA)+
    geom_point(data = meas_mass2, aes(date, meas, col = biomass_type))+
    facet_grid(site~year, scales = 'free_x') +
    scale_color_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 600))+
    xlab('Date') +
    ylab(expression('Algal Biomass (AFDM g '~ m^-2*')')) +
    theme_classic()+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, 'line'))
cc <- qq %>%
    filter(units == 'chla') %>%
    mutate(fit = case_when(fit < min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"] ~
                             min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"],
                           TRUE ~ fit),
           lower = case_when(lower < min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"] ~
                               min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"],
                             TRUE ~ lower),
           upper = case_when(upper < fit ~ fit,
                               TRUE ~ upper),
           year = year(date)) %>%
    ggplot(aes(date, fit, col = biomass_type)) +
    geom_line()+
    geom_ribbon(aes(ymax = upper, ymin = lower,
                    fill = biomass_type), alpha = 0.4, color = NA)+
    geom_point(data = meas_chl2, aes(date, meas, col = biomass_type))+
    facet_grid(site~year, scales = 'free_x') +
    scale_color_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_fill_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
    scale_y_log10(limits = c(0.3, 1000))+
    xlab('Date') +
    ylab(expression('Algal Chlorophyll (mg chl a '~ m^-2*')')) +
    theme_classic()+
    theme(panel.border = element_rect(fill = NA),
          panel.spacing = unit(0, 'line'))

png('figures/biomass_log_gamma_brms_gams_comb_zeros.png', width = 7.5, height = 5, units = 'in',
    res = 300)

    ggpubr::ggarrange(mm, cc, nrow = 1, common.legend = TRUE,
                      labels = c('(a)', '(b)'))
dev.off()

mm <- qq %>% filter(units == 'gm2') %>%
  mutate(fit = case_when(fit < min_vals$min[min_vals$biomass == "filamentous_gm2"] ~
                           min_vals$min[min_vals$biomass == "filamentous_gm2"],
                         TRUE ~ fit),
         lower = case_when(lower < min_vals$min[min_vals$biomass == "filamentous_gm2"] ~
                             min_vals$min[min_vals$biomass == "filamentous_gm2"],
                           TRUE ~ lower),
         upper = case_when(upper < fit ~ fit,
                           TRUE ~ upper),
         year = year(date)) %>%
  ggplot(aes(date, fit, col = biomass_type)) +
  geom_line()+
  geom_ribbon(aes(ymax = upper, ymin = lower,
                  fill = biomass_type), alpha = 0.4, color = NA)+
  geom_point(data = meas_mass2, aes(date, meas, col = biomass_type))+
  facet_grid(site~year, scales = 'free_x') +
  scale_color_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
  scale_fill_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
  xlab('Date') +ylim(0, 300) +
  ylab(expression('Algal Biomass (AFDM g '~ m^-2*')')) +
  theme_classic()+
  theme(panel.border = element_rect(fill = NA),
        panel.spacing = unit(0, 'line'))
cc <- qq %>%
  filter(units == 'chla') %>%
  mutate(fit = case_when(fit < min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"] ~
                           min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"],
                         TRUE ~ fit),
         lower = case_when(lower < min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"] ~
                             min_vals$min[min_vals$biomass == "epilithon_chla_mgm2"],
                           TRUE ~ lower),
         upper = case_when(upper < fit ~ fit,
                           TRUE ~ upper),
         year = year(date)) %>%
  ggplot(aes(date, fit, col = biomass_type)) +
  geom_line()+
  geom_ribbon(aes(ymax = upper, ymin = lower,
                  fill = biomass_type), alpha = 0.4, color = NA)+
  geom_point(data = meas_chl2, aes(date, meas, col = biomass_type))+
  facet_grid(site~year, scales = 'free_x') +
  scale_color_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
  scale_fill_discrete("Biomass Type", type = c('#1B9EC9', '#97BB43'))+
  xlab('Date') +
  ylab(expression('Algal Chlorophyll (mg chl a '~ m^-2*')')) +
  theme_classic()+
  theme(panel.border = element_rect(fill = NA),
        panel.spacing = unit(0, 'line'))


png('figures/biomass_gamma_brms_gams_comb_zeros.png', width = 7.5, height = 5, units = 'in',
    res = 300)

    ggpubr::ggarrange(mm, cc, nrow = 1, common.legend = TRUE,
                      labels = c('(a)', '(b)'))
dev.off()
qq3 <- qq2 %>%
  group_by(site, biomass_type, units) %>%
  mutate(meas = zoo::na.approx(meas, x = date, na.rm = F))
png("figures/ppdists_algal_biomass.png", width = 6, height = 5, units = "in", res = 300)
qq3 %>%
  pivot_longer(cols = c(fit, meas), names_to = 'group', values_to = 'val') %>%
  mutate(group = case_when(group == "fit" ~ "Modeled",
                           group == "meas" ~ "Measured"),
         units = case_when(units == "gm2" ~ "AFDM (g/m2)",
                           units == "chla" ~ "Chl a (mg/m2)"),
         biomass = paste(biomass_type, units, sep = "_")) %>%
  ggplot(aes(val, group = group, fill = group, col = group))+
  geom_density(adjust = 1, alpha = 0.2) +
  facet_grid(units ~ biomass_type, scales = "free")+
  xlab("Algal Biomass, Modeled and Measured") +
  ggtitle("Posterior predictive distributions of algal biomass")+
  theme_bw()
dev.off()
