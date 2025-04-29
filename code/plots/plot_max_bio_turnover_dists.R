library(tidyverse)
library(rstan)
library(lme4)
library(brms)
# Growth rate calculations ####

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
met <- left_join(met, select(q90, site, year, ARf = slope, ARf.se = slope.se), by = c('site', 'year')) %>%
    mutate(year = factor(year),
           NPP = GPP * (1-ARf),
           NPP.se = sqrt(GPP.se ^2 + ARf.se^2),
           # NPP_globalARf = GPP * (1 + beta[2,1]),
           AR = -GPP *(ARf),
           AR.se = sqrt(GPP.se^2 + ARf.se^2),
           HR = ER - AR,
           HR.se = sqrt(ER.se^2 + AR.se^2)) %>%
    select(-msgs.fit, -warnings, -errors, -K600, -DO_fit)

biogams <- read_csv('data/biomass_data/log_gamma_gam_fits_biomass.csv') %>%
    mutate(epilg = log(epil_gm2_fit),
           epilgse = abs(epil_gm2_se/epil_gm2_fit),
           filag = log(fila_gm2_fit + 1),
           filagse = abs(fila_gm2_se/(fila_gm2_fit+1)),
           epilc = log(epil_chla_mgm2_fit),
           epilcse = abs(epil_chla_mgm2_se/epil_chla_mgm2_fit),
           filac = log(fila_chla_mgm2_fit + 1),
           filacse = abs(fila_chla_mgm2_se/(fila_chla_mgm2_fit+1)))

# need to get light into the correct units! ####
light <- read_csv('data/site_data/daily_modeled_light_all_sites.csv') %>%
    mutate(site = case_when(site == 'BM' ~ 'BG',
                            TRUE ~ site))
ggplot(light, aes(date, PAR_bc_Jm2, col = site)) +
    geom_line()
ggplot(light, aes(date, LAI, col = site)) +
    geom_line()

# create data frame for models
bm_met <- select(biogams, site, date, epil_gm2_fit, fila_gm2_fit,
                 epil_chla_mgm2_fit, fila_chla_mgm2_fit,
                 epil_gm2_se, fila_gm2_se,
                 epil_chla_mgm2_se, fila_chla_mgm2_se,
                 epilg, epilgse, filag, filagse,
                 epilc, epilcse, filac, filacse) %>%
    rename_with(~gsub('_fit', '', .x)) %>%
    left_join(select(ungroup(met), site, date, year, GPP, ER, ARf, NPP,
                     GPP.se, ER.se, ARf.se, NPP.se),
              by = c('site', 'date')) %>%
    left_join(select(light, site, date, PAR_bc_Jm2)) %>%
    mutate(light = PAR_bc_Jm2/max(PAR_bc_Jm2),
           epilcl = epilc*light,
           epilcsel = epilcse*light,
           filacl = filac*light,
           filacsel = filacse*light,
           epil_chla_mgm2_l = epil_chla_mgm2*light,
           epil_chla_mgm2_sel = epil_chla_mgm2_se*light) %>%
    filter(!is.na(GPP))

# write_csv(bm_met, "data_for_model.csv")
# bm_met <- read_csv("data_for_model.csv")

# Simulate data ####
N = nrow(bm_met)
meanme <- c(40, 30)
sdme <- c(12, 8)
B_f = rgamma(N, shape = mean(bm_met$fila_chla_mgm2)^2/sd(bm_met$fila_chla_mgm2),
             rate = mean(bm_met$fila_chla_mgm2)/sd(bm_met$fila_chla_mgm2))
B_e = rgamma(N, shape = mean(bm_met$epil_chla_mgm2)^2/sd(bm_met$epil_chla_mgm2),
             rate = mean(bm_met$epil_chla_mgm2)/sd(bm_met$epil_chla_mgm2))

B_f_se = rgamma(N, 1/3, 1/3)
B_e_se = rgamma(N, 1/2, 1/2)

B_f_mod <- rgamma(N, shape = B_f^2/B_f_se, rate = B_f/B_f_se)
B_e_mod <- rgamma(N, shape = B_e^2/B_e_se, rate = B_e/B_e_se)

site = factor(bm_met$site)
light = bm_met$PAR_bc_Jm2

# define coefficients:
mu_f = 0.1
sigma_f = 0.02
mu_fj = rnorm(6, mu_f, sigma_f)
mu_e = 0.5
K_Bf = 0.2
K_I = 5
sigma = 1

# Simulate NPP:
NPP = mu_f * B_f + mu_e * B_e + rnorm(N, 0, sigma)
# NPP = (mu_f * B_f * (1/(1 + K_Bf * B_f)) + mu_e * B_e) * (light/(K_I + light)) +
NPP = (mu_fj[site] * B_f + mu_e * B_e) * (light/(K_I + light)) +
    rnorm(N, 0, sigma)


# test brms model:
dd <- data.frame(NPP = NPP,
                 B_e_mod = B_e_mod,
                 B_f_mod = B_f_mod,
                 B_e_se = B_e_se,
                 B_f_se = B_f_se)
ll <- lm(NPP ~ 0 + B_e_mod + B_f_mod, data = dd)
summary(ll)
bform <- bf(NPP ~ 0 + me(B_e_mod, B_e_se) + me(B_f_mod, B_f_se))
stancode(bform, data = dd)
bb <- brm(NPP ~ 0 + me(B_e_mod, B_e_se) + me(B_f_mod, B_f_se),
          data = dd,
          ncores = 4, nchains = 4)
summary(bb)


# Fit STAN model ####
# Prepare data for Stan
sim_list <- list(
    N = nrow(bm_met),
    NPP = NPP,
    Mme = 2,
    NCme = 1,
    Bf_mod = B_f_mod,
    Be_mod = B_e_mod,
    Bf_se = B_f_se,
    Be_se = B_e_se
    # light = light
)

# Compile the Stan model
stan_model <- stan_model(file = "code/model/stan_code/partition_NPP_brms_build.stan")
# stan_model <- stan_model(file = "code/model/stan_code/partition_NPP_error.stan")
# stan_model <- stan_model(file = "code/model/stan_code/partition_NPP_error2.stan")

# init_list = list(list(mu_f = 0.2,
#                       sigma_f = 0.1),
#                  list(mu_f = 0.3,
#                       sigma_f = 0.05),
#                  list(mu_f = 0.25,
#                       sigma_f = 0.12),
#                  list(mu_f = 0.22,
#                       sigma_f = 0.08))

# Run the Stan model
fit_sim_lb <- sampling(stan_model,
                    data = sim_list,
                    iter = 2000,
                    # init = init_list,
                    chains = 4,
                    cores = 4)

# Evaluate the output
summary(bb)
summary(fit_sim)
summary(fit_sim_lb, pars = c("mu_f", "mu_e", "meanme", "sdme", "sigma"))
traceplot(fit_sim_lb, pars = c("mu_f", "mu_e", "meanme", "sdme", "sigma"))
pairs(fit_sim, pars = c("mu_f", "mu_e", "meanme", "sdme", "sigma"))
summary(fit_sim, pars = c("mu_f", "sigma_f", "mu_e", "K_I", "sigma"))
traceplot(fit_sim, pars = c("mu_f", "mu_e", "K_I", "sigma"))
plot(fit_sim, pars = c("mu_f", "mu_e", "K_Bf", "K_I", "sigma"))
pairs(fit_sim, pars = c("mu_f", "mu_e", "K_Bf", "K_I", "sigma"))

sim_ppreds <- summary(fit_sim)$summary %>%
    data.frame() %>%
    mutate(param = rownames(summary(fit_sim)$summary))

plot(density(sim_ppreds$mean), xlim = c(0, 30), lty = 2)
lines(density(NPP))

# Fit model on real data####
lm(NPP/light ~ 0 + fila_chla_mgm2 + epil_chla_mgm2, data = bm_met)

dat_list <- list(
    N = nrow(bm_met),
    Mme = 2,
    NCme = 1,
    NPP = bm_met$NPP,
    B_f_mean = bm_met$fila_chla_mgm2,
    B_e_mean = bm_met$epil_chla_mgm2,
    B_f_se = bm_met$fila_chla_mgm2_se,
    B_e_se = bm_met$epil_chla_mgm2_se
    )
    # light = bm_met$PAR_bc_Jm2)

fit <- sampling(stan_model,
                data = dat_list,
                iter = 3000,
                control = list(max_treedepth = 15),
                chains = 4,
                cores = 4)

summary(fit)
summary(modmix_brma)
fit@stanmodel
traceplot(fit)
meanme <- c(7.2e-8, 25)
sdme <- c(7.46e-7, 8.19)

bfmix <- bf(NPP ~ 0 +  me(fila_chla_mgm2, fila_chla_mgm2_se) + me(epil_chla_mgm2, epil_chla_mgm2_se) )

get_prior(bfmix, data = bm_met, family = gaussian)

priorsa <- c(
    prior(normal(0, 0.1), class = "b", coef = "mefila_chla_mgm2fila_chla_mgm2_se"),
    prior(normal(0, 0.3), class = "b", coef = "meepil_chla_mgm2epil_chla_mgm2_se"),
    prior(normal(0, 1), class = "meanme"))

brms::stancode(bfmix, data = bm_met, family = gaussian, prior = priorsa)

modmix_brma <- brms::brm(bfmix,
                         data = bm_met,
                         prior = priorsa,
                         iter = 4000,
                         control = options(max_treedepth = 16,
                                           adapt_delta = 0.9),
                         chains = 4, cores = 4)


summary(modmix_brma)
mean(bm_met$epil_chla_mgm2)
mean(bm_met$fila_chla_mgm2)
Lme <- matrix(data = c(1,0, 0.7, 0.96), ncol = 2)
zme <- matrix(data = c(0.17,0.2,0.188,0.241,
                       1,1,1,1), ncol = 2)

matrix(c(rep(meanme[1], 4), rep(meanme[2], 4)), ncol = 2) + t((diag(sdme) %*% Lme) %*% t(zme))
bm_met %>% select(fila_chla_mgm2, epil_chla_mgm2)

summary(fit, pars = c("mu_f", "mu_e", "K_I", "sigma"))
traceplot(fit, pars = c("mu_f", "mu_e", "K_I", "sigma"))
traceplot(fit, pars = c("mu_f", "mu_e", "K_I", "sigma"))
summary(fit, pars = c("mu_f", "sigma_f", "mu_e", "sigma_e", "K_I", "sigma"))
traceplot(fit, pars = c("mu_f", "sigma_f", "mu_e", "sigma_e", "K_I", "sigma"))
pairs(fit, pars = c("mu_f", "mu_e", "K_I", "sigma"))

fit_ppreds <- summary(fit)$summary %>%
    data.frame() %>%
    mutate(param = rownames(summary(fit)$summary))

plot(density(fit_ppreds$mean), xlim = c(-5, 20), lty = 2)
lines(density(bm_met$NPP))

NPP_meas <- bm_met$NPP
shinystan::launch_shinystan(fit)

+ epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod1 <- lm(NPP ~ 0 + epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod2 <- lm(NPP ~ 0 + epil_chla_mgm2_l + fila_chla_mgm2_l, bm_met)
mod3 <- lm(NPP/light ~ 0 + epil_chla_mgm2 + fila_chla_mgm2, bm_met)
mod4 <- lm(NPP/light ~ 0 + epil_gm2 + fila_gm2, bm_met)

modmix1 <- lme4::lmer(NPP/light ~ 0 + epil_chla_mgm2 + fila_chla_mgm2 +
                         (0+epil_chla_mgm2+fila_chla_mgm2|site),
                     data = bm_met,
                     REML = FALSE,
                     control = lmerControl(optimizer ="Nelder_Mead"))
modmix2 <- lme4::lmer(NPP/light ~ 0 + epil_chla_mgm2 + filac +
                         (0+epil_chla_mgm2+filac|site),
                     data = bm_met,
                     REML = FALSE,
                     control = lmerControl(optimizer ="Nelder_Mead"))
modmix3 <- lme4::lmer(NPP/light ~ 0 + epilc + filac +
                         (0+epilc+filac|site),
                     data = bm_met,
                     REML = FALSE,
                     control = lmerControl(optimizer ="Nelder_Mead"))
modmix4 <- lme4::lmer(NPP ~ 0 + epilcl + filacl +
                         (0+epilcl+filacl|site),
                     data = bm_met,
                     REML = FALSE,
                     control = lmerControl(optimizer ="Nelder_Mead"))
modmix5 <- lme4::lmer(NPP ~ 0 + epil_chla_mgm2_l + filacl +
                         (0+epil_chla_mgm2_l+filacl|site),
                     data = bm_met,
                     REML = FALSE,
                     control = lmerControl(optimizer ="Nelder_Mead"))

summary(modmix1)
ranef(modmix)

ggplot(bm_met, aes(fila_chla_mgm2, NPP/light))+
    geom_point() +
    geom_errorbarh(aes(xmin = fila_chla_mgm2 - fila_chla_mgm2_se,
                       xmax = fila_chla_mgm2 + fila_chla_mgm2_se))

# build model using brms:
bm_met <- bm_met %>%
    mutate(NPPL = NPP/light,
           filacl = filac*light,
           epilcl = epilc*light,
           filacsel = filacse*light,
           epilcsel = epilcse*light)

hist(bm_met$NPP)
hist(log(bm_met$NPPL))
hist(bm_met$epil_chla_mgm2)
hist(bm_met$fila_chla_mgm2)
hist(bm_met$epilc)
hist(bm_met$epilcse)
hist(bm_met$filac)
hist(bm_met$filacse)

bm_met$siteyear = as.factor(paste0(bm_met$site, bm_met$year))


# Plot prior distributions:
median_epil <- mean(1/(bm_met$light*4))
max_epil <- mean(1/(bm_met$light*30))
min_epil <- mean(1/(bm_met$light*1))
median_epil <- mean(bm_met$epil_chla_mgm2/(bm_met$light*bm_met$epilc*0.1))
median_epil <- mean(bm_met$epil_chla_mgm2/(bm_met$light *bm_met$epilc*4))
max_fila <- mean(bm_met$fila_chla_mgm2/(bm_met$light*bm_met$filac*100))
median_fila <- mean(bm_met$fila_chla_mgm2/(bm_met$light *bm_met$filac*10))
min_fila <- mean(bm_met$fila_chla_mgm2/(bm_met$light*bm_met$filac*2))

png("figures/priors_for_NPP_partition_model.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,1), mar = c(4,4,1,1), oma = c(0,0,1,0))
x = seq(0, 12, by = 0.01)
plot(x, dlnorm(x, 0.4, 1), type = 'l', xlab = bquote("Value of "*mu ["f"]), ylab = "density")
mtext("Priors for Algal Growth Rates", adj = 0, line = 0.2)
mtext("Filamentous Algae", line = -1, adj = 0.9)
polygon(c(x[x<max_fila], max_fila, max_fila, rev(x[x<max_fila])),
        c(dlnorm(c(x[x<max_fila], max_fila), 0.4, 1), rep(0, length(x[x<max_fila]) + 1)),
        col = "purple")
polygon(c(min_fila, x[x>min_fila], rev(x[x>min_fila]), min_fila),
        c(dlnorm(c(min_fila, x[x>min_fila]), 0.4, 1), rep(0, length(x[x>min_fila]) + 1)),
        col = "purple")
text(1.45, 0.05, "p(turnover time \n> 100 days)", col = "purple", cex = 0.8)
text(3, 0.3, "Median expected \nturnover time \n= 10 days", cex = 0.8)
text(9.5, 0.04, "p(turnover time < 2 days)", col = "purple", cex = 0.8)
abline(v = median_fila, lty = 2)

x = seq(0, 3.5, by = 0.01)
plot(x, dlnorm(x, -0.9, 1), type = 'l', xlab = bquote("Value of "*mu ["e"]), ylab = "density")
mtext( "Epilithic Algae", line = -1, adj = 0.9)
polygon(c(x[x<max_epil], max_epil, max_epil, rev(x[x<max_epil])),
        c(dlnorm(c(x[x<max_epil], max_epil), -0.9, 1), rep(0, length(x[x<max_epil]) + 1)),
        col = "purple")
polygon(c(min_epil, x[x>min_epil], rev(x[x>min_epil]), min_epil),
        c(dlnorm(c(min_epil, x[x>min_epil]), -0.9, 1), rep(0, length(x[x>min_epil]) + 1)),
        col = "purple")
text(0.45, 0.18, "p(turnover time \n> 30 days)", col = "purple", cex = 0.8)
text(.85, 1.2, "Median expected \nturnover time \n= 4 days", cex = 0.8)
text(2.5, 0.18, "p(turnover time < 1 day)", col = "purple", cex = 0.8)
abline(v = median_epil, lty = 2)
dev.off()

example(stan_model, package = "rstan", run.dontrun = TRUE)
# Fit BRMS models: ####
# Model 1: NPPL = (epil) + (fila)) * L
bfmix <- bf(NPP ~ 0 +  me(fila_chla_mgm2, fila_chla_mgm2_se) + me(epil_chla_mgm2, epil_chla_mgm2_se) )
bfmix <- bf(NPPL ~ 0 +  me(epil_chla_mgm2, epil_chla_mgm2_se) + me(fila_chla_mgm2, fila_chla_mgm2_se) +
                (0 + me(epil_chla_mgm2, epil_chla_mgm2_se) + me(fila_chla_mgm2, fila_chla_mgm2_se)|siteyear))

get_prior(bfmix, data = bm_met, family = gaussian)
brms::stancode(bfmix, data = bm_met, family = gaussian, prior = priorsa)

priorsa <- c(
    prior(normal(0, 0.1), class = "b", coef = "mefila_chla_mgm2fila_chla_mgm2_se"),
    prior(normal(0, 0.3), class = "b", coef = "meepil_chla_mgm2epil_chla_mgm2_se"),
    prior(normal(0, 1), class = "meanme"))

priorsb <- c(
    prior(lognormal(-2.5, 1), class = "b", coef = "mefilaclfilacsel"),
    prior(lognormal(-1.2, 1), class = "b", coef = "meepilclepilcsel"),
    prior(normal(0, 0.2), class = "sd")
)

modmix_brma <- brms::brm(bfmix,
                        data = bm_met,
                        prior = priorsa,
                        iter = 4000,
                        control = options(max_treedepth = 16,
                                          adapt_delta = 0.9),
                        chains = 4, cores = 4)
# saveRDS(modmix_brma, "data/model_fits/brms_NPPL_partition_linElinF_gaussPrior.rds")
modL_linElinF_gauss <- readRDS("data/model_fits/model_fits/brms_NPPL_partition_linElinF_gaussPrior.rds")
summary(modL_linElinF_gauss)
# modmix_brmb <- brms::brm(bfmix,
#                         data = bm_met,
#                         family = gaussian(link = "identity"),
#                         prior = priorsb,
#                         iter = 6000,
#                         control = options(max_treedepth = 16,
#                                           adapt_delta = 0.9),
#                         chains = 4, cores = 4)
# saveRDS(modmix_brmb, "data/model_fits/brms_NPPL_partition_linElinF_logPrior.rds")
modL_linElinF_log <- readRDS("data/model_fits/model_fits/brms_NPPL_partition_linElinF_logPrior.rds")
summary(modL_linElinF_log)
# plot(modL_linElinF_log)
pp_check(modL_linElogF_log, ndraws = 100)
preds_linlin_gauss = predict(modL_linElinF_gauss, newdata = bm_met)
preds_linlin_log = predict(modL_linElinF_log, newdata = bm_met)

mods <- data.frame(
    response = c("NPP/L", "NPP/L"),
    preds = c("linElinF", "linElinF"),
    prior = c("gauss", "lnorm"),
    rmse = c(sqrt(mean((preds_linlin_gauss[,'Estimate'] - bm_met$NPPL)^2)),
             sqrt(mean((preds_linlin_log[,'Estimate'] - bm_met$NPPL)^2))),
    coef_epil = c(fixef(modL_linElinF_gauss)[1,1], fixef(modL_linElinF_log)[1,1]),
    coef_fila = c(fixef(modL_linElinF_gauss)[2,1], fixef(modL_linElinF_log)[2,1])
)

bm_met %>%
    bind_cols(preds_linlin_gauss) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP/light, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP/L ~ linE + linF, gauss prior")+
    theme_bw()
bm_met %>%
    bind_cols(preds_linlin_log) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP/light, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP/L ~ linE + linF, lnorm prior")+
    theme_bw()

# Model 2:  NPP/L = epil + log(fila) ####
bfmix <- bf(NPPL ~ 0 +  me(epil_chla_mgm2, epil_chla_mgm2_se) + me(filac, filacse) +
                (0 + me(epil_chla_mgm2, epil_chla_mgm2_se) + me(filac, filacse)|siteyear))

priors <- c(
    prior(normal(0, 0.5), class = "b", coef = "mefilaclfilacsel"),
    prior(normal(0, 0.5), class = "b", coef = "meepil_chla_mgm2epil_chla_mgm2_se"),
    prior(normal(0, 0.2), class = "sd")
)

# modmix_brmlinlog <- brms::brm(bfmix,
#                               data = bm_met,
#                               prior = priors,
#                               iter = 10000,
#                               control = options(max_treedepth = 15,
#                                                 adapt_delta = 0.9),
#                               chains = 4, cores = 4)
#
# saveRDS(modmix_brmlinlog, "data/model_fits/brms_NPPL_partition_linElogF.rds")
modL_linElogF_gauss <- readRDS("data/model_fits/model_fits/brms_NPPL_partition_linElogF.rds")
modL_linElogF_log <- readRDS("data/model_fits/brms_NPPL_partition_linElogF.rds")
# modmix_brmlinlog <- readRDS("data/model_fits/brms_NPPL_partition_linElogF.rds")
summary(modL_linElogF_gauss)
get_prior(modL_linElogF_gauss)
# plot(modL_linElogF_gauss)
# pp_check(modL_linElogF_gauss, ndraws = 100)

preds_linlog = predict(modL_linElogF_gauss, newdata = bm_met)
sqrt(mean((preds_linlog[,'Estimate'] - bm_met$NPPL)^2))
preds_linlog_log = predict(modL_linElogF_log, newdata = bm_met)
sqrt(mean((preds_linlog_log[,'Estimate'] - bm_met$NPPL)^2))

mods <- data.frame(
    response = c("NPP/L", "NPP/L"),
    preds = c("linElogF", "linElogF"),
    prior = c("gauss", "lnorm"),
    rmse = c(sqrt(mean((preds_linlog[,'Estimate'] - bm_met$NPPL)^2)),
             sqrt(mean((preds_linlog_log[,'Estimate'] - bm_met$NPPL)^2))),
    coef_epil = c(fixef(modL_linElogF_gauss)[1,1], fixef(modL_linElogF_log)[1,1]),
    coef_fila = c(fixef(modL_linElogF_gauss)[2,1], fixef(modL_linElogF_log)[2,1])
) %>%
    bind_rows(mods)

bm_met %>%
    bind_cols(preds_linlog) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP/light, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP/L ~ linE + logF, gauss prior")+
    theme_bw()

bm_met %>%
    bind_cols(preds_linlog_log) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP/light, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP/L ~ linE + logF, lnorm prior")+
    theme_bw()


# Model 2: NPPL = (log(epil) + log(fila))
bfmix <- bf(NPPL ~ 0 +  me(epilc, epilcse) + me(filac, filacse) +
                (0 + me(epilc, epilcse) + me(filac, filacse)|siteyear))
get_prior(bfmix, data = bm_met, family = gaussian)

# priorsa <- c(
#     prior(normal(0, 0.2), class = "b", coef = "mefilaclfilacsel"),
#     prior(normal(0, 1), class = "b", coef = "meepilclepilcsel"),
#     prior(normal(0, 0.2), class = "sd")
# )
# priorsb <- c(
#     prior(lognormal(-2, 1), class = "b", coef = "mefilaclfilacsel"),
#     prior(lognormal(0.7, 1), class = "b", coef = "meepilclepilcsel"),
#     prior(normal(0, 0.2), class = "sd")
# )

# modmix_brma <- brms::brm(bfmix,
#                         data = bm_met,
#                         family = gaussian(link = "identity"),
#                         prior = priorsa,
#                         iter = 6000,
#                         control = options(max_treedepth = 16,
#                                           adapt_delta = 0.9),
#                         chains = 4, cores = 4)
# saveRDS(modmix_brma, "data/model_fits/brms_NPPL_partition_logElogF_gaussPrior.rds")
modL_logElogF_gauss <- readRDS("data/model_fits/model_fits/brms_NPPL_partition_logElogF_gaussPrior.rds")
# modmix_brmb <- brms::brm(bfmix,
#                         data = bm_met,
#                         family = gaussian(link = "identity"),
#                         prior = priorsb,
#                         iter = 6000,
#                         control = options(max_treedepth = 16,
#                                           adapt_delta = 0.9),
#                         chains = 4, cores = 4)
# saveRDS(modmix_brmb, "data/model_fits/brms_NPPL_partition_logElogF_logPrior.rds")
modL_logElogF_log <- readRDS("data/model_fits/model_fits/brms_NPPL_partition_logElogF_logPrior.rds")
get_prior(modL_logElogF_gauss)
summary(modL_logElogF_gauss)
summary(modL_logElogF_log)
# pp_check(modL_logElogF_gauss, ndraws = 100)
# pp_check(modL_logElogF_log, ndraws = 100)
# plot(modL_logElogF_log)

preds_loglog_gauss = predict(modL_logElogF_gauss, newdata = bm_met)
preds_loglog_log = predict(modL_logElogF_log, newdata = bm_met)
sqrt(mean((preds_loglog_gauss[,'Estimate'] - bm_met$NPP/bm_met$light)^2))
sqrt(mean((preds_loglog_log[,'Estimate'] - bm_met$NPP/bm_met$light)^2))

mods <- data.frame(
    response = c("NPP/L", "NPP/L"),
    preds = c("logElogF", "logElogF"),
    prior = c("gauss", "lnorm"),
    rmse = c(sqrt(mean((preds_loglog_gauss[,'Estimate'] - bm_met$NPPL)^2)),
             sqrt(mean((preds_loglog_log[,'Estimate'] - bm_met$NPPL)^2))),
    coef_epil = c(fixef(modL_logElogF_gauss)[1,1], fixef(modL_logElogF_log)[1,1]),
    coef_fila = c(fixef(modL_logElogF_gauss)[2,1], fixef(modL_logElogF_log)[2,1])
) %>%
    bind_rows(mods)

bm_met %>%
    bind_cols(preds_loglog_gauss) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP/light, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP/L ~ logE + logF, gauss prior")+
    theme_bw()

bm_met %>%
    bind_cols(preds_loglog_log) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP/light, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP/L ~ logE + logF, lnorm prior")+
    theme_bw()

# Model 1.5: NPP = (log(epil) + log(fila)) * L
bfmix_np <- bf(NPP ~ 0 +  me(epil_chla_mgm2_l, epil_chla_mgm2_sel) + me(filacl, filacsel))
get_prior(bfmix_np, data = bm_met, family = gaussian)

priorsa <- c(
    prior(normal(0, 1), class = "b", coef = "mefilaclfilacsel"),
    prior(normal(0, 1), class = "b", coef = "meepil_chla_mgm2_lepil_chla_mgm2_sel")
)
priorsb <- c(
    prior(exponential(1), class = "b", coef = "mefilaclfilacsel"),
    prior(exponential(1), class = "b", coef = "meepil_chla_mgm2_lepil_chla_mgm2_sel")
)

modmix_brm <- brms::brm(bfmix_np,
                        data = bm_met,
                        family = gaussian(link = "identity"),
                        prior = priorsa,
                        iter = 10000,
                        control = list(max_treedepth = 16),
                        chains = 4, cores = 4)
saveRDS(modmix_brm, "data/model_fits/brms_NPP_pooled_linElogF_gaussprior.rds")
modmix_brmexp <- brms::brm(bfmix_np,
                        data = bm_met,
                        family = gaussian(link = "identity"),
                        prior = priorsb,
                        iter = 10000,
                        control = list(max_treedepth = 15),
                        chains = 4, cores = 4)
saveRDS(modmix_brm, "data/model_fits/brms_NPP_pooled_linElogF_expprior.rds")

#### set two
bfmix_np <- bf(NPP ~ 0 +  me(epil_chla_mgm2_l, epil_chla_mgm2_sel) + me(fila_chla_mgm2_l, fila_chla_mgm2_sel))
get_prior(bfmix_np, data = bm_met, family = gaussian)

priorsa <- c(
    prior(normal(0, 0.2), class = "b", coef = "mefila_chla_mgm2_lfila_chla_mgm2_sel"),
    prior(normal(0, .5), class = "b", coef = "meepil_chla_mgm2_lepil_chla_mgm2_sel")
)
priorsb <- c(
    prior(exponential(1), class = "b", coef = "mefila_chla_mgm2_lfila_chla_mgm2_sel"),
    prior(exponential(1), class = "b", coef = "meepil_chla_mgm2_lepil_chla_mgm2_sel")
)

modmix_brm <- brms::brm(bfmix_np,
                        data = bm_met,
                        family = gaussian(link = "identity"),
                        prior = priorsa,
                        iter = 10000,
                        control = list(max_treedepth = 15),
                        chains = 4, cores = 4)
saveRDS(modmix_brm, "data/model_fits/brms_NPP_pooled_linElinF_gaussprior.rds")
modmix_brmexp <- brms::brm(bfmix_np,
                        data = bm_met,
                        family = gaussian(link = "identity"),
                        prior = priorsb,
                        iter = 10000,
                        control = list(max_treedepth = 15),
                        chains = 4, cores = 4)
saveRDS(modmix_brm, "data/model_fits/brms_NPP_pooled_linElinF_expprior.rds")





summary(modmix_brm)
pairs(modmix_brm)
plot(modmix_brm)
summary(modmix_brmexp)
pairs(modmix_brmexp)
plot(modmix_brmexp)
saveRDS(modmix_brm_a, "data/model_fits/brms_NPP_partition_logElogF_gaussPrior.rds")
mod_logElogF_gauss <- readRDS("data/model_fits/model_fits/brms_NPP_partition_logElogF_gaussPrior.rds")
# modmix_brm_b <- brms::brm(bfmix_np,
#                         data = bm_met,
#                         family = gaussian(link = "identity"),
#                         prior = priorsb,
#                         iter = 6000,
#                         control = options(max_treedepth = 15,
#                                           adapt_delta = 0.9),
#                         chains = 4, cores = 4)
# saveRDS(modmix_brm_b, "data/model_fits/brms_NPP_partition_logElogF_logPrior.rds")
mod_logElogF_log <- readRDS("data/model_fits/model_fits/brms_NPP_partition_logElogF_logPrior.rds")
# get_prior(modmix_brm)
summary(mod_logElogF_gauss)
summary(mod_logElogF_log)
# plot(modmix_brm)
# pp_check(modmix_brm, ndraws = 100)

preds_ll_gauss = predict(mod_logElogF_gauss, newdata = bm_met)
preds_ll_log = predict(mod_logElogF_log, newdata = bm_met)
sqrt(mean((preds_ll_gauss[,'Estimate'] - bm_met$NPP)^2))
sqrt(mean((preds_ll_log[,'Estimate'] - bm_met$NPP)^2))

mods <- data.frame(
    response = c("NPP", "NPP"),
    preds = c("logElogF", "logElogF"),
    prior = c("gauss", "lnorm"),
    rmse = c(sqrt(mean((preds_ll_gauss[,'Estimate'] - bm_met$NPP)^2)),
             sqrt(mean((preds_ll_log[,'Estimate'] - bm_met$NPP)^2))),
    coef_epil = c(fixef(mod_logElogF_gauss)[1,1], fixef(mod_logElogF_log)[1,1]),
    coef_fila = c(fixef(mod_logElogF_gauss)[2,1], fixef(mod_logElogF_log)[2,1])
) %>%
    bind_rows(mods)

bm_met %>%
    bind_cols(preds_ll_gauss) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP ~ logE + logF, gauss prior")+
    theme_bw()

bm_met %>%
    bind_cols(preds_ll_log) %>%
    mutate(predsl = Estimate,
           predsl_q2.5 = Q2.5,
           predsl_q97.5 = Q97.5) %>%
    ggplot(aes(NPP, predsl, col = site)) +
    geom_errorbar(aes(ymin = predsl_q2.5, ymax = predsl_q97.5), alpha = 0.3)+
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    ggtitle("NPP ~ logE + logF, lnorm prior")+
    theme_bw()



# do calculations with bayesian model output:
coefs <- readRDS("data/full_biomass_NPP_partition_mod_chla_posterior_coefs.rds")
bm_met <- read_csv("data/biomass_data/log_gamma_brms_gam_fits_biomass.csv")

mu_epil <- coefs$coefs_df$mu_epil * 12/32
mu_fila <- coefs$coefs_df$mu_fila * 12/32

bm_met <- bm_met %>%
    mutate(light = PAR_bc_Jm2/max(PAR_bc_Jm2)) %>%
    select(site, date, doy, year, light,
           starts_with(c("NPP", "epilithon", "filamentous", "frac"))) %>%
    mutate(lscale = light/(light + 0.5),
           mu_epil_med = mu_epil[1]*lscale,
           mu_epil_low = mu_epil[2]*lscale,
           mu_epil_high = mu_epil[3]*lscale,
           mu_fila_med = mu_fila[1]*lscale,
           mu_fila_low = mu_fila[2]*lscale,
           mu_fila_high = mu_fila[3]*lscale,
           fila_prodgCd_med = filamentous_chla_mgm2 * mu_fila[1] * lscale,
           fila_prodgCd_low = filamentous_chla_lower * mu_fila[2] * lscale,
           fila_prodgCd_high = filamentous_chla_upper * mu_fila[3] * lscale,
           epil_prodgCd_med = epilithon_chla_mgm2 * mu_epil[1] * lscale,
           epil_prodgCd_low = epilithon_chla_upper * mu_epil[3] * lscale,
           epil_prodgCd_high = epilithon_chla_lower * mu_epil[2] * lscale,
           fila_turn_med = filamentous_gm2/fila_prodgCd_med/2,
           fila_turn_high = filamentous_gm2_lower/fila_prodgCd_high/2,
           fila_turn_low = filamentous_gm2_upper/fila_prodgCd_low/2,
           epil_turn_med = epilithon_gm2/epil_prodgCd_med/2,
           epil_turn_high = epilithon_gm2_lower/epil_prodgCd_high/2,
           epil_turn_low = epilithon_gm2_upper/epil_prodgCd_low/2
           )


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
bm_fila <- filter(bm_met, filamentous_gm2 > min_fila_gm2) %>%
    pivot_longer(cols = ends_with(c('med', 'low', 'high')),
                 values_to = 'value', names_to = c('Biomass', 'measure', 'stat'),
                 names_pattern = '(epil|fila|mu)_(.*)_(med|high|low)')  %>%
    filter(Biomass != 'mu') %>%
    left_join(bloom_df, by = c("site", "year")) %>%
    mutate(value = if_else(value > 10000, NA, value))

p1 <- bm_fila %>%
    filter(measure == 'turn') %>%
    ggplot(aes(x=value, group = Biomass, fill = Biomass)) +
    geom_density(adjust=1.5, alpha=.4) +
    xlim(0.1, 100)+
    scale_fill_manual(values = c('#1B9EC9', '#97BB43'))+
    scale_x_log10(limits = c(1, 1500))+
    ylab('Density')+
    xlab('Turnover time (d)')+
    theme_classic()+
    theme(legend.position = 'none',
          panel.border = element_rect(fill = NA))
p1 <- p1 + annotate(geom = 'text', x = 12, y = 0.25,
                    label="Epilithic", col = '#1B9EC9')
p1 <- p1 + annotate(geom = 'text', x = 32, y = 0.06,
                    label="Filamentous", col = '#97BB43')
p2 <- bm_fila %>%
    ggplot(aes(frac_fila, value))+
    xlab('Filamentous fraction of \ntotal biomass (%)')+
    ylab('Turnover time (d)')+
    geom_point(aes(pch = factor(fila_bloom), col = factor(fila_bloom)),
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
coeff = 200
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
    bm_met %>%
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
    ggplot(aes(date, prodgCd_med/gCm2, col = biomass))+
    geom_area(aes(date, gCm2/coeff, fill = biomass), color = NA, alpha = 0.4)+
    geom_line(size = 1.2)+
    geom_line(aes(y = light/2.5, lty = Light), col = 'grey20') +
    facet_grid(site~year, scales = 'free_x', ) +
    geom_text(data = ann_text2, aes(label = trophic), col = 'black') +
    geom_text(data = ann_text2.1, aes(label = sitename), col = 'black') +
    scale_y_continuous(
        expand = expand_scale(mult = c(0.05, 0.1)),
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

bm_sum <-  bm_met %>%
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

p1 <- bm_fila %>%
    filter(measure == 'turn') %>%
    ggplot(aes(x=value, group = Biomass, fill = Biomass)) +
    geom_density(adjust=1.5, alpha=.4) +
    xlim(0.1, 100)+
    scale_fill_manual(values = c('#1B9EC9', '#97BB43'))+
    scale_x_log10(limits = c(1, 1500))+
    ylab('Density')+
    xlab('Turnover time (d)')+
    theme_classic()+
    theme(legend.position = 'none',
          panel.border = element_rect(fill = NA))
p1 <- p1 + annotate(geom = 'text', x = 25, y = 3,
                    label="Epilithic", col = '#1B9EC9')
p1 <- p1 + annotate(geom = 'text', x = 80, y = 0.8,
                    label="Filamentous", col = '#97BB43')
png('figures/biomass_cumulative_and_turnover.png', width = 6.5, height = 3.5,
    units = 'in', res = 300)
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



