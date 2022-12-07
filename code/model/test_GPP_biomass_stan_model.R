# testing GPP biomass model

# remove.packages(c("rstan", "StanHeaders"))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
options(mc.cores = parallel::detectCores())
library(tidyverse)
dd <- read_csv('data/biomass_metab_model_data.csv')
dd <- dd %>%
    mutate(across(starts_with('epil_chla'),
                  function(x) x/max(dd$epil_chla_mgm2_fit, na.rm = T)),
           across(starts_with('epil_gm'),
                  function(x) x/max(dd$epil_gm2_fit, na.rm = T)),
           across(starts_with('fila_chla'),
                  function(x) x/max(dd$fila_chla_mgm2_fit, na.rm = T)),
           across(starts_with('fila_gm'),
                  function(x) x/max(dd$fila_gm2_fit, na.rm = T)))
# set params
n <- 600
gamma <- c(1,3,1)       # population level coefficients
tau <- 1                # standard deviation of regression coefficients
sites <- seq(1:6)       # unique sites
K <- length(gamma)-1    # number of parameters
S <- length(sites)      # number of sites
ss <- rep(sites, n/S)   # site of each observation

light <- abs(rnorm(n = n))
light <- light/max(light)
biomass <- abs(rnorm(n = n))

X <- matrix(c(light, biomass), ncol = K) # model matrix
sigma_epsilon <- 0.5

beta <- vector(mode = 'double', length = S)
for(i in 1:S){
    beta[i] <- rnorm(1, gamma[1], tau)
}

# simulate response vector
mu <- vector(mode = 'double', length = n)

for(i in 1:n){
    mu[i] <- beta[ss[i]] + X[i,]%*%gamma[2:3]
}

P <- vector(mode = "double", length = n)
P <- rnorm(n, mean = mu, sd = sigma_epsilon)

sim <- data.frame(site = ss, light = light, biomass = biomass, GPP = P)
ggplot(sim, aes(light, GPP, col = factor(site))) +
    geom_point()
# compile data for stan
datlist <- list(
    N = n, K = K, S = S,
    ss = ss, X = X,
    P = P
)

# compile stan model
P_mod <- stan_model("code/model/stan_code/GPP_biomass_model.stan")

# fit the model
mfit <- sampling(
    P_mod,
    data = datlist
)

shinystan::launch_shinystan(mfit)
plot(mfit)
# plot(mfit, pars = c('beta[1,1]','beta[2,1]','beta[3,1]',
#                     'beta[4,1]','beta[5,1]','beta[6,1]'))
# plot(mfit, pars = c('beta[1,2]','beta[2,2]','beta[3,2]',
#                     'beta[4,2]','beta[5,2]','beta[6,2]'))
# plot(mfit, pars = c('beta[1,3]','beta[2,3]','beta[3,3]',
#                     'beta[4,3]','beta[5,3]','beta[6,3]'))
traceplot(mfit)
pairs(mfit, pars = )
# test on dataset:
# plot the data:
dd %>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')),
           biomass = fila_gm2_fit + epil_gm2_fit,
           biomass_chla = fila_chla_mgm2_fit + epil_chla_mgm2_fit) %>%
    pivot_longer(cols = c('biomass', 'biomass_chla', 'light'),
                 values_to = 'value', names_to = 'driver') %>%
    ggplot(aes(value, GPP, col = site)) +
    geom_point()+
    geom_smooth(method = 'lm', se = FALSE)+
    facet_wrap(.~driver, ncol = 3, scales = 'free_x')
b_plot <- ggplot(dd, aes(epil_gm2_fit + fila_gm2_fit, GPP, col = factor(site)))+
    geom_point()
bc_plot <- ggplot(dd, aes(epil_chla_mgm2_fit + fila_chla_mgm2_fit,
                          GPP, col = factor(site)))+
    geom_point()
ggpubr::ggarrange(l_plot, b_plot, bc_plot, nrow = 1, common.legend = TRUE)

dat <- dplyr::filter(dd, !is.na(GPP), !is.na(light), !is.na(fila_chla_mgm2_fit))
X_dat <- matrix(c(dat$light, dat$fila_chla_mgm2_fit), ncol = 2)
# X_dat <- matrix(c(rep(1, nrow(dat)), dat$light, dat$fila_chla_mgm2_fit), ncol = 3)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 2, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)
dfit <- sampling(
    P_mod,
    data = datlist
)
plot(dfit)
plot(dfit, pars = c('beta', 'sigma'))
traceplot(dfit)
pairs(dfit)

shinystan::launch_shinystan(dfit)

# fit to different combinations of data
dat <- dplyr::filter(dd, !is.na(GPP), !is.na(light), !is.na(fila_chla_mgm2_fit))
X_dat <- matrix(c(dat$light, dat$fila_chla_mgm2_fit), ncol = 2)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 2, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)

fit_fila_chla <- sampling(
    P_mod,
    data = datlist
)

dat <- dplyr::filter(dd, !is.na(GPP), !is.na(light), !is.na(fila_gm2_fit))
X_dat <- matrix(c(dat$light, dat$fila_gm2_fit), ncol = 2)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 2, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)

fit_fila <- sampling(
    P_mod,
    data = datlist
)

dat <- dplyr::filter(dd, !is.na(GPP), !is.na(light), !is.na(epil_chla_mgm2_fit))
X_dat <- matrix(c(dat$light, dat$epil_chla_mgm2_fit), ncol = 2)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 2, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)

fit_epil_chla <- sampling(
    P_mod,
    data = datlist
)

dat <- dplyr::filter(dd, !is.na(GPP), !is.na(light), !is.na(epil_gm2_fit))
X_dat <- matrix(c(dat$light, dat$epil_gm2_fit), ncol = 2)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 2, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)

fit_epil <- sampling(
    P_mod,
    data = datlist
)

dat <- dplyr::filter(dd, !is.na(GPP), !is.na(light),
                     !is.na(epil_gm2_fit), !is.na(fila_gm2_fit))
X_dat <- matrix(c(dat$light, dat$epil_gm2_fit, dat$fila_gm2_fit), ncol = 3)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 3, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)

fit_epil_fila <- sampling(
    P_mod,
    data = datlist
)

dat <- dplyr::filter(dd, !is.na(GPP), !is.na(light),
                     !is.na(epil_chla_mgm2_fit), !is.na(fila_chla_mgm2_fit))
X_dat <- matrix(c(dat$light, dat$epil_chla_mgm2_fit, dat$fila_chla_mgm2_fit), ncol = 3)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 3, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)

fit_epil_fila_chla <- sampling(
    P_mod,
    data = datlist
)

print(fit_epil_fila_chla)
plot(fit_epil_fila_chla)
plot(fit_epil_chla, pars = c('gamma', 'tau', 'beta', 'sigma'))
pairs(fit_epil_chla, pars = c('gamma', 'tau', 'beta', 'sigma'))
plot(fit_epil, pars = c('gamma', 'tau', 'beta', 'sigma'))
pairs(fit_epil, pars = c('gamma', 'tau', 'beta', 'sigma'))
plot(fit_fila, pars = c('gamma', 'tau', 'beta', 'sigma'))
pairs(fit_fila, pars = c('gamma', 'tau', 'beta', 'sigma'))
plot(fit_fila_chla, pars = c('gamma', 'tau', 'beta', 'sigma'))+
    xlim(-4, 14)
pairs(fit_fila, pars = c('gamma', 'tau', 'beta', 'sigma'))
plot(fit_epil_fila_chla, pars = c('gamma', 'tau', 'beta', 'sigma')) +
    xlim(-4, 14)
plot(fit_epil_fila, pars = c('gamma', 'tau', 'beta', 'sigma'), xlim = c(-4, 14))

