# testing GPP biomass model

# next steps:
# - add grouping variable by site
# - add uncertainty draws for biomass
# - decide how to scale data

library(rstan)
dat <- read_csv('data/biomass_metab_model_data.csv')

# set params
n <- 100
b0 <- 1
b1 <- 3
b2 <- 1
light <- abs(rnorm(n = n))
light <- light/max(light)

biomass <- abs(rnorm(n = n))
sigma_epsilon <- 0.5

# simulate response vector
P <- vector(mode = "double", length = n)
P <- rnorm(n, mean = b0 + b1 * light + b2 * biomass, sd = sigma_epsilon)

# compile data for stan
datlist <- list(
    N = n,
    light = light,
    biomass = biomass,
    P = P
)

# compile stan model
P_mod <- stan_model("code/model/GPP_biomass_model.stan")

# fit the model
mfit <- sampling(
    P_mod,
    data = datlist
)


plot(mfit)
traceplot(mfit)

# test on dataset:
dat <- filter(dat, !is.na(GPP), !is.na(light), !is.na(epil_chla_mgm2_fit))
datlist <- list(
    N = nrow(dat),
    light = dat$light,
    biomass = dat$epil_chla_mgm2_fit,
    P = dat$GPP
)

dfit <- sampling(
    P_mod,
    data = datlist
)

plot(dfit)
traceplot(dfit)
