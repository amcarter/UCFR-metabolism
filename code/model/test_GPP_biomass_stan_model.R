# testing GPP biomass model

# next steps:
# - add grouping variable by site
# - add uncertainty draws for biomass
# - decide how to scale data

# remove.packages(c("rstan", "StanHeaders"))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
dat <- read_csv('data/biomass_metab_model_data.csv')

# set params
n <- 600
gamma <- c(1,3,1)       # population level coefficients
tau <- 1                # standard deviation of regression coefficients
sites <- seq(1:6)       # unique sites
K <- length(gamma)      # number of parameters
S <- length(sites)      # number of sites
ss <- rep(sites, n/S)   # site of each observation

light <- abs(rnorm(n = n))
light <- light/max(light)
biomass <- abs(rnorm(n = n))

X <- matrix(c(rep(1, n), light, biomass), ncol = 3) # model matrix
sigma_epsilon <- 0.5

beta <- matrix(NA, ncol = K, nrow = S)
for(i in 1:S){
    beta[i,] <- rnorm(3, gamma, tau)
}

# simulate response vector
mu <- vector(mode = 'double', length = n)

for(i in 1:n){
    mu[i] <- X[i,]%*%beta[ss[i],]
}

P <- vector(mode = "double", length = n)
P <- rnorm(n, mean = mu, sd = sigma_epsilon)

sim <- data.frame(site = ss, light = light, biomass = biomass, GPP = P)
ggplot(sim, aes(biomass, GPP, col = factor(site))) +
    geom_point()
# compile data for stan
datlist <- list(
    N = n, K = K, S = S,
    ss = ss, X = X,
    P = P
)

# compile stan model
P_mod <- stan_model("code/model/GPP_biomass_model.stan")

# fit the model
mfit <- sampling(
    P_mod,
    data = datlist
)


plot(mfit, pars = c('beta[1,1]','beta[2,1]','beta[3,1]',
                    'beta[4,1]','beta[5,1]','beta[6,1]'))
plot(mfit, pars = c('beta[1,2]','beta[2,2]','beta[3,2]',
                    'beta[4,2]','beta[5,2]','beta[6,2]'))
plot(mfit, pars = c('beta[1,3]','beta[2,3]','beta[3,3]',
                    'beta[4,3]','beta[5,3]','beta[6,3]'))
traceplot(mfit, pars = 'beta')

# test on dataset:
dat <- dplyr::filter(dat, !is.na(GPP), !is.na(light), !is.na(fila_chla_mgm2_fit))

X_dat <- matrix(c(rep(1, nrow(dat)), dat$light, dat$fila_chla_mgm2_fit), ncol = 3)
ss_dat <- as.numeric(factor(dat$site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

datlist <- list(
    N = nrow(dat),
    K = 3, S = 6,
    ss = ss_dat, X = X_dat,
    P = dat$GPP
)

dfit <- sampling(
    P_mod,
    data = datlist
)

plot(dfit, pars = c('beta', 'sigma'))
traceplot(dfit)
