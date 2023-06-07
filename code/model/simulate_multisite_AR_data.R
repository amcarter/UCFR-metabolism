# Simulate data from a hierarcical AR process:
library(mvtnorm)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

# define parameters:
S <- 6    # number of sites
N <- 200  # number of observations per site
K <- 3    # number of covariates

phi <- 0.5
beta <- rnorm(K+1)  # population level effects
tau <- 0.1
beta0 <- rnorm(S, beta[1], tau) # random intercept centered around gamma
sigma_proc <- 1
sigma_obs <- 0.1

new_year <- 101                 # index where a new year starts

X <- array(rnorm(N*K*S), dim = c(S, N, K))
# X <- matrix(rnorm(N*K), nrow = N)

theta <- vector()
theta[1] <- sigma_proc/sqrt(1-phi)
for(i in 2:N){
    if(i == new_year){
        theta[i] <- sigma_proc/sqrt(1-phi)
    }else{
        theta[i] = phi * theta[i-1] + rnorm(1,0,sigma_proc)
    }
}

mu <- matrix(NA, nrow = N, ncol = S)
for(s in 1:S){
    mu[,s] <- X[s,,] %*% beta[2:(K+1)] + theta + beta0[s]
}

y <- matrix(NA, nrow = N, ncol = S)

for(t in 1:N){
    y[t,] = mvtnorm::rmvnorm(n = 1, mu[t,], sigma_obs * diag(S))
}

plot(theta, ylim = c(min(y), max(y)), type = 'l', col = 2)
for(s in 1:S){
    lines(y[,s])
}

N_mis <- 5
ii_mis <- vector(mode = 'integer')
ii_mis = sample(1:N*S, N_mis)
y_mis <- y[ii_mis]

ar_mod <- stan_model("code/model/stan_code/AR1_multisite.stan")

datlist <- list(
    S = S,
    N = N, K = K,
    X = X, y = y,
    new_year = new_year,
    N_mis = N_mis,
    ii_mis = ii_mis)

sim_fit <- sampling(ar_mod,
                    datlist)

print(sim_fit, pars = c('phi', 'beta', 'sdp', 'sdo', 'beta0', 'tau', 'y_mis'))
print(paste0('phi = ', phi, '; beta = ', beta, '; sdp = ', sigma_proc,
             '; sdo = ', sigma_obs, '; beta0 = ', beta0, '; tau = ', tau,
             '; y_mis = ', y_mis))

# Test out on real dataset: ####
dd <- read_csv('data/model_fits/biomass_metab_model_data.csv')
sites <- factor(unique(dd$site))
S = length(sites)
K = 1
y = dd %>%
    mutate(log_GPP = log(GPP)) %>%
    pivot_wider(id_cols = date, names_from = site, values_from = log_GPP) %>%
    arrange(date) #%>% tail()
dates <- c(min(dd$date[!is.na(dd$GPP)]),
           max(dd$date[!is.na(dd$GPP)& dd$year == 2020]),
           min(dd$date[!is.na(dd$GPP)& dd$year == 2021]),
           max(dd$date[!is.na(dd$GPP)]))
dateseq <- data.frame(date = c(seq(dates[1],dates[2], by = 'day'),
                               seq(dates[3], dates[4], by = 'day')))
N <- nrow(dateseq)
new_year = which(dateseq$date == dates[3])
y <- left_join(dateseq, y) %>%
    select(-date) %>%
    as.matrix()

X <- array(NA_real_, dim = c(S, N, K))

for(s in 1:S){
    X[s,,] <- dd %>%
        filter(site == sites[s]) %>%
        right_join(dateseq, by = 'date') %>%
        pull(light)
}


X_fill <- apply(X, 2, function(x) mean(x, na.rm = T))
for(s in 1:S){
    w <- which(is.na(X[s,,]))
    X[s,w,] <- X_fill[w]
}
X2 <- array(NA_real_, dim = c(S,210, K))
X2[,,1] <- X[,1:210,]
X <- X2

y <- y[1:210,]
ii_mis <- which(is.na(y))
N_mis <- length(ii_mis)

y[ii_mis] <- 0
N = nrow(y)

datlist <- list(
    S = S,
    N = N, K = K,
    X = X, y = y,
    new_year = new_year,
    N_mis = N_mis,
    ii_mis = ii_mis)

mod_fit <- sampling(ar_mod,
                    datlist)

print(mod_fit, pars = c('phi', 'beta', 'sdp', 'sdo', 'beta0', 'tau', 'y_mis'))
pairs(mod_fit, pars = c('phi', 'beta', 'sdp', 'sdo', 'tau'))

yhat <- rstan::extract(mod_fit, pars = 'y_hat')
preds <- data.frame(apply(yhat$y_hat, c(2,3), mean))
colnames(preds) <- sites
preds %>%
    mutate(date = dateseq$date[1:210])%>%
    pivot_longer(cols = -date,
                 names_to = 'site', values_to = 'GPP_pred') %>%
    left_join(select(dd, date, site, GPP)) %>%
    mutate(GPP_pred = exp(GPP_pred),
           year = year(date)) %>%
    ggplot(aes(date, GPP)) +
    geom_line() +
    geom_line(aes(y = GPP_pred), col = 2)+
    facet_grid(site~year, scales = 'free_x')

ggplot(dd, aes(light, GPP, col = site)) +
    geom_point()

