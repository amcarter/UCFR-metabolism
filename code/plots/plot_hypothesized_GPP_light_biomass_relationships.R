library(tidyverse)
library(viridis)

dat <- read_csv('data/model_fits/biomass_metab_model_data.csv') %>%
    mutate(biomass = fila_chla_mgm2_fit + epil_chla_mgm2_fit)
# baseline, light only
pdat <- data.frame(light = seq(min(dat$light, na.rm = T),
                               max(dat$light, na.rm = T), length.out = 90),
                   biomass = rep(seq(min(dat$biomass, na.rm = T),
                                 max(dat$biomass, na.rm = T),
                                 length.out = 15), 6))
mod <- summary(lme4::lmer(GPP~light + (1|site), dat))$coefficients
pdat$linlight <- mod[1,1] + mod[2,1] * pdat$light

mod <- summary(lme4::lmer(log(GPP)~log(light) + (1|site), dat))$coefficients
pdat$loglight <- exp(mod[1,1] + mod[2,1] * log(pdat$light))

plot(pdat$light, pdat$linlight)
lines(pdat$light, pdat$loglight)


# biomass model
mod <- summary(lme4::lmer(GPP~light + biomass + (1|site), dat))$coefficients
pdat$linlightbio <- mod[1,1] + mod[2,1] * pdat$light + mod[3,1] * pdat$biomass


mod <- summary(lme4::lmer(log(GPP)~log(light) + log(biomass) + (1|site), dat))$coefficients
pdat$loglightbio <- exp(mod[1,1] + mod[2,1] * log(pdat$light) +
                            mod[3,1] * log(pdat$biomass))


ggplot(pdat, aes(light, linlightbio, col = biomass, group = factor(biomass)))+
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()
ggplot(pdat, aes(light, loglightbio, col = biomass, group = factor(biomass)))+
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()


mod <- summary(lme4::lmer(GPP~light * biomass + (1|site), dat))$coefficients
pdat$linlightbioint <- mod[1,1] + mod[2,1] * pdat$light + mod[3,1] * pdat$biomass +
    mod[4,1] * pdat$light * pdat$biomass


mod <- summary(lme4::lmer(log(GPP)~log(light) * log(biomass) + (1|site), dat))$coefficients
pdat$loglightbioint <- exp(mod[1,1] + mod[2,1] * log(pdat$light) +
                            mod[3,1] * log(pdat$biomass) +
                            mod[4,1]/9 * log(pdat$biomass)*log(pdat$light))


ggplot(pdat, aes(light, linlightbioint, col = biomass, group = factor(biomass)))+
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()
ggplot(pdat, aes(light, loglightbioint, col = biomass, group = factor(biomass)))+
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()

# interaction only
dat$interaction <- dat$light * dat$biomass
dat$loginteraction <- log(dat$light) * log(dat$biomass)
mod <- summary(lme4::lmer(GPP~ interaction + (1|site), dat))$coefficients
pdat$linlightbiointonly <- mod[1,1] + mod[2,1] * pdat$light * pdat$biomass


mod <- summary(lme4::lmer(log(GPP)~log(interaction) + (1|site), dat))$coefficients
pdat$loglightbiointonly <- exp(mod[1,1] + mod[2,1] * log(pdat$light* pdat$biomass))


# ggplot(pdat, aes(light, linlightbiointonly, col = biomass, group = factor(biomass)))+
#     geom_line(size = 1.5) +
#     scale_color_viridis(option = 'D') +
#     theme_bw()
ggplot(pdat, aes(light, loglightbiointonly, col = biomass, group = factor(biomass)))+
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()

pdat %>%
    pivot_longer(cols = starts_with(c('log', 'lin')),
                 values_to = 'GPP',
                 names_to = 'model') %>%
    filter(model %in% c('linlightbio','linlightbioint',
                        'loglightbio', 'loglightbioint')) %>%
    mutate(model = case_when(model == 'linlightbio' ~ 'linear model',
                             model == 'linlightbioint' ~ 'linear model with interaction',
                             model == 'loglightbio' ~ 'log model',
                             model == 'loglightbioint'~ 'log model with interaction')) %>%
    ggplot(aes(light, GPP, group = factor(biomass), col = biomass)) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D',
                        name = expression(paste('Chl a \n(mg', m^-2,')'))) +
    facet_wrap(.~model)+
    theme_bw()+
    xlab('Relative light')+
    ylab(expression(paste('Productivity (g ',O[2],  m^-2, d^-1, ')')))


png('figures/model_examples2.png', width = 6.5, height = 5, units = 'in', res = 300)

pdat %>%
    pivot_longer(cols = starts_with(c('log', 'lin')),
                 values_to = 'GPP',
                 names_to = 'model') %>%
    filter(model %in% c('loglight','loglightbioint',
                        'loglightbio', 'loglightbiointonly')) %>%
    mutate(model = factor(model, levels = c('loglight','loglightbio',
                                            'loglightbiointonly',
                                            'loglightbioint')),
           model = case_when(model == 'loglight' ~ '0. Baseline model (Light)',
                             model == 'loglightbiointonly' ~ '3. Light \U00D7 Biomass',
                             model == 'loglightbio' ~ '2. Light + Biomass',
                             model == 'loglightbioint'~ '4. Light + Biomass + Light \U00D7 Biomass'),
           cover = case_when(model == '0. Baseline model (Light)' ~ GPP,
                             TRUE ~ NA_real_)) %>%
    ggplot(aes(light, GPP, group = factor(biomass), col = biomass)) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D',
                        name = expression(paste('Chl a \n(mg', m^-2,')'))) +
    geom_line(aes(y = cover), col = 'grey', size = 1.5) +
    facet_wrap(.~model)+
    theme_bw()+
    xlab('Relative light')+
    ylab(expression(paste('Productivity (g ',O[2],  m^-2, d^-1, ')'))) +
    xlim(0,1) + ylim(0,16.5)

dev.off()

