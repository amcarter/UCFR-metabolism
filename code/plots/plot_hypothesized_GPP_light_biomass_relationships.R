library(tidyverse)
library(viridis)

# baseline, light only
beta1 <- 15

BL <- data.frame(light = seq(0,1, by = 0.01))
BL$GPP = beta1 * df$light


L <- ggplot(BL, aes(light, GPP)) +
    geom_line(size = 1.5) +
    theme_bw()
# linear no interaction
beta0 <- 2.5
beta1 <- 0.8
beta2 <- 0.05
biomass = seq(1, 10, by = 0.5)
PI_1 <- data.frame()

for(B in biomass){
    df <- data.frame(light = seq(0.13,1, by = 0.01),
                     biomass = B)

    df$GPP = beta0 + beta1 * log(df$light) + beta2 *(B)

    PI_1 <- bind_rows(PI_1, df)
}

# A <-
ggplot(PI_1, aes((light), exp(GPP), col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()

# linear + interaction
beta0 <- 0
beta1 <- 0.2
beta2 <- 0.1
beta3 <- 0.2
biomass = seq(1, 10, by = 0.5)
PI_2 <- data.frame()

for(B in biomass){
    df <- data.frame(light = seq(0,1, by = 0.01),
                     biomass = B)

    df$GPP = beta0 + beta1 * df$light + beta2 * B + beta3 * df$light * B

    PI_2 <- bind_rows(PI_2, df)
}

BB <-
ggplot(PI_2, aes(light, exp(GPP), col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()

#linear interaction term only
beta0 <- 0
beta1 <- 0.3
biomass = seq(1, 10, by = 0.5)
PI_3 <- data.frame()

for(B in biomass){
    df <- data.frame(light = seq(0,1, by = 0.01),
                     biomass = B)

    df$GPP = beta0 + beta1 * df$light * B

    PI_3 <- bind_rows(PI_3, df)
}

C <-
    ggplot(PI_3, aes(light, exp(GPP), col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()


# PI curve
alpha <- 5
Pmax <- 1.7
biomass = seq(1, 10, length.out = 20)
PI_4 <- data.frame()

for(B in biomass){
    df <- data.frame(light = seq(0,1, by = 0.01),
                     biomass = B)

    df$GPP = Pmax * B * tanh(alpha * df$light/Pmax)

    PI_4 <- bind_rows(PI_4, df)
}

D <-
    ggplot(PI_4, aes(light, GPP, col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()




ggpubr::ggarrange(A,BB,C,D, common.legend = TRUE, ncol = 2, nrow = 2)

PI_1$model = '1. Linear model'
PI_2$model = '2. Linear model with interaction'
PI_3$model = '3. Linear model only interaction'
PI_4$model = '4. PI curve'

PI <- bind_rows(PI_1, PI_2, PI_3, PI_4) %>%
    mutate(model = factor(model,
                          levels = c('1. Linear model',
                                     '2. Linear model with interaction',
                                     '3. Linear model only interaction',
                                     '4. PI curve')))
png('figures/model_examples.png', width = 5.5, height = 4, units = 'in', res = 300)
PI %>% rename(Biomass = biomass) %>%
ggplot( aes(light, GPP, col = Biomass, group = factor(Biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    facet_wrap(.~model, ncol = 2)+
    theme_bw()+
    ylab(expression(paste('Productivity (g ',O[2],  m^-2, d^-1, ')')))+
    xlab('Relative light')
dev.off()


BL$model = '0. Baseline model'
BL1 <- BL %>% mutate(model = 'a')
BL2 <- BL %>% mutate(model = 'b')
BL3 <- BL %>% mutate(model = 'c')
BL <- bind_rows(BL, BL2, BL1, BL3) %>%
    mutate(model = factor(model,
                          levels = c('a', 'b', '0. Baseline model', 'c')))
png('figures/baseline_model_examples.png', width = 4.5, height = 4, units = 'in', res = 300)
BL %>% ggplot( aes(light, GPP)) +
    geom_line(size = 1.5) +
    facet_wrap(.~model, ncol = 2)+
    theme_bw()+
    ylab(expression(paste('Productivity (g ',O[2],  m^-2, d^-1, ')')))+
    xlab('Relative light')
dev.off()
