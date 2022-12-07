# linear no interaction
beta0 <- -2
beta1 <- 2.7
beta2 <- 1.5
biomass = seq(1, 10, by = 0.5)
PI_1 <- data.frame()

for(B in biomass){
    df <- data.frame(light = seq(0,1, by = 0.01),
                     biomass = B)

    df$GPP = beta0 + beta1 * df$light + beta2 * B

    PI_1 <- bind_rows(PI_1, df)
}

A <- ggplot(PI_1, aes(light, GPP, col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()

# linear + interaction
beta0 <- -2
beta1 <- 1.07
beta2 <- 0.6
beta3 <- 1.06
biomass = seq(1, 10, by = 0.5)
PI_2 <- data.frame()

for(B in biomass){
    df <- data.frame(light = seq(0,1, by = 0.01),
                     biomass = B)

    df$GPP = beta0 + beta1 * df$light + beta2 * B + beta3 * df$light * B

    PI_2 <- bind_rows(PI_2, df)
}

BB <- ggplot(PI_2, aes(light, GPP, col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()

#linear interaction term only
beta0 <- 0
beta1 <- 1.5
biomass = seq(1, 10, by = 0.5)
PI_3 <- data.frame()

for(B in biomass){
    df <- data.frame(light = seq(0,1, by = 0.01),
                     biomass = B)

    df$GPP = beta0 + beta1 * df$light * B

    PI_3 <- bind_rows(PI_3, df)
}

C <- ggplot(PI_3, aes(light, GPP, col = biomass, group = factor(biomass))) +
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

D <- ggplot(PI_4, aes(light, GPP, col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    theme_bw()


ggpubr::ggarrange(A,BB,C,D, common.legend = TRUE, ncol = 2, nrow = 2)

PI_1$model = 'linear model'
PI_2$model = 'linear model with interaction'
PI_3$model = 'linear model only interaction'
PI_4$model = 'PI curve'

PI <- bind_rows(PI_1, PI_2, PI_3, PI_4) %>%
    mutate(model = factor(model,
                          levels = c('linear model',
                                     'linear model with interaction',
                                     'linear model only interaction',
                                     'PI curve')))

ggplot(PI, aes(light, GPP, col = biomass, group = factor(biomass))) +
    geom_line(size = 1.5) +
    scale_color_viridis(option = 'D') +
    facet_wrap(.~model, ncol = 2)+
    theme_bw()
