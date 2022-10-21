# fit GAMS to biomass data
# 10/2022

library(tidyverse)
library(bayesGAM)


biomass <- read_csv('data/biomass_data/biomass_working_data.csv') %>%
    mutate(site = case_when(site == 'BG' ~ 'BM',
                            site == 'WS' ~ 'PL',
                            TRUE ~ site),
           doy = as.numeric(format(date, '%j')))%>%
    mutate(site = factor(site, levels = c('PL', 'DL', 'GR', 'GC', 'BM', 'BN')))

glimpse(biomass)


# look at the data
biomass %>%
    select(date, site, ends_with('.om.area.g.m2')) %>%
    select(-cpom.om.area.g.m2)%>%
    filter(!is.na(site))%>%
    pivot_longer(cols = -c('date', 'site'), names_to = 'bm_category',
                 values_to = 'g_m2') %>%
    mutate(doy = format(date, "%j"),
           year = as.factor(substr(date, 1,4))) %>%
    ggplot(aes(doy, g_m2, col = year)) +
    geom_point() +
    facet_grid(site~bm_category, scales = 'free_y')


# There is a lot of missing data for the epiphyton and macrophyte categories

# test bayesGAM on the 2020 epilition data from deer lodge

dl_2020 <- biomass %>%
    mutate(year = as.factor(year))%>%
    filter(site == 'DL', year == 2020)

bg <- bayesGAM::bayesGAM(epil.om.area.g.m2 ~ doy, data = dl_2020)

preds <- data.frame(doy = seq(min(dl_2020$doy), max(dl_2020$doy), by = 1))

post <- bayesGAM::predict(bg, preds, draws = 200)

y <- data.frame(median = apply(post@pp$ypred, 2, FUN = median),
                y_2.5 = apply(post@pp$ypred, 2, FUN = quantile, probs = 0.025),
                y_97.5 = apply(post@pp$ypred, 2, FUN = quantile, probs = 0.975))

preds <- bind_cols(preds, y)
# try a non bayesian approach
library(mgcv)
fg <- gam(epil.om.area.g.m2 ~ s(doy),# by = year),
          data = dl_2020, method = "REML")

par(mfrow = c(2,1),
    mar = c(2,5,1,2), oma = c(3,0,1,0))
ylim = range(c(preds$y_2.5, preds$y_97.5))
plot(dl_2020$doy, dl_2020$epil.om.area.g.m2, ylim = ylim, ylab = 'BayesGAM fit')
lines(preds$doy, preds$median)
polygon(x = c(preds$doy, rev(preds$doy)),
        y = c(preds$y_2.5, rev(preds$y_97.5)),
        col = alpha('grey', .5), border = NA)

# par(mfrow = c(2,1))
plot(fg,
     residuals = TRUE,
     rug = FALSE,
     shade = TRUE,
     shift = coef(fg)[1],
     seWithMean = TRUE,
     ylab = 'mgcvGAM fit')

points(dl_2020$doy, dl_2020$epil.om.area.g.m2)

preds_2 <- mgcv::predict.gam(fg, preds, se.fit = TRUE)

preds$f_med <- preds_2$fit
preds$f_se <- preds_2$se.fit
ylim = range(c(preds$f_med + 1.96 * sqrt(nrow(dl_2020)) *preds$f_se,
               preds$f_med - 1.96* sqrt(nrow(dl_2020)) * preds$f_se))
plot(dl_2020$doy, dl_2020$epil.om.area.g.m2, ylim = ylim, ylab = 'mgcvGAM fit')

lines(preds$doy, preds$f_med)

summary(fg)
polygon(x = c(preds$doy, rev(preds$doy)),
        y = c(preds$f_med + 1.96 * sqrt(nrow(dl_2020)) *preds$f_se,
              rev(preds$f_med - 1.96 * sqrt(nrow(dl_2020)) * preds$f_se)),
        col = alpha('grey', .5), border = NA)
polygon(x = c(preds$doy, rev(preds$doy)),
        y = c(preds$f_med + preds$f_se,
              rev(preds$f_med -  preds$f_se)),
        col = alpha('grey', .5), border = NA)

# wy is this working better?
# fitting it with just time as a predictor is clearly not working very well
acf(preds$f_med)
