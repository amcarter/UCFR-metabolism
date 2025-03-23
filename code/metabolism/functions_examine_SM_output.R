# Functions to check stream metabolizer outputs
# A Carter
# 8/2022

library(tidyverse)
library(streamMetabolizer)

# inputs for function development:
# Stream Metabolizer input dataframe:
# dat <- read_csv('data/prepared_data/PL_2020.csv')
# Stream Metabolizer fit object:
# fit <- readRDS('data/metab_fits/PL_2020_kn_oipi.rds')

# Extract data from fit ####
# get metabolism and K600 fit data with Rhats and CIs
extract_metab <- function(fit, sitecode = NA, bad_days = NULL, mle = FALSE){ # bad days option allows for flagging bad DO fits
    if(mle){
        met <- fit@fit %>%
            select(date,
                   GPP = GPP.daily, GPP.sd = GPP.daily.sd,
                   ER = ER.daily, ER.sd = ER.daily.sd) %>%
            mutate(GPP.lower = GPP - 1.96*GPP.sd,
                   GPP.upper = GPP + 1.96*GPP.sd,
                   ER.lower = ER - 1.96*ER.sd,
                   ER.upper = ER + 1.96*ER.sd)
        met <- fit@data_daily %>% rename(K600 = K600.daily) %>%
            left_join(met, by = 'date')
    } else{

    met <- fit@fit$daily %>%
        select(date,
               GPP = GPP_daily_50pct, GPP.lower = GPP_daily_2.5pct,
               GPP.upper = GPP_daily_97.5pct, ER = ER_daily_50pct,
               ER.lower = ER_daily_2.5pct, ER.upper = ER_daily_97.5pct,
               K600 = K600_daily_50pct, K600.lower = K600_daily_2.5pct,
               K600.upper = K600_daily_97.5pct, GPP_Rhat = GPP_daily_Rhat,
               ER_Rhat = ER_daily_Rhat, K600_Rhat = K600_daily_Rhat)
    }

    m <- predict_metab(fit) %>%
        select(date, msgs.fit, warnings, errors)
    met <- left_join(met, m, by = 'date')

    if(!is.na(sitecode)){
        met <- mutate(met, site = sitecode) %>%
            relocate(site)
    }

    if(!is.null(bad_days)){
        met$DO_fit = 'good'
        if(length(bad_days)>0){
            met$DO_fit[met$date %in% bad_days]<- 'bad'
        }
    }

    return(met)
}

# Visually examine fits: ####
# functions that might be useful in addition to the built in functions:
# plot_DO_preds(fit, y_var = c('pctsat', 'conc', 'ddodt'), style = 'dygraphs')
# plot_metab_preds(fit)
# plot Rhats
plot_Rhats <- function(fit){
    rh <- fit@fit$daily %>%
        select(date, ends_with('daily_Rhat')) %>%
        rename_with(function(x) sub('_daily_', '_', x),
                    ends_with('Rhat'))
    ylim = range(c(rh$K600_Rhat, rh$GPP_Rhat, rh$ER_Rhat, 1.05),
                 na.rm = T)
    plot(x = rh$date, y = rh$K600_Rhat,
         ylab = "Rhat (convergence metric)", xlab = "date",
         type = "l", lwd = 1.5, ylim = ylim)
    lines(rh$date, rh$GPP_Rhat, col = "forestgreen", lwd = 1.5)
    lines(rh$date, rh$ER_Rhat, col = "sienna", lwd = 1.5)
    abline(h = 1.05, lty = 2, col = "red")
    mtext("Rhat below 1.05 is good", 3, 0, adj = 0, cex = .8)
    legend("topleft",
           legend = c("K600", "GPP", "ER"),
           col = c(1, "forestgreen", "sienna"),
           lty = 1, bty = "n", lwd = 1.5)
}

get_bad_Rhats <- function(fit, threshold = 1.05, vars = c('GPP', 'ER', 'K600')){
    rh <- fit@fit$daily %>%
        select(date, ends_with('daily_Rhat')) %>%
        rename_with(function(x) sub('_daily_Rhat', '', x),
                    ends_with('Rhat'))

    w <- which(colnames(rh) %in% vars)
    rh <- rh[,c(1,w)]

    high <- vector()
    for(i in 1:length(vars)){
        w <- which(rh[,i+1] > threshold)
        high <- append(high, w)
    }
    return(rh$date[high])
}


plot_kde_metab <- function(met, lim = NULL, col = "grey25"){

    kernel <- kde(na.omit(met[,c("GPP","ER")]))
    if(is.null(lim)) {
        lim <- quantile(c(met$GPP, -met$ER), .99, na.rm = T)
    }

    plot(kernel, xlab = "GPP (gO2m2d)", ylab = "ER (gO2m2d)",
         ylim = c(-lim, 0), xlim = c(0, lim),
         display = "filled.contour",
         cont=c(30,60,90), #lwd = 1,
         col = c('transparent',
                 alpha(col, .25),
                 alpha(col, .5),
                 alpha(col, .75)))

    abline(0,-1)
}

# K600 relationships ####
# plot K X ER relationship
plot_KxER <- function(fit, rm.bds = FALSE){

    if(inherits(fit, 'metab_bayes')){
        met <- fit@fit$daily

        p <- ggplot(met, aes(K600_daily_50pct, ER_daily_50pct)) +
            geom_point(size = 2) +
            xlab(expression(paste("K600 (d"^"-1"*")")))+
            ylab(expression(paste("Ecosystem Respiration (g"~O[2]~"m"^"-2"~"d"^"-1"*")")))+
            theme_bw()

        pcor <- cor(met$ER_daily_50pct, met$K600_daily_50pct,
                    method = 'pearson',
                    use = 'complete.obs')

        print(paste0('Pearson Correlation = ', pcor))

        return()
    }

    if(!inherits(fit, 'metab_bayes')){
        if('DO_fit' %in% colnames(fit)){

            p <- ggplot(fit, aes(K600, ER, col = DO_fit)) +
                geom_point(size = 2) +
                xlab(expression(paste("K600 (d"^"-1"*")")))+
                ylab(expression(paste("Ecosystem Respiration (g"~O[2]~"m"^"-2"~"d"^"-1"*")")))+
                theme_bw()
        } else {

            p <- ggplot(fit, aes(K600, ER)) +
                geom_point(size = 2) +
                xlab(expression(paste("K600 (d"^"-1"*")")))+
                ylab(expression(paste("Ecosystem Respiration (g"~O[2]~"m"^"-2"~"d"^"-1"*")")))+
                theme_bw()
        }

        if(rm.bds){
            fit$ER[fit$DO_fit == 'bad'] <- NA
            pcor <- cor(fit$ER, fit$K600,
                        method = 'pearson',
                        use = 'complete.obs')

            print(paste0('Pearson Correlation = ', pcor))

            return(p)
        }

        pcor <- cor(fit$ER, fit$K600,
                    method = 'pearson',
                    use = 'complete.obs')

        print(paste0('Pearson Correlation = ', pcor))

        return(p)
    }
}

# plot K x Q relationship
plot_KxQ <- function(fit, dat = NULL){

    if(inherits(fit, 'metab_bayes')){

        SM_day <- get_data_daily(fit) %>% select(date, discharge.daily)
        met <- fit@fit$daily %>%
            left_join(SM_day, by = 'date')

        p <- ggplot(met, aes(log(discharge.daily), K600_daily_50pct)) +
            geom_point(size = 2) +
            ylab(expression(paste("K600 (d"^"-1"*")")))+
            xlab(expression(paste("log Discharge (m"^"3"~"s"^"-2"*")")))+
            theme_bw()

        return(p)
    }

    if(!inherits(fit, 'metab_bayes')){
        if('solar.time' %in% colnames(dat)){
            daily <- dat %>%
                mutate(date = as.Date(solar.time)) %>%
                group_by(date) %>%
                summarize(across(-solar.time, mean, na.rm = T)) %>%
                ungroup()
        } else { daily <- dat}
        daily <- daily %>%
                left_join(fit, by = 'date')
        fit <- daily %>% select(date, discharge) %>% right_join(fit, by = 'date')

        if('DO_fit' %in% colnames(fit)){

            p <- ggplot(fit, aes(log(discharge), K600, col = DO_fit)) +
                geom_point(size = 2) +
                ylab(expression(paste("K600 (d"^"-1"*")")))+
                xlab(expression(paste("log Discharge (m"^"3"~"s"^"-2"*")")))+
                theme_bw()
        } else {

            p <- ggplot(fit, aes(log(discharge), K600)) +
                geom_point(size = 2) +
                ylab(expression(paste("K600 (d"^"-1"*")")))+
                xlab(expression(paste("log Discharge (m"^"3"~"s"^"-2"*")")))+
                theme_bw()
        }
        return(p)
    }
}

# plot KxQ relationship showing the kxq bins used when fitting and the posterior fit
plot_KxQ_bins <- function(fit, labs = TRUE, legend = TRUE){
    mm_fit <- get_fit(fit)

    SM_output <- mm_fit$daily
    SM_KQbin <-  mm_fit$KQ_binned
    SM_day <- get_data_daily(fit)
    SM_specs <- get_specs(fit)

    day <- data.frame(SM_day$discharge.daily,
                      SM_output$K600_daily_50pct,
                      SM_output$GPP_50pct,
                      SM_output$K600_daily_Rhat,
                      rep('daily', dim(SM_output)[1]))

    colnames(day)<-c('Q', 'K600', 'GPP','Rhat', 'Group')

    nodes<-data.frame(exp(SM_specs$K600_lnQ_nodes_centers),
                      exp(SM_KQbin$lnK600_lnQ_nodes_50pct),
                      exp(SM_KQbin$lnK600_lnQ_nodes_2.5pct),
                      exp(SM_KQbin$lnK600_lnQ_nodes_97.5pct),
                      exp(SM_specs$K600_lnQ_nodes_meanlog))
    prior_sd <- 1.96*exp(SM_specs$K600_lnQ_nodes_sdlog[1])
    colnames(nodes)<-c('Q', 'K600','K600_2.5', 'K600_97.5',  'K600_prior')
    xlim = c(min(c(log(day$Q), log(nodes$Q)), na.rm = T),
             max(c(log(day$Q), log(nodes$Q)), na.rm = T))
    ylim = c(min(c(day$K600, nodes$K600_2.5, nodes$K600_prior - prior_sd), na.rm = T),
             max(c(day$K600, nodes$K600_97.5, nodes$K600_prior + prior_sd), na.rm = T))

    if(labs){
        plot(log(day$Q), day$K600, type = 'n', xlim = xlim, ylim = ylim,
             xlab = expression(paste("log discharge (m"^"3"~"s"^"-1"*")")),
             ylab = expression(paste("K600 (d"^"-1"*")")))
    }
    if(!labs){
        par(mar = c(2, 2, 1, 1))
        plot(log(day$Q), day$K600, type = 'n', xlim = xlim, ylim = ylim,
             xlab = '', ylab = '')
        }
    if(legend){
        legend('topleft', legend = c("prior", "posterior", "data"),
               col = c( "brown3", "brown3", "grey25"), xpd = NA, inset = c(0, -0.1),
               pch = c(1, 19, 20), bty = 'n', ncol = 3)
    }
    polygon(log(c(nodes$Q, rev(nodes$Q))),
            c(nodes$K600_prior + prior_sd, rev(nodes$K600_prior - prior_sd)),
            col = alpha('brown3', 0.2), border=NA)
    points(log(day$Q), day$K600,
           col = "grey25", pch = 20, cex = 0.8)
    points(log(nodes$Q), nodes$K600_prior, col = "brown3", cex = 0.8)
    points(log(nodes$Q), nodes$K600, col = "brown3", pch = 19, cex = 0.8)
    arrows(log(nodes$Q), nodes$K600_2.5, log(nodes$Q), nodes$K600_97.5,
           length = 0, col = 'brown3', lwd = 2)

}



