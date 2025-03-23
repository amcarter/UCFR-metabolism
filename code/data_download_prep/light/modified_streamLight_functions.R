AppEEARS_proc2 <- function (unpacked_LAI, fit_method, plot = FALSE, write_output = FALSE,
          save_dir = NULL)
{
  site_proc <- function(Site, fit_method, plot, write_output,
                        save_dir) {
    SOI <- unpacked_LAI[[Site]]
    if (fit_method %in% c("AG", "Beck", "Elmore", "Gu",
                          "Klos", "Zhang")) {
      processed <- LAI_proc_phenofit2(SOI, fit_method = fit_method)
    }
    proc_na_rm <- na.omit(processed)
    r2 <- round(summary(lm(proc_na_rm[, "Lai"] ~ proc_na_rm[,
                                                            "LAI_proc"]))$adj.r.squared, 2)
    RMSE <- round(sqrt(mean((proc_na_rm[, "Lai"] - proc_na_rm[,
                                                              "LAI_proc"])^2, na.rm = TRUE)), 2)
    MAE <- round(mean(abs(proc_na_rm[, "Lai"] - proc_na_rm[,
                                                           "LAI_proc"]), na.rm = TRUE), 2)
    if (plot == TRUE) {
      ylabel <- expression(paste("LAI ", "(m"^{
        2
      }, "m"^{
        -2
      }, ")"))
      plot(processed[, "Date"], processed[, "Lai"], pch = 20,
           col = "grey60", ylim = c(0, 7), xlab = "Time (Days)",
           ylab = ylabel, main = paste("Site= ", Site))
      lines(processed[, "Date"], processed[, "LAI_proc"],
            col = "black", lwd = 2)
    }
    if (write_output == TRUE) {
      saveRDS(processed, paste0(save_dir, "/", Site, "_LAI_processed.rds"))
    }
    else {
      return(processed)
    }
  }
  if (write_output == TRUE) {
    lapply(names(unpacked_LAI), FUN = site_proc, fit_method = fit_method,
           plot = plot, write_output = write_output, save_dir = save_dir)
  }
  else {
    LAI_final <- lapply(names(unpacked_LAI), FUN = site_proc,
                        fit_method = fit_method, plot = plot, write_output = write_output,
                        save_dir = save_dir)
    names(LAI_final) <- names(unpacked_LAI)
    return(LAI_final)
  }
}


LAI_proc_phenofit2 <- function (Site, fit_method)
{
  ts <- Site
  ts$date <- as.Date(ts[, "pos_time"])
  ts$w <- QC_weights(SCF_QC = ts[, "FparLai_QC_SCF_QC"], wmin = 0.2,
                     wmid = 0.5, wmax = 1)
  date_diff <- diff(ts[, "date"], lag = 1)
  uniqv <- unique(date_diff)
  t_res <- uniqv[which.max(tabulate(match(date_diff, uniqv)))]
  ts_padded <- pad_ts(ts = ts, t_res = t_res)
  scene_count <- 365/t_res
  lai_input <- phenofit::check_input(t = ts_padded$date, y = ts_padded$Lai,
                                     w = ts_padded$w, nptperyear = scene_count, south = FALSE,
                                     maxgap = scene_count/4, alpha = 0.02, wmin = 0.2)
  breaks <- phenofit::season_mov(lai_input, FUN = smooth_wWHIT,
                                 wFUN = wTSM, r_min = 0.1, IsPlot = FALSE,
                                 IsPlot.OnlyBad = FALSE, print = FALSE, iters = 2, lambda = 1600)
  lai_fit <- phenofit::curvefits(lai_input, brks = breaks,
                                 list(methods = c(fit_method), wFUN = wTSM, nextend = 2, maxExtendMonth = 3,
                                 minExtendMonth = 1, minPercValid = 0.2, print = FALSE,
                                 verbose = FALSE, iters = 2, use.rough = FALSE, use.y0 = TRUE))
  l_param <- get_param(lai_fit)
  d_fit <- get_fitting(lai_fit)
  d_gof <- get_GOF(lai_fit)
  fit_bound <- setNames(data.frame(d_fit[d_fit$meth == fit_method,
  ]$t, d_fit[d_fit$meth == fit_method, ]$ziter2), c("date",
                                                    "Lai_fit"))
  full_dates <- setNames(data.frame(seq.Date(from = min(ts[,
                                                           "date"]), to = max(ts[, "date"]), by = "day")), "date")
  full_dates$Year <- as.numeric(format(full_dates[, "date"],
                                       "%Y"))
  full_dates$DOY <- as.numeric(format(full_dates[, "date"],
                                      "%j"))
  interpolated <- approx(fit_bound[, "date"], fit_bound[,
                                                        "Lai_fit"], xout = full_dates[, "date"])
  interpolated_df <- setNames(data.frame(interpolated$x, interpolated$y),
                              c("date", "LAI_proc"))
  dfs <- list(full_dates, ts[, c("date", "Lai")], interpolated_df)
  final_merged <- plyr::join_all(dfs, by = "date")
  colnames(final_merged)[1] <- "Date"
  return(final_merged)
}


AppEEARS_unpack_QC_2 <- function (zip_file, zip_dir, request_sites)
{
    request <- sub("\\..*", "", zip_file)
    grx_exp <- glob2rx(paste("*-results.csv", sep = ""))
    folder_files <- list.files(paste(zip_dir, "/", request,
                                     sep = ""))
    results_fn <- folder_files[grep(grx_exp, folder_files)][1]
    res_split <- strsplit(results_fn, "-")[[1]]
    prod_loc <- grep(glob2rx(paste("MCD*", sep = "")), res_split)
    raw <- data.frame(data.table::fread(paste0(zip_dir, "/",
                                               request, "/", results_fn)))
    heading <- paste(strsplit(results_fn, "-")[[1]][prod_loc:(prod_loc +
                                                                  1)], collapse = "_")
    QC_cols <- paste0(heading, c("_FparLai_QC_MODLAND", "_FparLai_QC_Sensor",
                                 "_FparLai_QC_DeadDetector", "_FparLai_QC_CloudState",
                                 "_FparLai_QC_SCF_QC"))
    reduced <- data.frame(raw[, c("ID", "Date", paste0(heading,
                                                       c("_Lai_500m", "_LaiStdDev_500m")))], apply(raw[, QC_cols],
                                                                                                   MARGIN = 2, FUN = function(x) {
                                                                                                       gsub(".*b", "", x)
                                                                                                   }))
    colnames(reduced)[3:9] <- c("Lai", "Lai_sd", "FparLai_QC_MODLAND",
                                "FparLai_QC_Sensor", "FparLai_QC_DeadDetector", "FparLai_QC_CloudState",
                                "FparLai_QC_SCF_QC")
    reduced$pos_time <- as.POSIXct(reduced[, "Date"], format = "%Y-%m-%d",
                                   tz = "UTC")
    reduced$Year <- as.numeric(strftime(reduced[, "pos_time"],
                                        format = "%Y", tz = "UTC"))
    reduced$DOY <- as.numeric(strftime(reduced[, "pos_time"],
                                       format = "%j", tz = "UTC"))
    VOI <- reduced[, c("ID", "pos_time", "Year", "DOY", "Lai",
                       "Lai_sd", "FparLai_QC_MODLAND", "FparLai_QC_Sensor",
                       "FparLai_QC_DeadDetector", "FparLai_QC_CloudState",
                       "FparLai_QC_SCF_QC")]
    site_split <- split(VOI, VOI[, "ID"])
    remove_null <- function(Site) {
        if (all(Site[, "FparLai_QC_SCF_QC"] == "100") == FALSE) {
            return(Site)
        }
    }
    split_filtered <- lapply(site_split, FUN = remove_null)
    split_na_rm <- split_filtered[!sapply(split_filtered, is.null)]
    for (i in 1:length(split_na_rm)) {
        site_name <- request_sites[gsub("[[:punct:]]", "", request_sites) %in%
                                       unique(split_na_rm[[i]][, "ID"])]
        split_na_rm[[i]][, "ID"] <- site_name
    }
    ID_DF <- setNames(data.frame(gsub("[[:punct:]]", "", request_sites),
                                 request_sites), c("no_punct", "Site_ID"))
    ID_merge <- merge(setNames(data.frame(names(split_na_rm)),
                               "no_punct"), ID_DF, by = "no_punct")
    rownames(ID_merge) <- ID_merge[, "no_punct"]
    names(split_na_rm) <- ID_merge[names(split_na_rm), "Site_ID"]
    missing <- request_sites[!(request_sites %in% names(split_na_rm))]
    if (length(missing) != 0) {
        message(paste("The following sites did not have LAI data in this request:",
                      paste(missing, sep = "", collapse = ", ")))
    }
    return(split_na_rm)
}
