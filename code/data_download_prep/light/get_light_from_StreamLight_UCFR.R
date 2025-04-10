# Get light data using the StreamLight package from P Savoy

# setup ####
#Use the devtools packge to install StreamLightUtils
# note this package is huge. you might have to adjust the timeout settings to install it
# devtools::install_github("psavoy/StreamLightUtils")
# devtools::install_github("psavoy/StreamLight")

library(StreamLightUtils)
library(StreamLight)
library(lubridate)
library(tidyverse)
setwd('C:/Users/alice.carter/git/UCFR-metabolism/')

source('code/data_download_prep/light/modified_streamLight_functions.R')
sitedat <- read_csv('data/site_data/site_data.csv') %>%
    dplyr::filter(!is.na(sitecode))%>%
    dplyr::rename(Lat = Latitude, Lon = Longitude) %>%
    mutate(startDate = as.Date('2020-06-01'),
           endDate = as.Date('2021-11-30'),
           Site_ID = paste0('nwis_', nwis_code))
# Download and Process NLDAS data for incoming radiation ####
# Set the download location (add your own directory)
base_dir <- "C:/Users/alice.carter/git/UCFR-metabolism"
working_dir <- paste0(base_dir, '/code/data_download_prep/light/NLDAS')

# for multiple sites
#Read in a table with initial site information
sites <- sitedat %>%
  select(Site_ID = sitecode, Lat, Lon, startDate) %>%
  mutate(startDate = as.character(date(startDate)),
         epsg_crs = rep(4269, nrow(sitedat))) %>%
  data.frame()

#Download NLDAS data at UCFR
NLDAS_DL_bulk(
  save_dir = working_dir,
  site_locs = sites
)

#List of successfully downloaded sites
NLDAS_list <- stringr::str_sub(list.files(working_dir), 1, -11)

#Processing the downloaded NLDAS data
NLDAS_processed <- StreamLightUtils::NLDAS_proc(read_dir = working_dir, NLDAS_list)

for(i in 1:length(NLDAS_processed)){
    NLDAS_processed[[i]]$sitecode <- names(NLDAS_processed)[i]
}

NLDAS_data <- Reduce(bind_rows, NLDAS_processed) %>%
    as_tibble() %>%
    mutate(date = as.Date(paste0(Year, '-', DOY), format = '%Y-%j')) %>%
    filter(date <= as.Date('2021-11-30'))

ggplot(NLDAS_data, aes(DOY, SW, col = sitecode)) +
    geom_line() +
    facet_wrap(.~Year, scales = 'free')

setwd('C:/Users/alice.carter/git/UCFR-metabolism/')
write_csv(NLDAS_data, 'data/site_data/sw_radiation_all_sites.csv')

# Download and process MODIS LAI ####
#Make a table for the MODIS request
request_sites <- sites[, c("Site_ID", "Lat", "Lon")]

#Export your sites as a .csv for the AppEEARS request
write.table(
  request_sites,
  "code/data_download_prep/light/MODIS/UCFR_sites.csv",
  sep = ",",
  row.names = FALSE,
  quote = FALSE,
  col.names = FALSE
)

# Save the zip file downloaded from AppEEARS
# Unpack the data after downloading
working_dir <- "C:/Users/alice.carter/git/UCFR-metabolism/code/data_download_prep/light/MODIS"
setwd(working_dir)

# do this manually because this function doesn't seem to work
MOD_unpack <- AppEEARS_unpack_QC_2(
  zip_file = "ucfr-sites.zip",
  zip_dir = working_dir,
  request_sites[, "Site_ID"]
)


# Process the downloaded data
png('LAI_focal_sites.png')
par(mfrow = c(2,3))
MOD_processed <- AppEEARS_proc2(
  unpacked_LAI = MOD_unpack,
  fit_method = "Gu",
  plot = TRUE
)
dev.off()

# Make an input datafile for streamlight ####
working_dir <- paste0(base_dir, "/code/data_download_prep/light/drivers")
make_driver(sites, NLDAS_processed, MOD_processed,
            write_output = TRUE, working_dir)

# add parameters to site file ####
site_parm <- sites %>%
  select(Site_ID , Lat, Lon)

# Calculate site azimuth using Google Maps. The angle I am using is the
# line that connects the sensor location with the point 1000 m upstream.
# angle is in degrees from north
n = nrow(site_parm)
site_parm <- site_parm %>%
  mutate(Width = rep(20, n),
         epsg_code = rep(4269, n),
         datum = rep('NAD83', n),
         Azimuth = c(0, 170, 135, 105, 95, 107),
         WL = rep(1, n), # water level
         BH = rep(0, n), # bank height
         BS = rep(100, n), # bankslope
         TH = rep(0, n), # tree height (m)
         x = rep(0, n),
         overhang = rep(0, n),
         overhang_height = rep(0, n)) %>%
  data.frame()

# get canopy data:
# extract_height(Site_ID = sites[,'Site_ID'], Lat = sites[,'Lat'],
#                Lon = sites[,'Lon'], site_crs = sites[,'epsg_code'])
# that didn't work
# from intuition, trees ~ 10 m tall, covering only 10-20% of stream length = 2m


# running streamlight ####
#Function for batching over multiple sites
batch_model <- function(Site, params, read_dir, save_dir){
  ss <- list.files(read_dir)
  if(!paste(Site, "_driver.rds", sep = "") %in% ss){
    return(NA)
  }
  #Get the model driver
  driver_file <- readRDS(paste(read_dir, "/", Site, "_driver.rds", sep = ""))

  #Get model parameters for the site
  site_p <- params[params[, "Site_ID"] == Site, ]

  #Run the model
  modeled <- stream_light(
    driver_file,
    Lat = site_p[, "Lat"],
    Lon = site_p[, "Lon"],
    channel_azimuth = site_p[, "Azimuth"],
    bottom_width = site_p[, "Width"],
    BH = site_p[, "BH"],
    BS = site_p[, "BS"],
    WL = site_p[, "WL"],
    TH = site_p[, "TH"],
    overhang = site_p[, "overhang"],
    overhang_height = site_p[, "overhang_height"],
    x_LAD = site_p[, "x"]
  )

  #Save the output
  saveRDS(modeled, paste(save_dir, "/", Site, "_predicted.rds", sep = ""))

} #End batch_model

#Applying the model to all sites
model_rd <- working_dir
model_sd <- working_dir

#Running the model
lapply(
  site_parm[, "Site_ID"],
  FUN = batch_model,
  params = site_parm,
  read_dir = model_rd,
  save_dir = model_sd
)

#Take a look at the output
predicted <- readRDS(paste(working_dir, '/BN_predicted.rds', sep = ''))

predicted %>%
  mutate(date = as.Date(local_time)) %>%
  group_by(date) %>%
  summarize(par_surface = sum(PAR_surface),
            par_inc = sum(PAR_inc)) %>%
  ungroup() %>%
  ggplot(aes(date, par_inc)) +
    geom_point(col = 'grey') +
    geom_point(aes(y = par_surface), col = 'gold') +
    theme_minimal()

# Combine the outputs into one file
dat <- data.frame()
dat_daily <- data.frame()
for(site in as.character(site_parm[,'Site_ID'])){
  fl <- list.files(working_dir)
  if(!paste(site, "_driver.rds", sep = "") %in% fl) next

  pred <- readRDS(paste(working_dir, '/', site, '_predicted.rds', sep = ''))
  ss <- pred %>%
  dplyr::select(local_time, LAI, PAR_inc, PAR_surface) %>%
  mutate(date = as.Date(local_time, tz = 'America/Denver'),
         site = site)
  dat_daily <- bind_rows(dat_daily, ss)

  ss <- ss %>%
  group_by(site, date) %>%
  summarize(light_hrs = length(which(PAR_inc>0)),
            LAI = mean(LAI, na.rm = T),
            PAR_inc = sum(PAR_inc),
            PAR_surface = sum(PAR_surface)) %>%
  ungroup()

  dat <- bind_rows(dat, ss)
}
dat %>%
ggplot(aes(date, PAR_surface, col = site)) +
  geom_point()

write_csv(dat, 'C:/Users/alice.carter/git/UCFR-metabolism/data/site_data/daily_modeled_light_all_sites.csv')
write_csv(dat_daily, 'C:/Users/alice.carter/git/UCFR-metabolism/data/site_data/hourly_modeled_light_all_sites.csv')
dat <- read_csv('C:/Users/alice.carter/git/UCFR-metabolism/data/site_data/daily_modeled_light_all_sites.csv')

