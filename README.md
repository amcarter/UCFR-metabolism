# UCFR-metabolism
Data and code for estimating DO metabolism and building biomass models in the Upper Clark Fork River for the CREWS project

**Contents**
  
1. [Data Sets](#data-sets-description)  
    - [Data files](#data-files)  
    - [Working Datasets](#working-datasets)  
2. [Workflow](#workflow)  
    - [Data Download Prep](#data-download-prep)  
    - [Metabolism](#metabolism)  
    - [Biomass](#biomass)  
    - [Model](#model)  
    - [Stan Code](#stan-code)  
    - [NPP Calculations](#npp-calculations)  


<!-- Data Sets description -->
## Data Sets Description

1. [DO data product](https://portal.edirepository.org/nis/metadataviewer?packageid=edi.1324.1) - Raw and cleaned dissolved oxygen data from 2020 and 2021 in the Upper Clark Fork River.

2. [Algal biomass data product](url) - 

<!-- Data files -->
### Data files
#### site_data

**a. site_data.csv**  -  Metadata for all Upper Clark Fork River sites.

**b. UCFR_depth_summary.csv** - average depth estimates taken at different discharge levels at sites along the Upper Clark Fork River.

**c. sw_radiation_all_sites.csv**, **hourly_modeled_light_all_sites.csv**, **daily_modeled_light_all_sites.csv** - light data downloaded from [NLDAS](https://appeears.earthdatacloud.nasa.gov/) using the [streamlight R package](https://github.com/psavoy/StreamLight) for all sites

**d. air_pressure_temp.csv** - air pressure and air temperature data downloaded from the DeerLodge airport. Released as part of DO data product (see above).

**e. discharge_UCFR_sites_2020** - discharge data for all sites downloaded from USGS gage stations.

#### raw_DO 
Local copies of raw dissolved oxygen data for running code. See above data product for official version.

#### prepared_data 
DO sensor data that has been cleaned, DO sat calculated based on water temperature and air pressure, paired with depth estimates based on discharge, paired with modeled light, and converted to standardized time units of solar time using the StreamMetabolizer R package. Data in the main folder have days with poor DO fits with StreamMetabolizer removed after an initial run.

**a. initial_run_prepped_data** - Prepared data without days removed. Used for initial Stream Metabolizer run.

**b. compiled_prepared_data.csv** - complete prepared datasets compiled into a single file.

#### biomass_data

**a. INTEGRATED_BIOMASS_FOR_AC_RFL.csv**  -  Raw biomass data with AFDM estimates and chlorophyll a split into different biomass categories.

**b. BIOMASS_DatasetAttributes.csv**  -  Reference list for raw biomass data file variable names.  

**c. biomass_working_data.csv** - a working datafile created by cleaning up the raw data and selecting the components relevant to this project.

**d. biomass_working_data_summary.csv** - a summary of the working dataset.

**e. log_gam_fits_biomass.csv** - estimates from HGAMs fit on log biomass data.

**f. log_gam_smoothness_parameter_checks.csv** - goodness of fit checks on HGAM fits on log data.

#### metabolism
Metabolism estimates from Stream Metabolizer

**a. metabolism_compiled_all_sites_2000iter_bdr_kss01.csv** - this file and those named similarly are compiled metabolism estimates from bayesian metabolism runs in stream metabolizer. The filename indicates the number of iterations of the chains that were run (2000 in this case), bdr indicates that bad DO fit days were removed before running, and kss shows the value of the prior on K600 sigma sigma, in this case, 0.01.

**b. metab_for_results.csv** - a cleaned up version of the finalized metabolism estimates that is used to generate the results in the latex manuscript and the data table in the SI.

**c. days_with_poor_DO_fits_initial_run.csv** - days where the modeled DO was a poor fit to the data in the initial runs.

**d. days_with_poor_DO_fits.csv** - days where the modeled DO was a poor fit to the data and days when the Rhat metric was greater than 1.1 in the final metabolism runs.

#### model_fits

**a. biomass_metab_model_data.csv** - compiled biomass, metabolism, light, and discharge for running biomass-metabolism models.

**b. brms_gpp_models.rds** - trial bayesian linear regression models.

**c. hierarcical_model_parameters.csv** - parameter estimates for hierarchical AR1 model runs

**d. hierarcical_model_parameters_logGPP.csv** - parameter estimates for hierarchical AR1 model runs on log GPP data.

**e. hierarcical_model_posterior_distribution_for_plots.csv** - full posterior distributions for hierarchical AR1 model runs.
    
    
<!-- Workflow -->
## Workflow

<!-- data-download-prep -->
### Data Download and Preparation

**1. download_air_temp_pres.R ** - download raw air temperature and air pressure data from the Deer Lodge Airport.

**2. light (folder)** Scripts and data for getting light from StreamLight

**2a. get_light_from_StreamLight_UCFR.R** - workflow to download light data from NLDAS and imagery from MODIS to calculate light at the stream surface based on the package [StreamLight](github.com/psavoy/StreamLight.com)

**2b. modified_streamLight_functions.R** - modified functions from the package StreamLight 

**3. clean_raw_DO.R ** - Clean raw dissolved oxygen datasets. See data product above for details.

**4. prep_metabolism_data.R** - reformats metabolism data for running in StreamMetabolizer by calculating DO sat, adding depth based on discharge estimates, modeling light and converting time to units of solar.time based on latitude and longitude.

<!-- metabolism --> 
### Metabolism

**1. UCFR depth model.R** - uses BRMS to generate a model of average depth as a function of discharge with estimates pooled across all six sites.

**2. UCFR_run_streamMetabolizer_batched_years.R** - run StreamMetabolizer on all six sites using data from the prepared data folder. script is set up to run with a bayesian model with normal pooling of K600.

**3. UCFR_run_streamMetabolizer_fixed_K.R** - run StreamMetabolizer on all six sites with K600 fixed at the median estimates generated using a normally pooled bayesian model. Script is set up to run a mle model.

**4. identify_poor_DO_fits.R** - visually identify days where the modeled DO is a poor fit to the data in the initial streamMetabolizer model runs. 

**5. compile_metabolism.R** - compile estimates of metabolism from a batch of runs. Code extracts estimates of GPP, ER, K600 and flags any days with visually identified poor DO fits or days with Rhat above 1.1.

**6. functions_examine_SM_output.R** - functions to interact with metabolism fit objects produced by streamMetabolizer.

**7. compare_across_k600_sigmasigmas.R** - Compile metabolism estimates run with different priors on K600 sigma sigma and compare the difference in metabolism estimates.

<!--biomass-->
### Biomass

**1. create_working_biomass_dataset.R** - clean up raw biomass datafile, select and summarize data relevant to this project. 

**2. compare_biomass_and_metab.R**- an exploratory comparison of raw biomass data with metabolism estimates.

**3. fit_HGAM_to_biomass.R** - fit hierarchical GAMs to biomass datasets using mcgv package. Compares across different model fits. Also generates plots of fitted biomass data.

<!--model-->
### Model

**1. prepare_biomass_metab_model_data.R** - compiles metabolism, biomass, light, discharge into a single file for running biomass metabolism models. 

**2. fit_linear_brms_models.R** - trial model fits using bayesian linear regression in the brms package.

**3. compare_ar1_linearmod.R** - compare coefficients generated using a basic linear model to those generated from an ar1 model.

**4. simulate_ar1_data.R** - simulate AR1 data for testing parameter recovery in stan models.

**5. run_hierarchical_ar1_models.R** 

<!--stan-code-->
### Stan Code

<!--npp-calculations-->
### NPP Calculations

**1. fit_quantile_regression.R** - fits a multilevel quantile regression model to the GPP and ER estimates from each site year to estimate the fraction of GPP that is respired as autotrophic respiration. Based on this calculation, it estimates NPP for each day. Part two of this code builds a second multilevel linear regression to estimate daily NPP as a function of chlorophyll in the epilithon and filamentous fraction with an intercept of zero. Based on these rates, we calculate turnover times of each biomass fraction. This script also makes several plots of these data. 
