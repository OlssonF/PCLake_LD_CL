#--------------------------------------#
## Project: Lake District pathways modelling
## Script purpose: Run physical lake model (FLake)
## Date: 2026-02-24
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

# For each of the individual lakes that we have appropriate data for we need to run the 
# physical lake model to derive an annual cycle of water temperatures that can be given to
# PCLake+.

# Model code: http://flake.igb-berlin.de/site/download (accessed 23rd Feb 2026)
# FLake is run twice with two light extinction values ("clear" and "turbid")

# The model is run for 10 years to obtain multiple years from which a more
# average condition can be obtained

# Uses the forcing data from E-OBS,v32.0 (or v30.0 for qq and fg), 0.1 degree grid
# e.g. for air temperature the data are downloaded in two chunks from here
# https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
# "We acknowledge the E-OBS dataset from the EU-FP6 project UERRA (https://www.uerra.eu) and 
# the Copernicus Climate Change Service, and the data providers in the ECA&D project (https://www.ecad.eu)"
# "Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: An Ensemble Version of 
# the E-OBS Temperature and Precipitation Datasets, J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200"

# Or use forcing data from ERA5-Land reanalysis dataset, 0.1 degree grid
# data obtained using the download_era5.py script from https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land

# Attribution:
# Generated using or contains modified Copernicus Climate Change Service information <2019>. 
# Neither the European Commission nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains.
# Muñoz Sabater, J. (2019): ERA5-Land hourly data from 1950 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). 
# DOI: 10.24381/cds.e2161bac (Accessed on 20-Mar-2026)

library(tidyverse)
library(mgcv)
source('R/flake_functions.R')

subset <- F # run all lakes or not
# the lakes portal data has all the basic info we need
lakes_portal_df <- read_csv('data/lakes4PCLake.csv', show_col_types = F)

if (subset == T) {
  # LakeIDs to loop through
  lakeIDs <- readxl::read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') |> 
    filter(!is.na(LAKE_LakesTour2021_Zooplankton.csv),
           !str_detect(LAKE_LakesTour2021_Zooplankton.csv, 'Loughrigg'))   # remove Loughrigg because it's not right
  
  lake_names_lookup <- lakeIDs$`WBID_Lake District_UKCEH Portal data_raw.xlsx`
  names(lake_names_lookup) <- str_extract(lakeIDs$LAKE_LakesTour2021_Zooplankton.csv, "....")
  
  lake_names_lookup[which(names(lake_names_lookup) == 'Wind')] <- 29233 
  # the NBAS and SBAS have seperate WBIDs but the one from the lakes portal has the combined one which is different
} else {
  
  lake_names_lookup <- lakes_portal_df$WBID
  names(lake_names_lookup) <- lakes_portal_df$NAME
  
}


# Run with "clear" light extinction -------------
for (i in 1:length(lake_names_lookup)) {
  # Select one lake at a time
  wbid <- lake_names_lookup[i]
  message('Running ', wbid, ' ', names(wbid))
  
  # Obtain the lake portal data (fetch, depth etc.)
  lakes_portal_subset <- lakes_portal_df |>
    filter(WBID == wbid)
  
  latitude <- lakes_portal_subset |>
    select(WBLAT) |> pull()
  
  longitude <- lakes_portal_subset |>
    select(WBLONG) |> pull()
  
  elev <- lakes_portal_subset |>
    select(WBALT) |> pull()
  
  # lake_fetch <- lakes_portal_subset$FETCH_KM * 1000 # fetch, convert from km to m - is this actually a good estimate (maximum fetch)
  lake_fetch <- sqrt(lakes_portal_subset$WBSAREA*10000) # convert from Ha->m2 (typical fetch)
  lake_depth <- lakes_portal_subset$MNDP # mean depth
  
  # Set up the drivers and nml files
  setup_FLake(wbid,
              latitude, longitude,
              elev, lake_depth, lake_fetch,
              lake_lightext = 'clear',
              calc_cc = F, make_met = T, use_annual = F,
              outputfile = file.path(wbid, paste0(wbid,'_clear.rslt')))
  
  # run FLake
  run_FLake(wbid)
  
}



# Run with "turbid" light extinction -------------
for (i in 1:length(lake_names_lookup)) {
  # Select one lake at a time
  wbid <- lake_names_lookup[i]
  message('Running ', wbid, ' ', names(wbid))
  
  # Obtain the lake portal data (fetch, depth etc.)
  lakes_portal_subset <- lakes_portal_df |>
    filter(WBID == wbid)
  
  latitude <- lakes_portal_subset |> 
    select(WBLAT) |> pull()
  
  longitude <- lakes_portal_subset |> 
    select(WBLONG) |> pull()
  
  elev <- lakes_portal_subset |> 
    select(WBALT) |> pull()
  
  # lake_fetch <- lakes_portal_subset$FETCH_KM * 1000 # fetch, convert from km to m - is this actually a good estimate (maximum fetch)
  lake_fetch <- sqrt(lakes_portal_subset$WBSAREA*10000) # convert from Ha->m2 (typical fetch)
  lake_depth <- lakes_portal_subset$MNDP # mean depth
  
  # Set up the drivers and nml files
  setup_FLake(wbid,
              latitude, longitude,
              elev, lake_depth, lake_fetch,
              lake_lightext = 'turbid',
              calc_cc = F, make_met = F, use_annual = F,
              outputfile = file.path(wbid, paste0(wbid,'_turbid.rslt')))
  
  # run FLake
  run_FLake(wbid)
  
}


#----------------------------------------------------------#

# Calculate a smoothed run -------------------------------------
flake_dir <- 'data/flake'

flake_results <- list.files(flake_dir, pattern = '*.rslt', recursive = T, full.names = T)
flake_nmls <- list.files(flake_dir, pattern = '*.nml', recursive = T, full.names = T)

WIND_WBID <- 29233

# Temps - fit GAM model -------------------------------------------

# Get some rough predictions of water temperature dynamics that can be used in PCLake

for (i in 1:length(lake_names_lookup)) {
  
  lake_ID_use <- lake_names_lookup[i]
  flake_IDs <- flake_results[str_detect(flake_results, as.character(lake_ID_use))]
  flake_nml <-  glmtools::read_nml(flake_nmls[str_detect(flake_nmls, as.character(lake_ID_use))])
  
  lake_name_use <- names(lake_names_lookup)[i]
  
  ## read the FLake results ------------------
  clear_df <- read_flake(flake_IDs[1]) |> 
    mutate(doy = yday(as_date(time)),
           year = year(as_date(time)))|> 
    # slice_tail(n = 365) |> 
    mutate(time = row_number())
  
  turbid_df <- read_flake(flake_IDs[2]) |> 
    mutate(doy = yday(as_date(time)),
           year = year(as_date(time)))|> 
    # slice_tail(n = 365) |> 
    mutate(time = row_number())
  
  #----------------------------------------#
  
  ## Calculate the mean of the two light extinction runs
  mean_df <- as.data.frame(Map(function(x, y) {(x + y) / 2}, clear_df, turbid_df))
  
  #-------------------------------------------------------#
  # ----------extract temps for PCLake --------------------
  #-------------------------------------------------------#
  
  # extract a typical seasonal pattern by fitting a GAM model 
  # Uses a cyclic basis function, smoother with DOY
  # fits the data on the mean_df (mean of the clear and turbid predictions)
  
  # Fit GAM with cyclic cubic spline
  gam_model_Ts <- gam(Ts ~ s(doy, bs = "cc"), data = mean_df)
  gam_model_Tb <- gam(Tb ~ s(doy, bs = "cc"), data = mean_df)
  
  # Create a sequence of DOY values to predict for
  newdata <- data.frame(doy = 1:365)
  
  # Extract fitted values for each DOY
  fitted_vals_Ts <- predict(gam_model_Ts, newdata = newdata)
  fitted_vals_Tb <- predict(gam_model_Tb, newdata = newdata)
  
  # Combine DOY and fitted values
  result <- data.frame(doy = newdata$doy,
                       Ts = round(fitted_vals_Ts, digits = 2),
                       Tb = round(fitted_vals_Tb, digits = 2))
  
  p_results <- ggplot(result, aes(x=doy)) + 
    geom_line(aes(y = Ts, colour = 'surface'))  + 
    geom_line(aes(y = Tb, colour = 'bottom')) +
    coord_cartesian(ylim = c(0,30)) +
    theme_bw() +
    labs(title = "Fitted GAMs", subtitle = paste(lake_ID_use, lake_name_use))
  
  ggsave(p_results ,filename = file.path(flake_dir, 'plots', paste0(lake_ID_use, '_predictions.png')),
         width = 15, height = 10, units = 'cm')
  
  # Write to file for PCLake
  write_csv(result, file = file.path(flake_dir, lake_ID_use, paste0(lake_ID_use,'_predictions.csv')))
  
}



# Calculate a timeseries of mixed depths -------------------------------------

# ML - fit GAM model -------------------------------------------

# Get some rough predictions of mixed depth dynamics that can be used in PCLake

for (i in 1:length(lake_names_lookup)) {
  
  lake_ID_use <- lake_names_lookup[i]
  flake_IDs <- flake_results[str_detect(flake_results, as.character(lake_ID_use))]
  flake_nml <-  glmtools::read_nml(flake_nmls[str_detect(flake_nmls, as.character(lake_ID_use))])
  
  lake_name_use <- names(lake_names_lookup)[i]
  
  mean_depth <- lakes_portal_df |>
    filter(WBID == lake_ID_use) |> 
    pull(MNDP)
  
  ## read the FLake results ------------------
  clear_df <- read_flake(flake_IDs[1]) |> 
    mutate(doy = yday(as_date(time)),
           year = year(as_date(time)))|> 
    # slice_tail(n = 365) |> 
    mutate(time = row_number())
  
  turbid_df <- read_flake(flake_IDs[2]) |> 
    mutate(doy = yday(as_date(time)),
           year = year(as_date(time)))|> 
    # slice_tail(n = 365) |> 
    mutate(time = row_number())
  
  #----------------------------------------#
  
  ## Calculate the mean of the two light extinction runs
  mean_df <- as.data.frame(Map(function(x, y) {(x + y) / 2}, clear_df, turbid_df))
  
  mean_df_medi <- mean_df |> 
    reframe(.by = doy, 
            h_ML = median(h_ML)) 
  
  mean_df |> 
    reframe(.by = doy, 
            h_ML = mean(h_ML)) |> 
    mutate(strat = ifelse(h_ML == mean_depth | h_ML == 0, F, T)) |> 
    ggplot(aes(x=doy, y = h_ML)) +
    geom_line(aes()) +
    geom_point(aes(colour = strat)) +
    geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cc",  k= 10)) +
    geom_hline(yintercept = mean_depth) +
    coord_cartesian(ylim = c(0, mean_depth), reverse = 'y')
  
  #-------------------------------------------------------#
  # ----------extract temps for PCLake --------------------
  #-------------------------------------------------------#
  
  # extract a typical seasonal pattern by fitting a GAM model 
  # Uses a cyclic basis function, smoother with DOY
  # fits the data on the mean_df_medi (median year of the mean of the clear and turbid predictions)
  
  # Fit GAM with cyclic cubic spline
  gam_model_ML <- gam(h_ML ~ s(doy, bs = "cc"), data = mean_df_medi)
  
  # Create a sequence of DOY values to predict for
  newdata <- data.frame(doy = 1:365)
  
  # Extract fitted values for each DOY
  fitted_vals_ML <- predict(gam_model_ML, newdata = newdata)
  
  # Combine DOY and fitted values
  result <- data.frame(doy = newdata$doy,
                       h_ML = round(fitted_vals_ML, digits = 2))
  
  # find the dates of the primary stratified period
  summer_strat <- get_summerstrat(result$h_ML, mean_depth)
  
  # omit stratification outside of these periods
  result <- result |> 
    mutate(h_ML = ifelse(between(doy, 
                                 summer_strat[1],
                                 summer_strat[2]),
                         h_ML, mean_depth)) 
  
  p_results <- ggplot(result, aes(x=doy)) + 
    geom_line(aes(y = h_ML)) +
    coord_cartesian(ylim = c(0,mean_depth), reverse = 'y') +
    geom_hline(yintercept = mean_depth, linetype = 'dashed') +
    theme_bw() +
    labs(title = "Fitted GAMs", subtitle = paste(lake_ID_use, lake_name_use))
  
  ggsave(p_results ,filename = file.path(flake_dir, 'plots', paste0(lake_ID_use, '_MLpredictions.png')),
         width = 15, height = 10, units = 'cm')
  

  # Write to files for PCLake
  # the mixed layer depth
  write_csv(result, file = file.path(flake_dir, lake_ID_use, paste0(lake_ID_use,'_MLpredictions.csv')))
  
  # the stratification period
  result |> 
    mutate(strat = ifelse(h_ML == mean_depth, 0, 1)) |> 
    select(doy, strat) |> 
    write_csv(file = file.path(flake_dir, lake_ID_use, paste0(lake_ID_use,'_stratpredictions.csv')))
  
}




