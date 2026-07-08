#--------------------------------------#
## Project: Lake District pathways modelling
## Script purpose: Run physical lake model (FLake)
## Date: 2026-02-24
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

# For each of the individual lakes that we have appropriate data for we need to run the 
# physical lake model to derive an average annual cycle of water temperatures that can be given to
# PCLake+.

# Model code: http://flake.igb-berlin.de/site/download (accessed 23rd Feb 2026)
# We will use FLake to derive the "perpetual year solution" based on a single year of driving data
# FLake is run twice with two light extinction values ("clear" and "turbid")

# A perpetual year solution represents the annual cycle of temperature and mixing in a given lake 
# that corresponds to a given annual cycle of input meteorological quantities. 
# Starting with arbitrary initial conditions, the year-long simulation is repeated, using one and 
# the same annual cycle of forcing. The initial conditions for the next year-long run are specified
# using the values of the lake-model at the end of the previous year-long run. 
# After a few model years, a periodic "perpetual year" solution is obtained. 

# A perpetual year solution obtained with FLake is useful to give an idea of the state of the lake if 
# no data from in the lake are available (but the atmospheric forcing can be specified in a rational way).

# Uses the forcing data from E-OBS,v32.0 (or v30.0 for qq and fg), 0.1 degree grid
# e.g. for air temperature the data are downloaded in two chunks from here
# https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php

# "We acknowledge the E-OBS dataset from the EU-FP6 project UERRA (https://www.uerra.eu) and 
# the Copernicus Climate Change Service, and the data providers in the ECA&D project (https://www.ecad.eu)"
# "Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: An Ensemble Version of 
# the E-OBS Temperature and Precipitation Datasets, J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200"
library(tidyverse)
source('R/flake_functions.R')

# the lakes portal data has all the basic info we need
lakes_portal_df <- read_csv('data/lakes4PCLake.csv')

# LakeIDs to loop through
lakeIDs <- readxl::read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') |> 
  filter(!is.na(LAKE_LakesTour2021_Zooplankton.csv),
         !str_detect(LAKE_LakesTour2021_Zooplankton.csv, 'Loughrigg'))   # remove Loughrigg because it's not right

lake_names_lookup <- str_extract(lakeIDs$LAKE_LakesTour2021_Zooplankton.csv, "....")
names(lake_names_lookup) <- lakeIDs$`WBID_Lake District_UKCEH Portal data_raw.xlsx`

names(lake_names_lookup)[which(lake_names_lookup == 'Wind')] <- 29233 
# the NBAS and SBAS have seperate WBIDs but the one from the lakes portal has the combined one which is different

# Run with "clear" light extinction -------------
for (i in 1:length(lake_names_lookup)) {
  # Select one lake at a time
  lake_name <- lake_names_lookup[i]
  message('Running ', lake_name)

  # Obtain the lake portal data (fetch, depth etc.)
  lakes_portal_subset <- lakes_portal_df |>
    filter(str_detect(NAME, lake_name))

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
  setup_FLake(lake_name,
              latitude, longitude,
              elev, lake_depth, lake_fetch,
              lake_lightext = 'clear',
              calc_cc = F, make_met = T, use_annual = F,
              outputfile = file.path(lake_name, paste0(lake_name,'_clear_obs.rslt')))

  # run FLake
  run_FLake(lake_name)

}



# Run with "turbid" light extinction -------------
for (i in 1:length(lake_names_lookup)) {
  # Select one lake at a time
  lake_name <- lake_names_lookup[i]
  message('Running ', lake_name)
  
  # Obtain the lake portal data (fetch, depth etc.)
  lakes_portal_subset <- lakes_portal_df |> 
    filter(str_detect(NAME, lake_name))
  
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
  setup_FLake(lake_name,
              latitude, longitude,
              elev, lake_depth, lake_fetch,
              lake_lightext = 'turbid',
              calc_cc = F, make_met = F, use_annual = F,
              outputfile = file.path(lake_name, paste0(lake_name,'_turbid_obs.rslt')))
  
  # run FLake
  run_FLake(lake_name)
  
}
