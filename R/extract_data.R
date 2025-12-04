#--------------------------------------#
## Project: Lake District nutrient critical limit estimation
## Script purpose: generate a dataframe that includes the details needed for each specific lake bifurcation
## Date: 2025-10-20
## Author: Freya Olsson
#--------------------------------------#

library(tidyverse)
setwd(here::here())
# Required lake specific parameters ------------------
# Lake specific parameters that we need:
# - mean depth (MNDP)
# - mixed depth - or assign a sine curve
# - water temperatures - use default values
# - fetch
# - light ts
# - wind ts
#-----------------------------------------------------#
lakes_portal <- readxl::read_xlsx('data/Lake District_UKCEH Portal data_raw.xlsx', sheet = 'Combined')
lakes_portal_rt <-  read_csv('data/Lake District_UKCEH Portal RT_data.csv') |> # this is extracted from the Lake District_UKCEH Portal data_raw.xlsx sheet RT Data
  rename(DISCHARGE_M3Y = `DISCHARGEm3/y`)

# need to have all of these
lakes_portal |>
  full_join(lakes_portal_rt, by = 'WBID') |>
  filter(!is.na(MNDP),
         !is.na(FETCH_KM),
         !is.na(DISCHARGE_M3Y)) |>
  write_csv('data/lakes4PCLake.csv')

# Meteorological data --------------------------------
#' which_grid
#'
#' @param x value to look up
#' @param vector vector that contains the grid values to be indexed 
#'
#' @returns index of vector that contains x
#' @export
#'
#' @examples will find the floor but to the resolution of the data e.g. x = 54.344, will take index 2 of vector = c(54.15, 54.25, 54.35, 54.45)
which_grid <- function(x, vector) {
  which(vector > x)[1] -1
}


#' get_EOBS_ts
#'
#' @param nc_file path and file name of nc_file to extract grid value from
#' @param var_name variable name that can be found in the nc_file
#' @param latitude 
#' @param longitude 
#'
#' @returns a timeseries of values
#' @export
#'
#' @examples

get_EOBS_ts <- function(nc_file, var_name, latitude, longitude) {
  library(ncdf4)
  library(tidyverse)
  library(data.table)
  
  # Open the .nc file
  our_nc_data <- nc_open(nc_file)
  time <- ncvar_get(our_nc_data, "time")
  nt <- dim(time) # how long is the time series
  
  lat_index <- which_grid(latitude, our_nc_data$dim$latitude$vals)
  lon_index <- which_grid(longitude, our_nc_data$dim$longitude$vals)
  
  var_df <- data.frame(date = as.Date(time, origin = "1950-01-01")) |>  ## make sure this is right!
    mutate("{var_name}" := ncvar_get(our_nc_data, varid = var_name, 
                                     start = c(lon_index,lat_index, 1), 
                                     count = c(1,1, nt)))
  
  return(var_df)
}
