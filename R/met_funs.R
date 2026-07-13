#--------------------------------------#
## Project: Lake District nutrient critical limit estimation
## Script purpose: generate a dataframe that includes the details needed for each specific lake bifurcation
## Date: 2025-10-20
## Author: Freya Olsson
#--------------------------------------#

library(tidyverse)
library(ncdf4)
library(data.table)
# setwd(here::here())
# Required lake specific parameters ------------------
# Lake specific parameters that we need:
# - mean depth (MNDP)
# - mixed depth - or assign a sine curve
# - water temperatures - use default values
# - fetch
# - light ts
# - wind ts
#-----------------------------------------------------#

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

  # Open the .nc file
  our_nc_data <- nc_open(nc_file)
  time <- ncvar_get(our_nc_data, "time")
  nt <- dim(time) # how long is the time series
  
  lat_index <- which_grid(latitude, our_nc_data$dim$latitude$vals)
  lon_index <- which_grid(longitude, our_nc_data$dim$longitude$vals)
  
  var_df <- data.frame(date = as.Date(time, origin = "1950-01-01")) |>  ## make sure this is right! our_nc_data$dim$time$units
    mutate("{var_name}" := ncvar_get(our_nc_data, varid = var_name, 
                                     start = c(lon_index,lat_index, 1), 
                                     count = c(1,1, nt)))
  
  return(var_df)
}



#' Extract timeseries from the downloaded ERA5-Land data
#'
#' @param var_name variable name in ERA
#' @param latitude 
#' @param longitude 
#' @param daily calculate the daily means?
#'
#' @returns
#' @export
#'
#' @examples
get_era5_ts <- function(var_name, latitude, longitude, daily = T) {
  
  file_var <- if (var_name %in% c('u10', 'v10')) {
    "wind"
  } else if (var_name %in% c("d2m", "t2m")) {
    "temperature"
  } else if (var_name %in% c("ssrd")) {
    "radiation"
  } else {
    stop('var_name is not in the era5 files (u10, v10, d2m, t2m, ssrd')
  }
  
  missing <- 0 # starts with no data and then checks each time
  lat_round <- formatC(floor(latitude / 0.1) * 0.1, format = "f", digits = 2)
  long_round <- formatC(floor(longitude / 0.1) * 0.1, format = "f", digits = 2)
  
  # Uses the while loop because some points are in grid cells with no data, near the coast
  while (missing == 0) {

    nc_path <- paste(lat_round, long_round, sep ="_")
    
    nc_file <- list.files(file.path("data", "era5_grid", nc_path, "nc"),
                          full.names = T)[str_detect(list.files(file.path("data", "era5_grid", nc_path, "nc")), file_var)]
    
    message("Looking for ERA5-Land at ", nc_path)
    
    if (length(nc_file) == 0) {
      stop('No downloaded ERA5-Land data is available for this location. Check that the lat-lon are right.')
    }
    
    # Open the .nc file
    our_nc_data <- nc_open(nc_file)
    time <- ncvar_get(our_nc_data, "valid_time") 
    
    nt <- dim(time) # how long is the time series
    
    var_df <- data.frame(datetime = as_datetime(time * 3600, tz = 'UTC')) |> 
      # it's in hours since 1970, convert to secs |>  
      ## make sure this is right! our_nc_data$dim$valid_time$units
      mutate("{var_name}" := ncvar_get(our_nc_data, varid = var_name))
    
    # checks if this grid has data
    missing <- var_df |> 
      filter(!if_any(var_name, is.na)) |> 
      nrow()
    
    if (missing == 0) {
      message('That grid square has no data, using adjacent')
      
    }
    
    # reset the lat-lon to look for data
    long_round <- formatC(as.numeric(long_round) + 0.1, 
                          format = "f", digits = 2) # move to the E by 0.1 and check again
  }
  
  
  if (var_name == 'ssrd') {
    # deaccumulate radiation values
    var_df <- var_df |> 
      mutate(across(all_of(var_name), ~ .x / 3600))
  }
  
  if (daily) {
    
    var_df_daily <- var_df |> 
      mutate(date = as_date(datetime)) |> 
      reframe(.by = 'date',
              across(all_of(var_name), mean))
    
    return(var_df_daily)
  } else {
    
    return(var_df)
    
  }
  
}

