#--------------------------------------#
## Project: Lake District pathways modelling
## Script purpose: Functions to run physical lake model (FLake)
## Date: 2026-02-23
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
source('R/extract_data.R')
dir.create('data/FLake', showWarnings = F)

#' Set up the FLake simulations for LD PCLake runs
#'
#' @param lake_name short code lake name
#' @param latitude 
#' @param longitude 
#' @param lake_depth mean depth 
#' @param lake_fetch 
#' @param lake_lightext numeric value or "clear" or "turbid"
#' @param met_dir 
#' @param run_dir 
#' @param elev 
#'
#' @returns
#' @export
#'
#' @examples
setup_FLake <- function(lake_name, 
                        latitude, 
                        longitude,
                        elev,
                        lake_depth,
                        lake_fetch,
                        lake_lightext = 'clear', 
                        calc_cc = T,
                        outputfile = file.path(lake_name, paste0(lake_name,'.rslt')),
                        met_dir = 'data/FLake', run_dir = 'data/FLake',
                        make_met = T) {
  
  if (!dir.exists(file.path(run_dir, lake_name))) {
    dir.create(file.path(run_dir, lake_name),
               recursive = T)
  }
  
  
  
  # Generate driver file (meteorology) ------------------------------------------
  # extract the data from the E-OBS files
  if (make_met) {
    met_vars <- c('qq', 'fg', 'tg', 'hu')
    
    met_obs <- data.frame(date = NA)
    for (var in met_vars) {
      eobs_files <- list.files('data/E-OBS', pattern = var, full.names = T)
      met_obs <- map(eobs_files, get_EOBS_ts, var_name = var, 
                     latitude = latitude, longitude = longitude) |> 
        list_rbind() |> 
        filter(between(date, 
                       as_date('1998-01-01'),
                       as_date('2024-01-01')),
               # yday(date) != 366)
        ) |>  # remove leap year days
        full_join(met_obs, by = 'date')
      
    }
    
    met_obs <- met_obs |> 
      filter(!is.na(date)) |> # fill NAs
      mutate(across(all_of(met_vars), zoo::na.approx))|> 
      # Convert relative humidity (%) and temperature (°C)
      # into actual vapor pressure (millibars, mb)
      mutate(vp = calc_vaporpressure(temp = tg, relhum = hu))
    
    if (calc_cc) {
      # Estimate cloud cover - start by calculating the clear sky for every hour
      met_obs_hourly <- tibble(datetime = as.POSIXct(seq(ymd_hm(format(met_obs$date[1], '%Y-%m-%d %H:%M')), 
                                                         ymd_hm(format(met_obs$date[nrow(met_obs)], '%Y-%m-%d %H:%M')) + hours(23), 
                                                         by = '1 hour'), format ="%Y-%m-%d %H:%M")) |>
        mutate(date = as_date(datetime))  |>
        full_join(select(met_obs, date, tg, hu), by = 'date') |> filter(!is.na(datetime)) |> 
        mutate(qq_clear = calc_clearskyrad(latitude, longitude, elev, datetime, temp = tg, relhum = hu))
      
      # compare the clear sky with observed to estimate cloud cover (as a fraction 0-1)
      met_obs1 <- met_obs_hourly |> 
        reframe(.by = date, qq_clear = mean(qq_clear)) |> 
        full_join(met_obs, by = 'date') |> 
        mutate(cc = calc_cc(clearsky = qq_clear, obs = qq))
    } else {
      met_obs1 <- met_obs |>  
        mutate(cc = 0.75)
    }
    
    
    # get the annual cycle
    met_obs_annual <- met_obs1 |> 
      mutate(doy = yday(date)) |> 
      reframe(.by = doy, across(all_of(c('hu', 'qq', 'tg', 'vp', 'fg', 'cc')), .fns = median)) |> # annual average
      
      # Generate a repeating time series
      slice(rep(1:365, each = 10)) |> 
      group_by(doy) |> 
      mutate(rep = row_number()) |> ungroup() |> 
      arrange(rep, doy) |> 
      mutate(seqnum = row_number()) |> 
      # select the met vars that go in the drivers
      # order = Sequential number	Solar Radiation (W/m2)	Air Temperature (oC)	Air Humidity (mb)	Wind Speed (m/s)	Cloudiness (0-1)
      select(seqnum, qq, tg, vp, fg, cc) |> 
      mutate(across(where(is.double), ~round(.x, 3))) # make it easier to read
    
    
    ## write the drivers file -----------------
    met_file <- file.path(run_dir, lake_name, paste0(lake_name, '_met.dat'))
    write_delim(met_obs_annual, file = met_file, col_names = F)
  }
  
  # Generate nml file ----------------------------------
  lake_nml <- glmtools::read_nml('data/FLake/example.nml')
  lake_nml$SIMULATION_PARAMS$del_time_lk <- 86400
  lake_nml$SIMULATION_PARAMS$time_step_number <- 365*10
  
  ## Meteorology -------------
  lake_nml$METEO$meteofile     = file.path(lake_name, paste0(lake_name,'_met.dat'))
  lake_nml$METEO$outputfile    = outputfile
  lake_nml$METEO$`z_wind_m(1)` = 2 # TRY THIS
  
  ## lake specific parameters -----------
  lake_nml$LAKE_PARAMS$depth_w_lk  = lake_depth     # Lake depth [m]
  lake_nml$LAKE_PARAMS$fetch_lk    = lake_fetch     # Typical wind fetch [m] 
  lake_nml$LAKE_PARAMS$sediments_on = TRUE   	    # .FALSE. if the sediments layer is switched off
  lake_nml$LAKE_PARAMS$depth_bs_lk =  4.0          # Depth of the thermally active layer of the bottom sediments [m]
  lake_nml$LAKE_PARAMS$T_bs_lk     =  6.0           # Temperature at the outer edge of the thermally active layer of the bottom sediments [C]
  lake_nml$LAKE_PARAMS$latitude_lk = latitude       # Geographical latitude [dgr]
  
  ## lake transparency ------------ 
  if (!is.numeric(lake_lightext)) {
    if (lake_lightext == 'clear') {
      lake_lightext <- 0.5
    } else if (lake_lightext == 'turbid') {
      lake_lightext <- 2
    } else {
      stop('Light extinction should be numeric or "clear" or "turbid"')
    }
  }
  
  lake_nml$TRANSPARENCY$extincoef_optic = lake_lightext   # Extinction coefficients 
  
  glmtools::write_nml(lake_nml, file = file.path(run_dir, lake_name, paste0(lake_name, '.nml')))

}

#' Run flake exe
#'
#' @param lake_name 
#' @param run_dir where is the flake model?
#'
#' @returns
#' @export
#'
#' @examples
run_FLake <- function(lake_name, run_dir = 'data/FLake') {
  start_wd <- getwd()
  setwd(run_dir)  
  system(paste0('./flake ./', lake_name, '/', lake_name, '.nml'))
  setwd(start_wd)
}


#' calculate vapour pressure from relative humidity and temperature using Tetens formula
#'
#' @param temp air temperature in celcius
#' @param relhum relative humidity (%)
#'
#' @returns
#' @export
#'
#' @examples
calc_vaporpressure <- function(temp, relhum) {
  # Saturation vapor pressure using Tetens formula (mb)
  e_s <- 6.112 * exp((17.67 * temp) / (temp + 243.5))
  
  # Actual vapor pressure
  e <- (relhum / 100) * e_s
  
  return(e)
}


#' Calculate cloud cover from daily mean solar radiation 
#' @description Uses the daily mean observed met (temp, humidity, and radiation) to estimate cloud cover
#' based on the maximum clear sky radiation on a given day/time at a specific lat/lon
#' @param lat latitude, degrees
#' @param lon longitude, degrees
#' @param elev elevation or height (m)
#' @param datetime yyyy-mm-dd HH:MM as a posixct
#' @param temp air temperature in degree C
#' @param relhum relative humidity %
#' @param solrad oberved solar radiation
#'
#' @returns numeric
#' @export
#'
#' @examples
#' calc_cc(lat, lon, elev = 45, datetime = as.POSIXct('2020-08-01 14:00'), temp = 25, relhum = 55, solrad = 600)

calc_cc <- function(clearsky, obs) {
  
    #if the simulated is 0, measured must be zero
  #conditional statement to change the measured SR data
  qq_corr <- ifelse(clearsky == 0, 0, obs)
  
  # cloud cover as a fraction of observed out of total
  cc <- 1 - (qq_corr/clearsky)
  
  #if there is no clear sky radiation
  cc <- ifelse(cc == 1 & clearsky == 0, NA, cc)
  # if cc falls outside bounds
  cc <- ifelse(cc > 1, 1, cc)
  cc <- ifelse(cc < 0, 0, cc)
  
  
  return(as.vector(cc))
}

#' Calculate maximum clear sky radiation 
#' @description Uses the observed met (temp, humidity, and radiation) to model the 
#' maximum clear sky radiation on a given day/time at a specific lat/lon
#' @param lat latitude, degrees
#' @param lon longitude, degrees
#' @param elev elevation or height (m)
#' @param datetime yyyy-mm-dd HH:MM as a posixct
#' @param temp air temperature in degree C
#' @param relhum relative humidity %
#'
#' @returns numeric
#' @export
#'
#' @examples
#' calc_clearskyrad(lat, lon, elev = 45, datetime = as.POSIXct('2020-08-01 14:00'), temp = 25, relhum = 55)

calc_clearskyrad <- function(lat, lon, elev, datetime, temp, relhum) {
  
  #extract the julian day - using the JD function from insol is the only one that seems to work?
  metjd <- insol::JD(datetime) 
  
  # calculate sun angle
  sunv <- insol::sunvector(metjd,latitude = lat,longitude = lon, timezone = 0) 
  
  zenith <- insol::sunpos(sunv)[,2] #zenith angle at each timestep
  
  # Compute direct and diffuse beam irradiance 
  Idirdif <- insol::insolation(zenith = zenith, jd = metjd, height = elev,
                               visibility = 90, RH = relhum, 
                               tempK = temp + 273.15,
                               O3 = 0.003, alphag = 0.55) #03 is ozone thickness and alphag is albedo (between 0 and 1)
  
  # modify for angle of incidence on horizontal surface (pyranometer)
  cos_inc_sfc <- sunv%*%as.vector(insol::normalvector(0,0)) ## or sum(sunv*normalvector(0,0))
  
  # set to zero values with no indicent light
  cos_inc_sfc[cos_inc_sfc<0]=0
  
  # Add direct and diffuse simulated radiation on horizontal surface
  Isim  = Idirdif[,1] * cos_inc_sfc + Idirdif[,2]
  
  
  # if the simulated clear sky radiation is less than or equal to 
  # the diffuse radiation then make equal to 0 it is night time?
  qq_clear <- ifelse(Isim[,1] <= Idirdif[,2], 0, Isim[,1])
  
  return(as.vector(qq_clear))
  
}



read_flake <- function(filename) {
  
  # Read all lines, skipping the first header line
  datastr <- readLines(filename)
  
  # Extract variable names from the first data line
  vars <- strsplit(datastr[2], "\\s+")[[1]][-1]
  
  # Read the numeric data starting from line 3
  df <- read.table(filename,
                   header = FALSE,
                   skip = 2,
                   col.names = vars)
  
  return(df)
}


plot_flake <- name <- function(filename, var = 'Ts', ylim, xlim) {
  plot <- read_flake(filename) |> 
    ggplot(aes(x = time, y = get(var))) +
    geom_line() +
    ggtitle(basename(filename)) +
    theme_bw() +
    labs(y = var) +
    coord_cartesian(ylim = ylim,
                    xlim = xlim)
  
  return(plot)
}

