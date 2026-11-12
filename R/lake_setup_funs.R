#--------------------------------------#
## Project: Pathway optimisation - Lake District
## Script purpose: functions for setting lake specific parameters (morphology, 
#   SAGIS, FLAKE, met) - updates DATM file
## Date: 2026-11-12
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#
# Contents -----------------------------
# report_pathway()
# log_to_long()
# summarise_achievement()#' Function that updates the lDATM file using lake specific variables for LD modelling
#'
#' @param WBID_use site, water basin ID
#' @param lDATM_SETTINGS_use, the lDATM list 
#' @param change_sets, which sets to change
#'
#' @returns
#' @export
#'
#' @examples
set_lake_specifics <- function(WBID_use, lDATM_SETTINGS_use, change_sets) {
  source('R/met_funs.R')
  
  lakes_portal_use <- read_csv('data/lakes4PCLake.csv', show_col_types = F) |> 
    mutate(RET_TIMEyrs = ifelse(WBID == 29222, # this is Elterwater, something weird with the WRT, set to 20 days, 20/365 = 0.055 yrs
                                0.055, RET_TIMEyrs)) |> 
    filter(WBID == WBID_use)
  
  lake_id <- lakes_portal_use$WBID
  lake_name <- lakes_portal_use$NAME
  
  years_run <- lDATM_SETTINGS_use$run_settings['dReady', 'Set2']
  
  #-------------------------------------------------------------#
  # Forcing timeseries -----------------------------------------#
  # ------------------------------------------------------------#
  
  #-----------------------------------------------------#
  ## Extract and set WB specific parameters ---------------------
  #-----------------------------------------------------#
  latitude_use <- lakes_portal_use |> 
    select(WBLAT) |> pull()
  
  longitude_use <- lakes_portal_use |> 
    select(WBLONG) |> pull()
  
  Qin_use <- lakes_portal_use |> 
    mutate(RET_TIMEdays = RET_TIMEyrs * 365) |> # this value is in years, convert to days
    mutate(QIn = (MNDP*1000)/RET_TIMEdays) |> # convert MNDP to mm
    select(QIn) |> 
    pull() 
  
  # Set WB specific parameters -------------------------#
  #-----------------------------------------------------#
  lDATM_SETTINGS_use$params['cFetch', change_sets] <- sqrt(lakes_portal_use$WBSAREA*10000) # convert from Ha->m2 (typical fetch)
  lDATM_SETTINGS_use$params['cLAT', change_sets] <- latitude_use # latitude
  lDATM_SETTINGS_use$params['cDepthWInit0', change_sets] <- lakes_portal_use$MNDP # mean depth
  
  #-----------------------------------------------------#
  # Qin--------------------------------------------------
  #-----------------------------------------------------#
  
  Qin_use <- lakes_portal_use |> 
    mutate(RET_TIMEdays = RET_TIMEyrs * 365) |> # this value is in years, convert to days
    mutate(QIn = (MNDP*1000)/RET_TIMEdays) |> # convert MNDP to mm
    select(QIn) |> 
    pull() 
  
  # scale the Qin using teh seasonal scaling factors
  discharge_scaling <- read_csv('data/EArivers_archive/month_scaling_factors.csv', show_col_types = F) |> 
    distinct(month, discharge_scaling) |>  # only take the discharge scaling
    mutate(Qin_scaled = Qin_use * discharge_scaling) 
  
  daily_Qin <- discharge_scaling |> 
    mutate(Date = dmy(paste0('01-',month, '-2026')),
           n_days = days_in_month(Date)) |> 
    select(Date, Qin_scaled) |> 
    full_join(data.frame(Date = as_date(seq.Date(as_date('2026-01-01'), as_date('2026-12-31'), 'day'))),
              by = join_by(Date)) |> 
    arrange(Date) |> 
    mutate(Qin_scaled = zoo::na.locf(Qin_scaled),
           day = yday(Date))
  
  # Repeat Qin
  pclake_Qin <- daily_Qin |>  
    select(day, Qin_scaled) |> 
    reframe(.by = all_of(c('day')),
            value = rep(Qin_scaled, years_run),
            year = 1:years_run) |> 
    arrange(year, day) |> 
    mutate(time = row_number())
  
  lDATM_SETTINGS_use$forcings$sSet2$mQInEpi$value <- pclake_Qin$value[c(1:nrow(pclake_Qin), nrow(pclake_Qin))] 
  lDATM_SETTINGS_use$forcings$sSet2$mQInHyp$value <- 0
  
  lDATM_SETTINGS_use$forcings$sSet3$mQInEpi$value <- pclake_Qin$value[c(1:nrow(pclake_Qin), nrow(pclake_Qin))] 
  lDATM_SETTINGS_use$forcings$sSet3$mQInHyp$value <- 0
  
  # -----------------------------------------------------#
  ## SAGIS (nutrient loading) ---------------------
  #------------------------------------------------------#
  sagis_loads <- read_csv(file.path(project_location, 'data/FW_ Request for information - Ref_ EIR2026_11092GC/lake_loads_daily.csv'), show_col_types = F) |> 
    filter(WBID == lake_id) |> 
    full_join(select(lakes_portal_use, WBID, WBSAREA), by = join_by(WBID)) |> 
    mutate(value_kg_m2_day = value_kg_day/(WBSAREA*10000), # convert from Ha to m2
           value_g_m2_day = value_kg_m2_day * 1000) # convert from kg to g
  
  # Repeat nutrient Loads
  loads_pclake <- sagis_loads |>  
    select(variable, day, value_g_m2_day) |> 
    reframe(.by = all_of(c('variable', 'day')),
            value = rep(value_g_m2_day, years_run),
            year = 1:years_run) |> 
    arrange(variable, year, day) |> 
    mutate(.by = variable, time = row_number())
  
  PLoad <- loads_pclake |> 
    filter(variable == 'phosphate')
  
  NLoad  <- loads_pclake |> 
    filter(variable == 'nitrate')
  
  lDATM_SETTINGS_use$forcings$sSet2$mPLoadEpi$value <- PLoad$value[c(1:nrow(PLoad), nrow(PLoad))] 
  lDATM_SETTINGS_use$forcings$sSet2$mNLoadEpi$value <- NLoad$value[c(1:nrow(NLoad), nrow(NLoad))] 
  
  lDATM_SETTINGS_use$forcings$sSet3$mPLoadEpi$value <- PLoad$value[c(1:nrow(PLoad), nrow(PLoad))] 
  lDATM_SETTINGS_use$forcings$sSet3$mNLoadEpi$value <- NLoad$value[c(1:nrow(NLoad), nrow(NLoad))] 
  # FYI; repeat the last row, for some reason that I don't know the timeseries is longer in PCLake
  
  
  # -----------------------------------------------------#
  ## FLake (water temperatures) --------------------------
  #------------------------------------------------------#
  flake_temps <- read_csv(file.path(project_location, 'data', 'flake', lake_id, paste0(lake_id, '_predictions.csv')),
                          show_col_types = F)
  flake_mixdepth <- read_csv(file.path(project_location, 'data', 'flake', lake_id, paste0(lake_id, '_MLpredictions.csv')),
                             show_col_types = F)
  flake_strat <- read_csv(file.path(project_location, 'data', 'flake', lake_id, paste0(lake_id, '_stratpredictions.csv')),
                          show_col_types = F)
  
  # Repeat temperatures
  flake_temps <- flake_temps |> 
    full_join(flake_mixdepth, by = 'doy') |> 
    full_join(flake_strat, by = 'doy') |> 
    reframe(.by = 'doy',
            across(all_of(c('Ts', 'Tb', 'h_ML', 'strat')),
                   ~rep(.x, years_run)),
            year = 1:years_run) |> 
    arrange(year, doy) |> 
    mutate(time = row_number())
  
  lDATM_SETTINGS_use$forcings$sSet2$mTempEpi$value <- flake_temps$Ts[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS_use$forcings$sSet2$mTempHyp$value <- flake_temps$Tb[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS_use$forcings$sSet2$mMixDepth$value <- flake_temps$h_ML[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS_use$forcings$sSet2$mStrat$value <- flake_temps$strat[c(1:nrow(flake_temps), nrow(flake_temps))] 
  
  lDATM_SETTINGS_use$forcings$sSet3$mTempEpi$value <- flake_temps$Ts[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS_use$forcings$sSet3$mTempHyp$value <- flake_temps$Tb[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS_use$forcings$sSet3$mMixDepth$value <- flake_temps$h_ML[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS_use$forcings$sSet3$mStrat$value <- flake_temps$strat[c(1:nrow(flake_temps), nrow(flake_temps))] 
  # FYI; repeat the last row, for some reason that I don't know the timeseries is longer in PCLake
  
  
  # -----------------------------------------------------#
  ## ERA5-Land observations ------------------------------
  #------------------------------------------------------#
  met_vars <- c('u10', 'v10', 'ssrd')
  names(met_vars) <- c('u10', 'v10', 'sr')
  
  era5_met <- map(met_vars, get_era5_ts, 
                  latitude = latitude_use, longitude = longitude_use) |> 
    reduce(full_join, by = 'date') |>
    filter(between(date, 
                   as_date('2000-01-01'), # extract only 25 years
                   as_date('2024-12-31'))) |>
    mutate(day = day(date),
           month = month(date)) |> 
    filter(!(day==29 & month ==2)) |>     
    #Need to remove the 29th February from all years
    pivot_longer(cols = all_of(unname(met_vars))) |> 
    reframe(.by = all_of(c('name', 'date')),
            value = rep(value, 2), # two repeats of the 25 years
            year = 1:2) |> 
    arrange(year, date) |> 
    pivot_wider(names_from = name, values_from = value) |> 
    rename(!!!met_vars) |> 
    mutate(time = row_number(), 
           ws = sqrt(u10^2 + v10^2)) 
  
  
  lDATM_SETTINGS_use$forcings$sSet2$mVWind$value <- era5_met$ws[c(1:nrow(era5_met), nrow(era5_met))][c(1:((years_run*365)+1))] 
  lDATM_SETTINGS_use$forcings$sSet2$mLOut$value <- era5_met$sr[c(1:nrow(era5_met), nrow(era5_met))] [c(1:((years_run*365)+1))]
  
  lDATM_SETTINGS_use$forcings$sSet3$mVWind$value <- era5_met$ws[c(1:nrow(era5_met), nrow(era5_met))][c(1:((years_run*365)+1))] 
  lDATM_SETTINGS_use$forcings$sSet3$mLOut$value <- era5_met$sr[c(1:nrow(era5_met), nrow(era5_met))] [c(1:((years_run*365)+1))]
  # FYI; repeat the last row, for some reason that I don't know the timeseries is longer in PCLake
  
  return(lDATM_SETTINGS_use)
}
