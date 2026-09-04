#--------------------------------------#
## Project: Lake District nutrient critical limit estimation
## Script purpose: Baseline PCLake+ runs 
## Date: 2025-12-03
## Author: Freya Olsson
#--------------------------------------#

library(tidyverse)
library(here)
library(doSNOW)
library(parallelly)
source('R/met_funs.R')

## Order of actions to run PCLake in R
##   1. Making folder structure for running the model
##   2. Load DATM file 
##   < Make adjustments to the model > 
##   4. Make cpp files
##   5. Compile model
##   6. Run model

## 0. Settings --------------------------
##---------------------------------------#

## Global settings
options(scipen = 999) ## no scientific notation
subset <- T

## 1. Directory settings ---------------------------------------------------------
## using relative paths in which the project and script is saved in the work_cases
## "scripts" contains only the PCLake functions

project_location <- here::here()
dirHome <- str_split(project_location,  "(?=PCModel1350)", simplify = T)[1,1]	# location of the PCModel1350 folder
dirShell <- str_split(project_location,  "(?<=PCShell)", simplify = T)[1,1]	#  PCShell folder path
dirCpp_root <- list.dirs(dirHome)[which(str_detect(list.dirs(dirHome),"3.01/PCLake_plus"))] # location of C++ code
nameWorkCase <- tail(str_split_1(project_location, "/"), n = 1) # workcase name
fileDATM <- list.files(list.dirs(dirHome)[which(str_detect(list.dirs(dirHome), "PCLake\\+/6.13.16"))], "PL613162PLUS_LakeDistrict.xls", full.names = T)
dirSave <- dirShell

## load external functions from the scripts folder
source(file.path(dirShell, "scripts", "R_system", "functions.R"))
source(file.path(dirShell, "scripts", "R_system", "functions_PCLake.R")) 
# use the functions from the pathway optimisation workcase for now
source(file.path(dirShell, "work_cases", "PCLake_pathway_optimisation", "R/optim_functions.R"))
## 2. Loading DATM file ------------------
lDATM_SETTINGS <- PCModelReadDATMFile_PCLakePlus(fileXLS  = fileDATM,folderTXT = 'input',
                                                 locDATM = "excel",
                                                 locFORCING = "txt")
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 4. Optional: adjust model settings------------------------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## For example, change the sediment settings of a lake
# lDATM_SETTINGS$params <- adjustSedimentParamSettings_inclBank(lDATM_SETTINGS$params, 
#                                                               paramset = 2, 
#                                                               sediment_type = "peat")

# Adjust lake specific parameters and look up the forcings
# Required lake specific parameters ------------------
# Lake specific parameters that we need:
# - mean depth (MNDP)
# - water temperatures - from FLake
# - fetch 
# - light ts
# - wind ts

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

# report the validation states
# plus some other states to be plotted
rep_vars <- c(
  'oChlaEpi', 'oO2WEpi', 'oPTotWEpi','oPTotWHyp', 'oNTotWEpi',  'oNTotWHyp', 'aSecchiT',
  'uTmEpi', 'uTmHyp', 'uDepthMixMeas', 'uDepthWEpi', 'uDepthWHyp', 'aStrat', 
  'oPO4WHyp',  'oPO4WEpi', 'aDVeg', 'oO2WHyp',  'tPDifPO4Hyp', 'aDFiAd', 'aDFiJv', 'aDPisc', 
  'aDError','aNError','aPError','aO2Error', 'aDepthError',# model errors
  'uVWind', 'uLOut','uPLoadEpi', 'uNLoadEpi' # forcings
)
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% rep_vars)] <- 1 # report these in the output



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5.  Make and adjust cpp files ----------------------------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## The nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 2) # using set2

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 6.  Compile the model ------------------------------------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PCModelCompileModelWorkCase(dirSHELL = dirShell,
                            nameWORKCASE = nameWorkCase)


years_run <- lDATM_SETTINGS$run_settings['dReady', 'Set3']
change_sets <-c('sSet2', 'sSet3')

setwd(project_location)

# i=1
for (i in 1:length(lake_names_lookup)) {
  
  # Select one lake at a time
  lake_id <- lake_names_lookup[i]
  lake_name <- names(lake_names_lookup[i])
  
  # Obtain the lake portal data (fetch, depth etc.)
  lakes_portal_use <- lakes_portal_df |> 
    filter(WBID ==lake_id) 
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
  lDATM_SETTINGS$params['cFetch', change_sets] <- sqrt(lakes_portal_use$WBSAREA*10000) # convert from Ha->m2 (typical fetch)
  lDATM_SETTINGS$params['cLAT', change_sets] <- latitude_use # latitude
  lDATM_SETTINGS$params['cDepthWInit0', change_sets] <- lakes_portal_use$MNDP # mean depth

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
    mutate(Qin_scaled = imputeTS::na_locf(Qin_scaled, na_remaining = 'rev'),
           day = yday(Date))
  
  # Repeat nutrient Loads
  pclake_Qin <- daily_Qin |>  
    select(day, Qin_scaled) |> 
    reframe(.by = all_of(c('day')),
            value = rep(Qin_scaled, years_run),
            year = 1:years_run) |> 
    arrange(year, day) |> 
    mutate(time = row_number())
  
  lDATM_SETTINGS$forcings$sSet3$mQInEpi$value <- pclake_Qin$value[c(1:nrow(pclake_Qin), nrow(pclake_Qin))] 
  lDATM_SETTINGS$forcings$sSet3$mQInHyp$value <- 0
  
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
  
  lDATM_SETTINGS$forcings$sSet2$mPLoadEpi$value <- PLoad$value[c(1:nrow(PLoad), nrow(PLoad))] 
  lDATM_SETTINGS$forcings$sSet2$mNLoadEpi$value <- NLoad$value[c(1:nrow(NLoad), nrow(NLoad))] 
  
  lDATM_SETTINGS$forcings$sSet3$mPLoadEpi$value <- PLoad$value[c(1:nrow(PLoad), nrow(PLoad))] 
  lDATM_SETTINGS$forcings$sSet3$mNLoadEpi$value <- NLoad$value[c(1:nrow(NLoad), nrow(NLoad))] 
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
  
  lDATM_SETTINGS$forcings$sSet2$mTempEpi$value <- flake_temps$Ts[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS$forcings$sSet2$mTempHyp$value <- flake_temps$Tb[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS$forcings$sSet2$mMixDepth$value <- flake_temps$h_ML[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS$forcings$sSet2$mStrat$value <- flake_temps$strat[c(1:nrow(flake_temps), nrow(flake_temps))] 
  
  lDATM_SETTINGS$forcings$sSet3$mTempEpi$value <- flake_temps$Ts[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS$forcings$sSet3$mTempHyp$value <- flake_temps$Tb[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS$forcings$sSet3$mMixDepth$value <- flake_temps$h_ML[c(1:nrow(flake_temps), nrow(flake_temps))] 
  lDATM_SETTINGS$forcings$sSet3$mStrat$value <- flake_temps$strat[c(1:nrow(flake_temps), nrow(flake_temps))] 
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
    pivot_wider(names_from = name, values_from = value) |> 
    rename(!!!met_vars) |> 
    mutate(time = row_number(), 
           ws = sqrt(u10^2 + v10^2)) 
  
  
  lDATM_SETTINGS$forcings$sSet2$mVWind$value <- era5_met$ws[c(1:nrow(era5_met), nrow(era5_met))][c(1:((years_run*365)+1))] 
  lDATM_SETTINGS$forcings$sSet2$mLOut$value <- era5_met$sr[c(1:nrow(era5_met), nrow(era5_met))] [c(1:((years_run*365)+1))]
  
  lDATM_SETTINGS$forcings$sSet3$mVWind$value <- era5_met$ws[c(1:nrow(era5_met), nrow(era5_met))][c(1:((years_run*365)+1))] 
  lDATM_SETTINGS$forcings$sSet3$mLOut$value <- era5_met$sr[c(1:nrow(era5_met), nrow(era5_met))] [c(1:((years_run*365)+1))]
  # FYI; repeat the last row, for some reason that I don't know the timeseries is longer in PCLake
  
  
  # Run PCLake+ --------------------------------------------------
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 7.  Initialize model  ------------------------------------------------
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Make all initial states according to the run settings
  InitStates <- PCModelInitializeModel(lDATM = lDATM_SETTINGS,
                                       dirSHELL = dirShell,
                                       nameWORKCASE = nameWorkCase)
  
  message('Running turbid initials')
  run_turbid <- PCmodelSingleRun(lDATM = lDATM_SETTINGS,
                                 nRUN_SET = 2,
                                 dfSTATES = InitStates,
                                 integrator_method = "rk45ck", #euler
                                 dirHOME = dirHome,
                                 nameWORKCASE = nameWorkCase) 
  message('Running clear initials')
  run_clear <- PCmodelSingleRun(lDATM = lDATM_SETTINGS,
                                nRUN_SET = 3,
                                dfSTATES = InitStates,
                                integrator_method = "rk45ck", #euler
                                dirHOME = dirHome,
                                nameWORKCASE = nameWorkCase) 
  

  dir.create(file.path(project_location, 'output'), showWarnings = F)
  run_turbid |> 
    write_csv(file = file.path(project_location, 'output', 
                               paste0('baseline_', lake_id, '-turbid.csv')),
              progress = F)
  
  run_clear |> 
    write_csv(file = file.path(project_location, 'output', 
                               paste0('baseline_', lake_id, '-clear.csv')),
              progress = F)
  
  message('Finished ', lake_id, ' ', lake_name)
}
