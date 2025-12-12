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

## 1. Directory settings ---------------------------------------------------------
## using relative paths in which the project and script is saved in the work_cases
## "scripts" contains only the PCLake functions

project_location <- here::here()
dirHome <- str_split(project_location,  "(?=PCModel1350)", simplify = T)[1,1]	# location of the PCModel1350 folder
dirShell <- str_split(project_location,  "(?<=PCShell)", simplify = T)[1,1]	#  PCShell folder path
dirCpp_root <- list.dirs(dirHome)[which(str_detect(list.dirs(dirHome),"3.01/PCLake_plus"))] # location of C++ code
nameWorkCase <- tail(str_split_1(project_location, "/"), n = 1) # workcase name
fileDATM <- list.files(list.dirs(dirHome)[which(str_detect(list.dirs(dirHome), "PCLake\\+/6.13.16"))], "PL613162PLUS_bifurcation_testing.xls", full.names = T)
dirSave <- dirShell

## load external functions from the scripts folder
source(file.path(dirShell, "scripts", "R_system", "functions.R"))
source(file.path(dirShell, "scripts", "R_system", "functions_PCLake.R")) 
# use the functions from the pathway optimisation workcase for now
source(file.path(dirShell, "work_cases", "PCLake_pathway_optimisation", "scripts/optim_functions.R"))
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
# - mixed depth - or assign a sine curve
# - water temperatures - use default values
# - fetch
# - light ts
# - wind ts


# lakes_portal <- readxl::read_xlsx('data/Lake District_UKCEH Portal data_raw.xlsx', sheet = 'Combined')
# lakes_portal_rt <-  read_csv('data/Lake District_UKCEH Portal RT_data.csv') |> # this is extracted from the Lake District_UKCEH Portal data_raw.xlsx sheet RT Data
#   rename(DISCHARGE_M3Y = `DISCHARGEm3/y`)
# 
# # need to have all of these
# lakes_portal |>
#   full_join(lakes_portal_rt, by = 'WBID') |>
#   filter(!is.na(MNDP),
#          !is.na(FETCH_KM),
#          !is.na(DISCHARGE_M3Y)) |>
#   write_csv('data/lakes4PCLake.csv')


# the lakes portal data has all the basic info we need
lakes_portal_df <- read_csv('data/lakes4PCLake.csv')

# LakeIDs to loop through
lakeIDs <- read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') |> 
  filter(!is.na(LAKE_LakesTour2021_Zooplankton.csv),
         !str_detect(LAKE_LakesTour2021_Zooplankton.csv, 'Loughrigg'))   # remove Loughrigg because it's not right

lake_names_lookup <- str_extract(lakeIDs$LAKE_LakesTour2021_Zooplankton.csv, "....")
names(lake_names_lookup) <- lakeIDs$`WBID_Lake District_UKCEH Portal data_raw.xlsx`

names(lake_names_lookup)[which(lake_names_lookup == 'Wind')] <- 29233 
# the NBAS and SBAS have seperate WBIDs but the one from the lakes portal has the combined one which is different

# use these functions:
source('R/extract_data.R')

for (i in 1:length(lake_names_lookup)) {
  
  # Select one lake at a time
  lake_name <- lake_names_lookup[i]
  
  # Obtain the lake portal data (fetch, depth etc.)
  lakes_portal_subset <- lakes_portal_df |> 
    filter(str_detect(NAME, lake_name))
  
  # Forcing timeseries ---------------------------------
 
  ## Example using Elterwater
  latitude_use <- lakes_portal_subset |> 
    select(WBLAT) |> pull()
  
  longitude_use <- lakes_portal_subset |> 
    select(WBLONG) |> pull()
  
  ## E-OBS (met) ------------------------------------------
  # extract the data from the E-OBS files
  qq_files <- list.files('data/E-OBS', pattern = 'qq', full.names = T)
  qq_ts <- map(qq_files, get_EOBS_ts, var_name = 'qq', latitude = latitude_use, longitude = longitude_use) |> 
    list_rbind() |> 
    filter(between(date, 
                   as_date('1998-01-01'),
                   as_date('2024-01-01')),
           yday(date) != 366)
  
  fg_files <- list.files('data/E-OBS', pattern = 'fg', full.names = T)
  fg_ts <- map(fg_files, get_EOBS_ts, var_name = 'fg', latitude = latitude_use, longitude = longitude_use) |> 
    list_rbind()|> 
    filter(between(date, 
                   as_date('1998-01-01'),
                   as_date('2024-01-01')),
           yday(date) != 366) # remove leap year days
  

  fg_ts |> 
    # repeat the first year
    slice(rep(1:365, each = 25)) |> 
    group_by(date) |> 
    mutate(rep = row_number()) |> ungroup() |> 
    arrange(rep, date) |> 
    # combine with the original timeseries
    bind_rows(fg_ts) |> 
    mutate(dTime = row_number()-1,
           dValue = as.numeric(fg)) |>
    select(dTime, dValue) |>

    write_delim(file = 'input/mVWind.txt', delim = '\t', eol = '\t-1\n', append = F)
  
  data.frame('-1') |>
    write_delim(file = 'input/mVWind.txt', delim = '\t', append = T)
  
  qq_ts  |> 
    # repeat the first year
    slice(rep(1:365, each = 25)) |> 
    group_by(date) |> 
    mutate(rep = row_number()) |> ungroup() |> 
    arrange(rep, date) |> 
    # combine with the original timeseries
    bind_rows(qq_ts) |> 
    mutate(dTime = row_number()-1,
           dValue = as.numeric(qq)) |> 
    select(dTime, dValue) |>  
    na.omit() |> 
    write_delim(file = 'input/mLOut.txt', delim = '\t', eol = '\t-1\n', append = F)
  
  data.frame('-1') |>
    write_delim(file = 'input/mLOut.txt', append = T)
  
  ## SAGIS-SIMCAT (nutrient loading) -----------------------
  nut_conc <- read_csv('data/SAGIS/data/InflowNutrientConcs_SAGIS.csv', show_col_types = F) |> 
    filter(str_detect(NAME_short, lake_name)) |> 
    group_by(Nutrient) |> 
    summarise_if(is.numeric, mean)
  
  ### Calculate the loads -----------------------------------
  # use the conc from SAGIS and the lakes port RT (as discharge)
  nut_loads <- lakes_portal_subset |> 
    select(WBID, 
           RET_TIMEyrs,# residence time in years
           WBSAREA, # area in Ha
           VOL, # volume in m3
           DISCHARGE_M3Y) |> 
    full_join(nut_conc, by = join_by('WBID')) |> # combine with the nutrient data
    mutate(DISCHARGE_M3D = DISCHARGE_M3Y / 365, 
           AREA_M2 = WBSAREA * 10000,
           LOAD_GD = UCLimMnCon * DISCHARGE_M3D,
           AREALLOAD_GDM2 = LOAD_GD/AREA_M2)
  
  PLoad <- nut_loads |> filter(Nutrient == 'P') |> pull(AREALLOAD_GDM2) 
  # ignore N for now and just scale with P
  
  # Lake specific parameters ------------------------------
  change_sets <- which(colnames(lDATM_SETTINGS$params) %in% c('sSet2', 'sSet3'))
  #-----------------------------------------------------#
  lDATM_SETTINGS$params[which(rownames(lDATM_SETTINGS$params) == 'cFetch'), change_sets] <- lakes_portal_subset$FETCH_KM * 1000 # fetch, convert from km to m - is this actually a good estimate
  lDATM_SETTINGS$params[str_detect(rownames(lDATM_SETTINGS$params), 'cLAT'), change_sets] <- lakes_portal_subset$WBLAT # latitude
  lDATM_SETTINGS$params[str_detect(rownames(lDATM_SETTINGS$params), 'cDepthWInit0'), change_sets] <- lakes_portal_subset$MNDP # mean depth
  lDATM_SETTINGS$forcings$sSet2$mPLoadEpi$value <- PLoad
  lDATM_SETTINGS$run_settings[which(rownames(lDATM_SETTINGS$run_settings) == "dReady"),] <- 50
  
  # report the validation states
  val_states <- c('oChlaEpi', 'oO2WEpi', 'oPTotWEpi', 'oNTotWEpi', 'aSecchiT')
  lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% val_states)] <- 1 
  # plus some other states to be plotted
  diag_states <- c('uTmEpi', 'uTmHyp', 'uDepthMixMeas', 'uPLoadEpi', 'uNLoadEpi')
  lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% diag_states)] <- 1 # report these in the output
  
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
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 7.  Initialize model  ------------------------------------------------
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Make all initial states according to the run settings
  InitStates <- PCModelInitializeModel(lDATM = lDATM_SETTINGS,
                                       dirSHELL = dirShell,
                                       nameWORKCASE = nameWorkCase)
  
  
  # Run PCLake+ --------------------------------------------------
  setwd(here())
  # readxl::read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') |> 
  #   filter(!is.na(LAKE_Lakes_Tour_Chem_TeOx.xlsx)) |> 
  #   glimpse() 
  
  run1 <- PCmodelSingleRun(lDATM = lDATM_SETTINGS,
                           nRUN_SET = 2,
                           dfSTATES = InitStates,
                           integrator_method = "rk45ck",
                           dirHOME = dirHome,
                           nameWORKCASE = nameWorkCase) 
  
  
  # 
  # run1 |> 
  #   select(any_of(c('time', val_states, diag_states))) |> 
  #   pivot_longer(cols = any_of(c(val_states, diag_states)),
  #                names_to = 'variable', values_to = 'value') |> 
  #   ggplot(aes(x=time,y=value)) +
  #   geom_line()+
  #   facet_wrap(~variable, scales = 'free_y')
  
  dir.create(file.path(project_location, 'output'), showWarnings = F)
  write_csv(run1, 
            file = file.path(project_location, 'output', 
                             paste0('baseline_', lake_name, '.csv')),
            progress = F)
  
  message('Finished ', lake_name)
}
