#--------------------------------------#
## Project: Lake District nutrient critical limit estimation
## Script purpose: run a bifurcation analysis for a single lake (Elterwater)
## Date: 2025-10-15
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
#-----------------------------------------------------#

# Start with Elterwater 
# the lakes portal data has all the basic info we need
lakes_portal_df <- read_csv('data/lakes4PCLake.csv')

lakes_portal_subset <- lakes_portal_df |> 
  filter(str_detect(NAME, 'Elterwater'))
  
# Lake specific parameters
change_sets <- which(colnames(lDATM_SETTINGS$params) %in% c('sSet2', 'sSet3'))

lDATM_SETTINGS$params[str_detect(rownames(lDATM_SETTINGS$params), 'cFetch$'), change_sets] <- lakes_portal_subset$FETCH_KM * 1000 # fetch, convert from km to m - is this actually a good estimate
lDATM_SETTINGS$params[str_detect(rownames(lDATM_SETTINGS$params), 'cLAT'), change_sets] <- lakes_portal_subset$WBLAT # latitude
lDATM_SETTINGS$params[str_detect(rownames(lDATM_SETTINGS$params), 'cDepthWInit0'), change_sets] <- lakes_portal_subset$MNDP # mean depth

# Report restart variables
restart_states <- read_table(file.path(project_location, 'restart_states.txt'), col_names = 'state', show_col_types = F)
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% restart_states$state)] <- 1 # report these in the output

# report the used P Load
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) == 'uPLoadEpi')] <- 1 # report these in the output
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) == 'uDepthMixMeas')] <- 1 # report these in the output
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) == 'oChlaEpi')] <- 1 # report these in the output
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5.  Make and adjust cpp files ----------------------------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## The nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 3) # using set3

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


# Bifurcation analysis --------------------------------------------------
# Using a for loop to run the model with a set value of nutrient loading to identify tipping points.
# Each nutrient load is run from a clear state and for a turbid state until equilibrium.
# Then extract the states to be evaluated (summer chla ?) and combine

mPLoadEpi <- seq(0.0001, 0.01, 0.0001) # we will loop of these values

# set to a high value
lDATM_SETTINGS$params$sSet2[which(rownames(lDATM_SETTINGS$params) == 'mPLoadEpi')] <- max(mPLoadEpi) 
# set to a low value
lDATM_SETTINGS$params$sSet3[which(rownames(lDATM_SETTINGS$params) == 'mPLoadEpi')] <- min(mPLoadEpi) 

# run to equilibrium for clear states?
PCModel_run_baseline1 <- PCmodelSingleRun(lDATM = lDATM_SETTINGS,
                                          nRUN_SET = 2,
                                          dfSTATES = InitStates,
                                          integrator_method = "rk45ck",
                                          dirHOME = dirHome,
                                          nameWORKCASE = nameWorkCase)

# extract the restart variables from the end of the baseline run
turbid_state <- prepInitials(listPCModelRun = PCModel_run_baseline1, 
                             day =  lDATM_SETTINGS$run_settings['dReady','Set2'] * 365)


# run to equilibrium for turbid states?
PCModel_run_baseline2 <- PCmodelSingleRun(lDATM = lDATM_SETTINGS,
                                          nRUN_SET = 3,
                                          dfSTATES = InitStates,
                                          integrator_method = "rk45ck",
                                          dirHOME = dirHome,
                                          nameWORKCASE = nameWorkCase)

# extract the restart variables from the end of the baseline run
clear_state <- prepInitials(listPCModelRun = PCModel_run_baseline2, 
                            day =  lDATM_SETTINGS$run_settings['dReady','Set3'] * 365)



n_threads <- parallel::detectCores() - 2 ## leave one free for other tasks
snowCLUSTER <- makeCluster(n_threads)
clusterExport(snowCLUSTER, list("lDATM_SETTINGS", 
                                "PCModelInitializeModel", 
                                "dirShell", "nameWorkCase", 'dirHome',
                                "PCmodelSingleRun", "RunModel", 'run_pathway', 
                                'turbid_state', 'clear_state', 'mPLoadEpi'))

registerDoSNOW(snowCLUSTER)

# stuff to show progress bar
pb <- txtProgressBar(0,length(mPLoadEpi), style = 3)
progress <-function(n){setTxtProgressBar(pb, n)}
opts <- list(progress = progress)


# Use foreach to run each P loading, first from a turbid state and then from clear
bifurcation_results <- foreach(i = 1:length(mPLoadEpi),
                               .combine = 'bind_rows', 
                               .multicombine = TRUE,
                               .packages=c("tidyverse", "deSolve"),
                               .options.snow = opts) %dopar% {
                                 
                                 # Run PCLake from turbid and clear conditions
                                 res1 <- run_pathway(val_pars = mPLoadEpi[i],
                                                     name_pars = 'mPLoadEpi',
                                                     initial_conditions = turbid_state)
                                 res2 <- run_pathway(val_pars = mPLoadEpi[i],
                                                     name_pars = 'mPLoadEpi',
                                                     initial_conditions = clear_state)
                                 
                                 # sumamrise the results
                                 bind_rows(res1 |>
                                             mutate(year = floor((time-1)/365) + 1,
                                                    doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
                                             filter(year == max(year), # filters to summer in the last year of the simulation
                                                    doy %in% 121:244) |> 
                                             select(oChlaEpi, uPLoadEpi) |> 
                                             summarise(across(any_of(everything()), mean)) |> 
                                             mutate(init_state = 'turbid'),
                                           
                                           
                                           res2 |>
                                             mutate(year = floor((time-1)/365) + 1,
                                                    doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
                                             filter(year == max(year), # filters to summer in the last year of the simulation
                                                    doy %in% 121:244) |> 
                                             select(oChlaEpi, uPLoadEpi) |> 
                                             summarise(across(any_of(everything()), mean)) |> 
                                             mutate(init_state = 'clear')
                                 )
                               }

parallel::stopCluster(snowCLUSTER)

ggplot(bifurcation_results, aes(x=uPLoadEpi, y=oChlaEpi, colour = init_state)) + geom_point() + geom_line()
