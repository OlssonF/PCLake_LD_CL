#--------------------------------------#
## Project: Pathway Optimisation Framework - multi-objective
## Script purpose: Apply the optimsation framework using PCLake for multiple objectives measures for Windermere
##                  The method was extended to also include optimisation of SAGIS inputs
## Date: 2026-01-10
## Author: Freya Olsson
#--------------------------------------#

library(tidyverse)
library(here)
library(DEoptim)
library(doSNOW)
library(parallelly)
library(ggh4x)
library(yaml)

## 0. Settings --------------------------
##---------------------------------------#

## Global settings
options(scipen = 999) ## no scientific notation

## Optimisation configuration (run settings, targets, measures, optimiser) is in a YAML file
source(here('R', 'config_funs.R'))
config_path <- here('input', 'optimisation_config.yml')
config <- load_config(config_path)

save_output <- config$run$save_output
pathway_name <- config$run$pathway_name
run_initials <- config$run$run_initials
archive_initials <- F

## 1. Directory settings ---------------------------------------------------------
## using relative paths in which the project and script is saved in the work_cases
## "scripts" contains only the PCLake functions

project_location <- here()
DATM_file <- "PL613162PLUS_LakeDistrict.xls"

dirHome <- str_split(project_location,  "(?=PCModel1350)", simplify = T)[1,1]	# location of the PCModel1350 folder
dirShell <- str_split(project_location,  "(?<=PCShell)", simplify = T)[1,1]	#  PCShell folder path
dirCpp_root <- list.dirs(dirHome)[which(str_detect(list.dirs(dirHome),"3.01/PCLake_plus"))] # location of C++ code
nameWorkCase <- tail(str_split_1(project_location, "/"), n = 1) # workcase name
fileDATM <- file.path(str_remove(dirShell, "/[^/]+$"), DATM_file)
dirSave <- dirShell
# ----------------------------------------------------------------------------- #

## load external functions from the scripts folder
source(file.path(dirShell, "scripts", "R_system", "functions.r"))
source(file.path(dirShell, "scripts", "R_system", "functions_PCLake.r")) 

source('R/temp funs.R')
source('R/sagis_funs.R')

## Order of actions to run PCLake in R
##   1. Making folder structure for running the model
##   2. Load DATM file 
##   < Make adjustments to the model > 
##   4. Make cpp files
##   5. Compile model
##   6. Run model

## For the optimisation of pathways it is more like:
##   1. Making folder structure for running the model - this is likely already done as I am using a project structure
##   2. Load DATM file 
##   < Make adjustments to the model > 
##   4. Make cpp files
##   5. Compile model
##   6. Run optimisation - runs in parallel kind of in a forloop
#   a) dataframe of parameter ranges (model parameters) and lags
#   b) define objective function (comparison of desired future and "current" state)
#   c) identify top pathways

## 2. Load DATM file  -------------------             
lDATM_SETTINGS <- PCModelReadDATMFile_PCLakePlus(fileXLS = fileDATM,
                                                 locDATM = "excel",
                                                 readAllForcings = F)
##----------------------------------------#

## Modifications can be made to the DATM file here! (those that need to be compiled)
## Might be a good idea to make sure the default lake parameters are loaded

# Report restart variables
restart_states <- read_table(file.path(project_location,'input', 'restart_states.txt'), col_names = 'state', show_col_types = F)
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% restart_states$state)] <- 1 # report these in the output

# report the validation states
# plus some other states to be plotted
rep_vars <- c(
  'oChlaEpi', 'oO2WEpi', 'oPTotWEpi','oPTotWHyp', 'oNTotWEpi',  'oNTotWHyp', 'aSecchiT', "oChlaBlueEpi",
  'uTmEpi', 'uTmHyp', 'uDepthMixMeas', 'uDepthWEpi', 'uDepthWHyp', 'aStrat', "aDSubVeg", "aDFish",
  'oPO4WHyp',  'oPO4WEpi', 'aDVeg', 'oO2WHyp',  'tPDifPO4Hyp', 'aDFiAd', 'aDFiJv', 'aDPisc', 
  'aDError','aNError','aPError','aO2Error', 'aDepthError',# model errors
  'uVWind', 'uLOut','uPLoadEpi', 'uNLoadEpi', 'uQInEpi'# forcings
)
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% rep_vars)] <- 1 # report these in the output

# Make it Windermere specific ------------------------
WBID_use <- 29233 # Windermere
lDATM_SETTINGS <- set_lake_specifics(WBID_use = WBID_use, change_sets = c('sSet2', 'sSet3'), 
                                     lDATM_SETTINGS_use = lDATM_SETTINGS)

#----------------------------------------------------#
## 3.4.A Make and adjust cpp files ----------       
#  - nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 2) # use turbid stratified case
##----------------------------------------#

## 5.A Compile model -----------------------
PCModelCompileModelWorkCase(dirSHELL = dirShell,
                            nameWORKCASE = nameWorkCase)
##----------------------------------------#
setwd(project_location)

if (run_initials) {
  ## Find equilibrium values ----------------#
  # Run the model to an equilibrium before starting the optimisation
  InitStates <- PCModelInitializeModel(lDATM = lDATM_SETTINGS,
                                       dirSHELL = dirShell,
                                       nameWORKCASE = nameWorkCase)
  
  message('Running turbid initials')
  windermere_baseline <- PCmodelSingleRun(lDATM = lDATM_SETTINGS,
                                          nRUN_SET = 2,
                                          dfSTATES = InitStates,
                                          integrator_method = "lsodes", #euler, #rk45ck
                                          dirHOME = dirHome,
                                          nameWORKCASE = nameWorkCase) 
  
  # extract the restart variables from the end of the baseline run
  equilibrium_states <- prepInitials(listPCModelRun = windermere_baseline, 
                                     day =  lDATM_SETTINGS$run_settings['dReady','Set2'] * 365)
  
  if (archive_initials) {
    write_csv(equilibrium_states, file.path(project_location, 'input/initialised_states/initials_windermere.csv'))
  }
  
} else {
  equilibrium_states <- read_csv(file.path(project_location, 'input/initialised_states/initials_windermere.csv'))
}

## 6.A Set optimisation values ---------------------------

### a. Define parameter values -------------------
# dataframe of parameter ranges (model parameters), lags, and rates

# Define the parameter values to be optimised (upper and lower) as well as
# the "unchanged" value (before the measure is in place) - could also be a timeseries I guess?
possible_measures <- build_possible_measures(config)


# see if there are other things required to run the measure optimisation
for (i in 1:length(possible_measures$parameter)) {
  if (!is.na(possible_measures$note[i])) {
    message(possible_measures$parameter[i], ' has requirements: ', possible_measures$note[i])
  } 
}

# check for associated parameters required for some optimisations - see possible_measures notes
lDATM_SETTINGS <- apply_fixed_params(config, lDATM_SETTINGS)


### b. Define the desired future ------------------
# What is the objective
# Ecosystem services -> indicators (PCLake states) -> targets, summaries, weights and aggregation
# all come from the 'targets' section of the config

source(file.path(project_location, "R/optim_functions.R")) 
# functions for running and evaluating the pathways

desired_es    <- build_es_config(config)
eval_settings <- build_eval_settings(config)
optim_states  <- es_states(desired_es)

# check for unknown states
unknown_states <- setdiff(optim_states, rownames(lDATM_SETTINGS$auxils))
if (length(unknown_states) > 0) {
  stop('Target states not found in the PCLake output variables: ', paste(unknown_states, collapse = ', '))
}

# how much does each indicator / state actually count in the objective?
print(weight_audit(desired_es))

## Update the DATM file and recompile the model ---------------#
# Report variables
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% restart_states$state)] <- 0 # these can be turned off
lDATM_SETTINGS$auxils$iReport[which(rownames(lDATM_SETTINGS$auxils) %in% optim_states)] <- 1 # report optim vars
lDATM_SETTINGS$params$iReport[which(rownames(lDATM_SETTINGS$params) %in% possible_measures$parameter)] <- 1 # report measure params
lDATM_SETTINGS$auxils[which(rownames(lDATM_SETTINGS$auxils) %in% 'uPLoadEpi'), ] <- 1 # also report the auxillary variable for the PLoadEpi

# make sure the pathway is for 30 years
# lDATM_SETTINGS$run_settings['dReady',] <- 30

# forcing variables --------------#
# anything that has a lag or an implementation rate needs to be in the forcings before compilation
# (SAGIS variables are not: they modify mPLoadEpi)
# (optimised lags / rates, plus any lag / rate fixed in the config)
timed_measures <- c(possible_measures$measure[possible_measures$role %in% c('lag', 'rate')],
                    get_fixed_timing(possible_measures)$measure)
timed_measures <- timed_measures[!str_detect(timed_measures, 'SAGIS')]

for (m in unique(forcing_key(timed_measures))) {
  if (is.null(lDATM_SETTINGS$forcings$sSet2[[m]])) {
    lDATM_SETTINGS$forcings$sSet2[[m]] <- data.frame(time = 0:(lDATM_SETTINGS$run_settings["dReady", "Set2"]*365),
                                                     value = -99999) 
    message('adding ', m, ' forcing')
  }
}


## 3.4. B Remake and adjust cpp files ----------       
#  - nRUN_SET determines which forcings are switched on
PCModelAdjustCPPfiles(dirSHELL = dirShell,
                      nameWORKCASE = nameWorkCase,
                      lDATM = lDATM_SETTINGS,
                      nRUN_SET = 2)
##----------------------------------------#

## 5.B Compile model -----------------------
PCModelCompileModelWorkCase(dirSHELL = dirShell,
                            nameWORKCASE = nameWorkCase)
##----------------------------------------#

# The compiled model needs to have all of the forcing variables in it (e.g. marsh area, P-loading)
# otherwise the .dll and DATM file won't match



### c. Define the objective function -----------------
# 1) run model
# 2) extract output
# 3) return a minimised value by comparing with desired future (DEOptim will make the value as negative as possible)

source(file.path(project_location, "R/optim_functions.R")) # functions for running and evaluating the pathways

#' Define the objective function with the output to be optimised
#'
#' @param val_pars a vector of parameter values to be evaluated
#' @param name_pars a vector of the parameter names, same order as val_pars
#' @param es_config ecosystem services with their indicators / targets (from build_es_config(cfg))
#' @param eval_settings evaluation settings: days, year, nonneg, bound (from build_eval_settings(cfg))
#' @description the function is split into two parts - run the pathway returns PCLake output, this is then compared with the ecosystem service targets
#' @returns an single evaluation value to be optimised
#' @export
#'
#' @examples obj_function(val_pars = val_pars, name_pars = name_pars, es_config = desired_es, eval_settings = eval_settings)
obj_function <- function(val_pars, name_pars, es_config, eval_settings) {
  
  # # For debugging ----------------- #
  # val_pars <- possible_measures$upper_bound
  # name_pars <- possible_measures$parameter
  # es_config <- desired_es
  #--------------------------------#
  
  model_output <- run_pathway(val_pars, name_pars, current_val = possible_measures$current_val,
                              measure_map = possible_measures,
                              initial_conditions = equilibrium_states)
  
  eval_output <- evaluate_pathway_es(model_output, es_config,
                                     eval_settings, # are common across ES
                                     return_detail = T)
  
  
  # log the objective value plus which ESs / indicators were met and why not
  report_pathway(val_pars, name_pars, obj_val = eval_output$total, log_dir = log_dir,
                 extra = flatten_es_result(eval_output))
  
  return(eval_output$total)
}


## 6.B Run optimisation ---------------
log_dir <- tempfile("deoptim_log_")

# where to save the output from the objective function
if (dir.exists(log_dir)) {
  unlink(log_dir, recursive = TRUE)
  message("Cleaned up log_dir: ", log_dir)
} else {
  warning("log_dir is empty or missing — log_dir NOT deleted, check ", log_dir)
}

# make temp dir
dir.create(log_dir)


## Objects / functions the workers need
cluster_export_names <- c("lDATM_SETTINGS", 'possible_measures', 'equilibrium_states',
                          "PCModelInitializeModel", 
                          "above_obj", "below_obj", "exact_obj", "range_obj",
                          "dirShell", "nameWorkCase", 'dirHome', "log_dir",
                          "PCmodelSingleRun", "RunModel", 
                          'run_pathway', 'evaluate_pathway_es', 'prepare_es_config', 'prepare_eval_settings', 'es_states',
                          'flatten_es_result', 'report_pathway')

## Parallelisation with FOREACH package
{
  nC <- parallelly::availableCores() 
  cl <- makeSOCKcluster(nC- config$run$reserve_cores)
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(deSolve))
  clusterExport(cl, cluster_export_names)
  
  doSNOW::registerDoSNOW(cl)
  
  # DEoptim settings (NP, CR, F, itermax, ...) come from the 'optimiser' section of the config
  deoptim_control <- build_deoptim_control(config, n_measures = nrow(possible_measures))
  
  # Several conditions can cause the optimization process to stop:
  # if the maximum number of iterations is reached (itermax), or
  # if a number (steptol) of consecutive iterations are unable to reduce the best function value by a certain amount (reltol * (abs(val) + reltol)). 
  # 100*reltol is approximately the percent change of the objective value required to consider the parameter set an improvement over the current best member.
  
  set.seed(1234)
  
  opt_pathway <- DEoptim::DEoptim(lower = possible_measures$lower_bound,
                                  upper = possible_measures$upper_bound, 
                                  fn = obj_function,
                                  control = deoptim_control, 
                                  name_pars = possible_measures$parameter,
                                  es_config = desired_es,
                                  eval_settings = eval_settings)
  
  parallel::stopCluster(cl)
}

# Read in the logged objective function output
log_files <- list.files(log_dir, full.names = TRUE)
log_list <- lapply(log_files, read.csv)
log_df <- do.call(rbind, log_list)
# order by write time (embedded in filename) and chunk into generationslog_df <- log_df[order(files), ]
NP <-  deoptim_control$NP
log_df$iteration <- rep(0:(nrow(log_df)/NP - 1), each = NP)


# The output of DEoptim is based on members, iterations, and populations
# iteration is a generation of a population
# one population is a collection of members
# one member is a single run of the objective function with specific parameter values
# rename all the missing names
names(opt_pathway$member$lower) <- possible_measures$parameter
names(opt_pathway$member$upper) <- possible_measures$parameter
colnames(opt_pathway$member$bestmemit) <- possible_measures$parameter
names(opt_pathway$optim$bestmem) <- possible_measures$parameter

## Visualise the best member output from each iteration
n_iter <- opt_pathway$optim$iter
iteration_summary <- as.data.frame(opt_pathway$member$bestmemit) |>
  mutate(fn_out = opt_pathway$member$bestvalit[1:n_iter])


# write output for later?------------------------------
if (save_output) {
  file.copy(config_path, file.path(project_location, 'pathways_output', paste0('config_', pathway_name, '.yml')), overwrite = TRUE) # keep the config this run used
  write_rds(c(list(possible_measures = possible_measures, desired_es = desired_es, eval_settings = eval_settings), opt_pathway$optim),
            file = file.path(project_location, 'pathways_output', paste0('summary_',pathway_name, '.RData'))) # summary of the optimisation example # summary of the optimisation example
  write_csv(iteration_summary, file = file.path(project_location, 'pathways_output',  paste0('bestmemit_',pathway_name, '.csv'))) # best pathway for each population
  
  # instead use the log
  log_df |> 
    write_csv(file = file.path(project_location, 'pathways_output',  paste0('allpops_',pathway_name, '.csv')))
}

## 7. Re-run last iteration -----------------
### 7.a Exract objective values -----------
# The DEoptim does not output the results of the function (valit) only the best member per iteration and then the parameter values
# To get the output values we can rerun the last iteration of the optimisation
last_iteration <- as_tibble(opt_pathway$member$pop)
colnames(last_iteration) <- colnames(opt_pathway$member$bestmemit)

last_iteration$runID <- rownames(last_iteration)

# set up parallelisation
n_threads <- parallel::detectCores() - 2 ## leave one free for other tasks
snowCLUSTER <- makeCluster(n_threads)
clusterExport(snowCLUSTER, c(cluster_export_names, 'obj_function')) ## here you can push certain objects to the clusters
registerDoSNOW(snowCLUSTER)

# stuff to show progress bar
pb <- txtProgressBar(0, nrow(last_iteration), style = 3)
progress <-function(n){setTxtProgressBar(pb, n)}
opts <- list(progress = progress)

# Use foreach to run each member of the last iteration in the function
last_iteration$fn_out <- foreach(i = 1:nrow(last_iteration),
                                 .combine = 'c', 
                                 .multicombine = TRUE,
                                 .packages=c("tidyverse", "deSolve"),
                                 .options.snow = opts) %dopar% {
                                   
                                   val_pars <- last_iteration[i,] |> 
                                     select(-runID) |> unlist()
                                   
                                   name_pars <- last_iteration[i,] |> 
                                     select(-runID) |> names()
                                   
                                   save <- obj_function(val_pars = val_pars, name_pars = name_pars,
                                                        es_config = desired_es, eval_settings = eval_settings)
                                 }


### 7.b Extract the indicator values -----------
# For each member of the final population: how each ES / indicator did. This is evaluated exactly as in the
# optimisation (same states, summaries, targets and window - all from the config), so what is reported is
# what was optimised.
ES_opt <- foreach(i = 1:nrow(last_iteration),
                  .combine = bind_rows, 
                  .multicombine = TRUE,
                  .packages=c("tidyverse", "deSolve"),
                  .options.snow = opts) %dopar% {
                    
                    val_pars <- last_iteration[i, possible_measures$parameter] |> unlist()
                    name_pars <- possible_measures$parameter
                    
                    model_output <- run_pathway(val_pars, name_pars, possible_measures$current_val,
                                                initial_conditions = equilibrium_states)
                    
                    res <- evaluate_pathway_es(model_output, desired_es, eval_settings, return_detail = TRUE)
                    
                    res$detail |> 
                      mutate(ID = i, objective = res$total) |> 
                      bind_cols(as_tibble_row(val_pars))
                    
                  }



if (save_output) {
  write_csv(last_iteration, file = file.path(project_location, 'pathways_output',  paste0('lastpop_',pathway_name,'.csv')))  # each member values of the final population
  
  # one row per member x ES x indicator: summary value, target, met / not met, gap, weights + the member's parameter values
  write_csv(ES_opt, file = file.path(project_location, 'pathways_output',  paste0('lastpopstate_',pathway_name, '.csv')))
}

### 7.c Extract the time series ---------
# above we extracted the last year (the year being optimised) but in this case we want to 
# look at the full pathway to acheive the final target
# (if a state is used by several indicators, its time series uses the first indicator listed for it)

# what is the "pathway" to the optimised value?
state_pathways <- foreach(i = 1:nrow(last_iteration),
                          .combine = bind_rows, 
                          .multicombine = TRUE,
                          .packages=c("tidyverse", "deSolve"),
                          .options.snow = opts) %dopar% {
                            
                            val_pars <- last_iteration[i, possible_measures$parameter] |> unlist()
                            name_pars <- possible_measures$parameter
                            
                            model_output <- run_pathway(val_pars, name_pars, possible_measures$current_val, 
                                                        initial_conditions = equilibrium_states)
                            
                            n_years <- max(floor((model_output$time - 1) / 365) + 1)
                            
                            # evaluate every simulated year with the config's summaries and window
                            yearly <- purrr::map_dfr(seq_len(n_years), \(yr) {
                              
                              d <- evaluate_pathway_es(model_output,
                                                       desired_es,
                                                       utils::modifyList(eval_settings, list(eval_year = yr)),
                                                       return_detail = TRUE)$detail |>
                                dplyr::distinct(state, .keep_all = TRUE) |> # filters out duplicated states
                                dplyr::select(state, summary_value) |> 
                                tidyr::pivot_wider(names_from = state, values_from = summary_value) |> 
                                dplyr::bind_cols(as_tibble_row(val_pars)) |> 
                                dplyr::mutate(ID = i, 
                                              year = yr)
                              
                            })
                            
                          }


if (save_output) {
  write_csv(state_pathways, file = file.path(project_location, 'pathways_output',  paste0('lastpoppathways_',pathway_name,'.csv')))  # the annual output of the final population
}


# Plotting ---------------------------------------
library(tidyverse)

# Aesthetics --------------------------------------------------------------
{labels_measures <- c(expression(atop("Vegetation removed", "(fraction)")),
                      expression(atop("Start vegetation", "removal (year)")),
                      expression(atop("Sewage P load reduction", "(%)")),
                      expression(atop("Septic tanks P load reduction", "(%)")),
                      expression(atop("Start sewage P load", "reduction (year)")),
                      expression(atop("Start ST P load", "reduction (year)")),
                      expression(atop("Agriculture P load reduction", "(%)")),
                      expression(atop("Start Ag P load", "reduction (year)")),
                      expression(atop("Average start P load", "reduction (year)")),
                      expression(atop("Total reduction in P load", "(%)")))

names(labels_measures) <- c(possible_measures$parameter, "SAGIS_total",'SAGIS_lag_mean')
labels_measures <- labels_measures[sort(names(labels_measures))]
labels_measures <- labels_measures[order(grepl("_lag", names(labels_measures)))]
}

cols_measures <- c("#6DCD59FF",
                   "#3E4A89FF","#3E4A89FF",  "#3E4A89FF","#3E4A89FF",
                   "#6DCD59FF",
                   "#3E4A89FF","#3E4A89FF",  "#3E4A89FF","#3E4A89FF")

cols_states <- c("#C7EF34FF", "#36AAF9FF", "#CB2A04FF")

labels_states <- c(oChlaEpi =expression(atop('Chlorophyll-a concentration', (mu*g~L^-1))),
                   aDSubVeg =expression(atop('Submerged vegetation biomass', (gDW~m^-2))),
                   aSecchiT   =expression(atop("Total Secchi disk depth", (m))))

labels_states_str <- c(oChlaEpi = "atop('Chlorophyll-a concentration', (mu*g~L^-1))",
                       aDSubVeg = "atop('Submerged vegetation biomass', (gDW~m^-2))",
                       aSecchiT   = "atop('Total Secchi disk depth', (m))")


last_iteration <- read_csv('pathways_output/lastpop_windermere_example.csv')
state_opt <- read_csv('pathways_output/lastpopstate_windermere_example.csv')
state_pathways <- read_csv('pathways_output/lastpoppathways_windermere_example.csv')
allpops <- read_csv('pathways_output/allpops_windermere_example.csv')

possible_measures <- read_rds('pathways_output/summary_windermere_example.RData')$possible_measures


p1 <- allpops |> 
  filter(obj == 0)  |>
  ggplot(aes(x=SAGIS_sewage *100, y = SAGIS_agriculture*100, colour = SAGIS_STLoad*100)) +
  geom_point(alpha = 0.8) +
  coord_cartesian(xlim = c(25, 40),
                  ylim = c(25, 40)) +
  scale_colour_viridis_c(limits = c(0, 10), 
                         option = 'C', begin = 1, end = 0,
                         name = 'Septic tank P load\nreduction (%)\n') +
  theme_minimal(base_size = 16) +
  labs(x= 'Sewage P load\nreduction (%)',
       y= 'Agricultural P load\nreduction (%)',
       # title = "Optimised SAGIS Load Reductions",
       subtitle = "Final population members achieving objective = 0") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  ) + 
  scale_x_continuous(breaks = seq(0,100, 5)) +
  scale_y_continuous(breaks = seq(0,100, 5))



p2 <- allpops |> 
  filter(obj == 0) |> mutate(ID = row_number()) |> 
  pivot_longer(cols = -c(ID), names_to = 'parameter', values_to = 'val_par') |> 
  right_join(select(possible_measures, parameter, upper_bound)) |> 
  mutate(of_max = (val_par/ upper_bound) * 100) |> 
  select(ID, parameter, of_max)  |> 
  pivot_wider(names_from = parameter, values_from = of_max) |>
  
  ggplot(aes(x=SAGIS_sewage, y = SAGIS_agriculture, colour = SAGIS_STLoad)) +
  geom_point(alpha = 0.8) +
  coord_cartesian(xlim = c(62.5, 100),
                  ylim = c(62.5, 100)) +
  scale_colour_viridis_c(limits = c(0, 100), 
                         option = 'C', begin = 1, end = 0,
                         name = 'Septic tank P load\nreduction efficiency\n(% of maximum)') +
  theme_minimal(base_size = 16) +
  labs(x= 'Sewage P load reduction\nefficiency (% of maximum)',
       y= 'Agricultural P load reduction\nefficiency (% of maximum)',
       # title = "Optimised SAGIS Load Reductions",
       subtitle = "Final population members achieving objective = 0") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = seq(0,100, 10)) +
  scale_y_continuous(breaks = seq(0,100, 10))


p3 <- allpops |>
  mutate(
    success = obj == 0,
    SAGIS_total = SAGIS_sewage + SAGIS_STLoad + SAGIS_agriculture,
    SAGIS_lag_mean = (SAGIS_sewage_lag + SAGIS_STLoad_lag + SAGIS_agriculture_lag) / 3
  ) |>
  ggplot(aes(
    x = SAGIS_lag_mean,
    y = SAGIS_total * 100,
    colour = obj,
    shape = success
  )) +
  geom_point(size = 3, alpha = 0.9) +
  
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 4)) +
  
  scale_colour_viridis_c(
    begin = 0.85, end = 0,
    option = "viridis", 
    name = "Objective value\n(0 = target achieved)"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    legend.position = "right",
    legend.title.position = "top",
    legend.title = element_text(face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12)
  ) +
  
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 100)) +
  
  guides(
    shape = "none",
    size = "none",
    colour = guide_colourbar(
      direction = "vertical",
      barheight = unit(6, "lines"),
      barwidth  = unit(0.6, "lines")
    )
  ) +
  
  labs(
    x = "Mean timing of P load reduction (year)",
    y = "Total reduction in P load (%)",
    title = "Relationship between total P load reductions and the timings",
    subtitle = "Colour indicates objective value; crosses mark successful solutions"
  )


## example pathways 
p4 <- state_pathways |>
  filter(ID %in% round(runif(n=4, min=1, max=80))) |>
  mutate(across(
    c('fManVeg', "SAGIS_sewage", "SAGIS_STLoad", "SAGIS_agriculture"),
    ~ ifelse(get(paste0(cur_column(), "_lag")) > year, 0, .x)
  )) |>
  select(-ends_with("lag")) |>
  select(-aDSubVeg) |>
  pivot_longer(cols = !any_of(c('ID','year')), names_to = 'opt_var') |>
  mutate(opt_var = factor(opt_var, levels = c(names(labels_states), names(labels_measures))),
         ID = paste0("Pathway ", ID)) |>
  
  ggplot(aes(
    x = year,
    y = value,
    colour = opt_var,
    linetype = opt_var
  )) +
  geom_line(linewidth = 1, alpha = 0.9) +
  
  ggh4x::facet_nested_wrap(
    vars(ID, opt_var),
    scales = 'free_y',
    nrow = 6,
    dir = 'v',
    nest_line = element_line(colour = 'grey95', linewidth = 0.6),
    strip = strip_split(
      position = c('right', 'top'),
      text_y = element_text(size = 10, face = "bold"),
      text_x = element_blank(),
      background_y = element_rect(fill = "grey55", colour = NA),
      background_x = element_blank()
    )
  ) +
  
  # Modern theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    panel.spacing.y = unit(0.4, "lines"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "grey95", colour = NA)
  ) +
  
  scale_colour_manual(
    values = c(cols_measures, cols_states),
    name = "Lake state",
    breaks = c(names(labels_measures), names(labels_states)),
    labels = c(labels_measures, labels_states)
  ) +
  
  scale_linetype_manual(
    values = c(
      "SAGIS_sewage"      = "dotted",
      "SAGIS_STLoad"      = "longdash",
      "SAGIS_agriculture" = "dotdash",
      "fManVeg"           = "solid",
      "aSecchiT"          = "solid",
      "oChlaEpi"          = "solid"
    ),
    name = "Lake state",
    breaks = c(names(labels_measures), names(labels_states)),
    labels = c(labels_measures, labels_states)
  ) +
  
  # Longer legend line segments
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(linewidth = 1, alpha = 1, size = 10)
    ),
    linetype = guide_legend(
      ncol = 1,
    )
  ) +
  
  scale_x_continuous(expand = c(0.01, 0.01)) +
  
  facetted_pos_scales(
    y = list(
      opt_var == "aSecchiT" ~ scale_y_continuous(limits = c(0.5, 2.5)),
      opt_var == "oChlaEpi" ~ scale_y_continuous(limits = c(0, 75)),
      str_detect(opt_var, 'SAGIS') ~ scale_y_continuous(limits = c(0, 0.5), n.breaks = 2)
    )
  ) +
  
  labs(
    x = "Year",
    y = "",
    title = "Four example pathways to achieve Swimming ES",
    subtitle = "Annual progression of state variables and management measures"
  )


ggsave(p1, filename = 'pathways_output/plots/p1.png', 
       width = 24, height = 12,  units = 'cm', dpi = 200)
ggsave(p2, filename = 'pathways_output/plots/p2.png', 
       width = 24, height = 12,  units = 'cm', dpi = 200)
ggsave(p3, filename = 'pathways_output/plots/p3.png', 
       width = 24, height = 12,  units = 'cm', dpi = 200)
ggsave(p4, filename = 'pathways_output/plots/p4.png', 
       width = 24, height = 12,  units = 'cm', dpi = 200)
