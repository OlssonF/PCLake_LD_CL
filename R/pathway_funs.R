#--------------------------------------#
## Project: Pathway optimisation - Lake District
## Script purpose: functions for pathway implementation. Translation of measures
  # to time varying forcings, and model execution
## Date: 2026-11-12
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

# Contents -----------------------------
# forcing_key()
# ramp_fraction()
# build_ramp_path()
# get_fixed_timing()
# infer_measure_map()
# apply_ramp_to_pathways()
# run_pathway()
#

# ---------------------------------------------------------------------------
# Helpers for the timing of measures: lag (when implementation starts) and
# rate (years taken to go from the current value to the full value)
# ---------------------------------------------------------------------------

# measures that share one forcing (mPLoadEpiMult, mPLoadEpiMult2) -> same key
forcing_key <- function(x) gsub("^\\d+|\\d+$", "", x)

# Fraction of a measure implemented at time t: 0 before `start`, rising linearly
# to 1 over `duration`, 1 afterwards. duration = 0 is an instant step at `start`.
# Vectorised over t, start and duration (all in the same units, e.g. years or days).
ramp_fraction <- function(t, start, duration) {
  start    <- rep_len(start, length(t))
  duration <- rep_len(duration, length(t))
  ifelse(duration > 0,
         pmin(pmax((t - start) / duration, 0), 1),
         as.numeric(t >= start))
}


#' Daily time series of a forcing for one group of measures
#'
#' @param days    simulation days (0:n)
#' @param current value before any measure (days before the first lag)
#' @param lag_yr
#' @param rate_yr 
#' @param value   the parameter value reached at the end of the ramp
#' @details Steps are applied in order of lag. Each ramps from wherever the path is when it
#'   starts (so a step that starts before the previous one finishes simply takes over).
#'   With rate = 0 this is the original step change.
build_ramp_path <- function(days, current, lag_yr, rate_yr, value) {
  ord <- order(lag_yr)                              # stable for ties
  lag_yr <- lag_yr[ord]; rate_yr <- rate_yr[ord]; value <- value[ord]
  
  path <- rep(current, length(days))
  for (k in seq_along(value)) {
    t0  <- floor(lag_yr[k] * 365)                   # start day, as before
    dur <- rate_yr[k] * 365                         # days to reach the full value
    idx <- which(days >= t0)
    if (length(idx) == 0) next                      # starts after the end of the run
    v0 <- path[idx[1]]                              # where the path is when this step starts
    frac <- ramp_fraction(days[idx], t0, dur)
    # once fully implemented use the exact target value (no floating point drift, so rate = 0 gives
    # exactly the original step change)
    path[idx] <- ifelse(frac >= 1, value[k], v0 + (value[k] - v0) * frac)
  }
  path
}

# lags / rates fixed in the config (not optimised): data frame measure, role, value.
# They travel as an attribute of the measure_map (possible_measures).
get_fixed_timing <- function(measure_map) {
  ft <- attr(measure_map, "fixed_timing")
  if (is.null(ft)) {
    ft <- data.frame(measure = character(0), role = character(0), value = numeric(0), stringsAsFactors = FALSE)
  }
  ft
}

# Parameter roles inferred from names, used only if run_pathway() is not given a measure_map:
# <measure>_lag / <measure>_lag2 and <measure>_rate / <measure>_rate2
infer_measure_map <- function(name_pars) {
  role <- ifelse(grepl("_rate\\d*$", name_pars), "rate",
                 ifelse(grepl("lag", name_pars), "lag", "measure"))
  measure <- ifelse(role == "measure", name_pars, sub("_(lag|rate)(\\d*)$", "\\2", name_pars))
  data.frame(parameter = name_pars, role = role, measure = measure, stringsAsFactors = FALSE)
}

# Show measures as implemented in annual pathway output (plotting): held at the current value
# until the lag, then ramped over the rate. pairs: from measure_lag_pairs(cfg)
apply_ramp_to_pathways <- function(df, pairs) {
  for (k in seq_len(nrow(pairs))) {
    m   <- pairs$main[k]
    # optimised (a column in the pathway output) > fixed in the config > none (0)
    lg  <- if (!is.na(pairs$lag[k]))        df[[pairs$lag[k]]]   else if (!is.na(pairs$lag_fixed[k]))  pairs$lag_fixed[k]  else 0
    rt  <- if (!is.na(pairs$rate[k]))       df[[pairs$rate[k]]]  else if (!is.na(pairs$rate_fixed[k])) pairs$rate_fixed[k] else 0
    cur <- if (is.na(pairs$current[k])) 0 else pairs$current[k]
    df[[m]] <- cur + (df[[m]] - cur) * ramp_fraction(df$year, lg, rt)
  }
  df
}


#' Run PCLake management pathway with modified model parameters
#'
#' @param val_pars values of the pathway parameters
#' @param name_pars names of the pathway parameters
#' @param initial_conditions dataframe of initial values (from prepInitials() or use default)
#' @param current_val 
#' @param measure_map the possible_measures df, from build_possible_measures()
#'
#' @returns a dataframe of PCLake output as defined in the DATM file
#' @export
#'
#' @examples

run_pathway <- function(val_pars,
                        name_pars, 
                        current_val, 
                        initial_conditions = NULL, 
                        measure_map = NULL) {
  
  # For debugging ----------------- #
  # val_pars <- possible_measures$lower_bound
  # name_pars <- possible_measures$parameter
  # current_val <- possible_measures$current_val
  # initial_conditions <- equilibrium_states
  #--------------------------------#
  lDATM_SETTINGS_obj <- lDATM_SETTINGS
  val_pars <- as.numeric(val_pars)
  
  # ---- what is each parameter: a measure, a lag or a rate? ----
  if (is.null(measure_map)) {
    measure_map <- infer_measure_map(name_pars)
  }
  idx <- match(name_pars, measure_map$parameter)
  
  if (anyNA(idx)) {
    stop("measure_map has no entry for: ", paste(name_pars[is.na(idx)], collapse = ", "))
  }
  
  role     <- measure_map$role[idx]
  measure  <- measure_map$measure[idx]
  fixed_timing <- get_fixed_timing(measure_map)   # lags / rates fixed in the config
  is_sagis <- grepl("SAGIS", name_pars)
  
  days <- 0:(lDATM_SETTINGS_obj$run_settings["dReady", "Set2"]*365)
  
  # time series for one group of measures (e.g. mPLoadEpiMult + mPLoadEpiMult2 share a forcing)
  path_for_group <- function(g) {
    mains <- name_pars[role == "measure" & forcing_key(name_pars) == g]
    
    cur <- unique(current_val[match(mains, name_pars)])
    if (length(cur) != 1 || is.na(cur)) {
      stop("Measures ", paste(mains, collapse = ", "), " need one common, non-NA current_val.")
    }
    
    # the lag / rate belonging to each measure (0 if it has none)
    timing <- function(m, r) {
      i <- which(role == r & measure == m)
      if (length(i) > 0) return(val_pars[i[1]])                    # optimised
      fx <- fixed_timing$value[fixed_timing$measure == m & fixed_timing$role == r]
      if (length(fx) > 0) fx[1] else 0                              # fixed in the config, or none
    }
    build_ramp_path(days, cur,
                    lag_yr  = vapply(mains, timing, numeric(1), r = "lag"),
                    rate_yr = vapply(mains, timing, numeric(1), r = "rate"),
                    value   = val_pars[match(mains, name_pars)])
  }
  
  
  # -------- SAGIS parameters -----------
  # Are there any SAGIS parameters
  sagis_groups <- unique(forcing_key(name_pars[role == "measure" & is_sagis]))
  
  if (length(sagis_groups) > 0) {
    
    forcing_df <- NULL
    
    for (g in sagis_groups) {
      forcing_df_temp <- data.frame(time = days, value = path_for_group(g), sagis_soure = g)
      forcing_df <- bind_rows(forcing_df, forcing_df_temp)
    }
    
    P_df <- forcing_df |> 
      pivot_wider(names_from = sagis_soure, values_from = value) |> 
      mutate(total = 1 - (rowSums(across(all_of(sagis_groups))))) |> 
      # sum the differences for each of the SAGIS parameters
      # 1 minus to get the right multiplier
      mutate(total = ifelse(total < 0, 0.01, total)) # multiplier cannot be less than 1
    full_join(lDATM_SETTINGS_obj$forcings$sSet2$mPLoadEpi, by = 'time') |> 
      mutate(value = value * total) |> 
      # what is the change in P load
      select(time, value) |> 
      as.data.frame()
    
    rownames(P_df) <- P_df$time + 1
    
    lDATM_SETTINGS_obj$forcings$sSet2$mPLoadEpi <- P_df
    
    # also need to take out the N load forcings?
    lDATM_SETTINGS_obj$params['ReadNLoad', c('sSet2')] <- 0
    lDATM_SETTINGS_obj$forcings$sSet2$mNLoadEpi$value <- -9999
  } else {
    sagis_pars <- NULL
    sagis_pars_lag <- NULL
  }
  
  
  # -------- PCLake parameters -----------
  #  NOTE: a measure with a lag or rate only works if there is already something in the forcings and you
  # are just modifying it, otherwise you have to recompile the model
  pc_main     <- name_pars[role == "measure" & !is_sagis]
  has_timing  <- function(m) any(role %in% c("lag", "rate") & measure == m) || any(fixed_timing$measure == m)
  done_groups <- character(0)
  
  for (m in pc_main) {
    g          <- forcing_key(m)
    grp_mains  <- pc_main[forcing_key(pc_main) == g]
    timed      <- any(vapply(grp_mains, has_timing, logical(1)))
    
    if (!timed && length(grp_mains) == 1) {
      # constant for the whole run: change the forcing if there is one, otherwise the parameter
      v <- val_pars[name_pars == m]
      if (m %in% names(lDATM_SETTINGS_obj$forcings$sSet2)) {
        lDATM_SETTINGS_obj$forcings$sSet2[[m]]$value <- v
      } else {
        lDATM_SETTINGS_obj$params[m, "sSet2"] <- v
      }
      
    } else if (!(g %in% done_groups)) {
      # implemented over time via the forcings in the DATM file
      if (!(g %in% names(lDATM_SETTINGS_obj$forcings$sSet2))) {
        warning("'", g, "' is not in the DATM forcings - the model needs recompiling for it to have an effect.")
      }
      lDATM_SETTINGS_obj$forcings$sSet2[[g]] <- data.frame(time = days, value = path_for_group(g))
      done_groups <- c(done_groups, g)
    }
  }
  
  # Initialise and run model with these parameters
  
  if (!is_null(initial_conditions)) {  
    lDATM_SETTINGS_obj$states[initial_conditions$variable, "sDefaultSetTurbid0"] <- initial_conditions$value 
    # these could be extracted from a equilibrium model run
    message("Using provided initial conditions not DATM file values")
  }
  # or just reinitialise based on the DATM file
  
  InitStates_01 <- PCModelInitializeModel(lDATM = lDATM_SETTINGS_obj,
                                          dirSHELL = dirShell,
                                          nameWORKCASE = nameWorkCase)
  
  PCModel_run <- PCmodelSingleRun(lDATM = lDATM_SETTINGS_obj,
                                  nRUN_SET = 2, # turbid stratified
                                  dfSTATES = InitStates_01,
                                  integrator_method = "lsodes",#rk45ck
                                  dirHOME = dirHome,
                                  nameWORKCASE = nameWorkCase)
  return(PCModel_run)
  
}