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