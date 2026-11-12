#--------------------------------------#
## Project: Pathway optimisation - Lake District
## Script purpose: functions to build and configure measures, desired ESs and
#  other things that go into the pathway optimisation. Parses the optimisation.yaml
## Date: 2026-11-12
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#
# Contents -----------------------------
# load_config()
# 
## Helpers -------
# .is_on()
# .num()
# .chr()
# .lag_name()
# .rate_name()
# parse_days()
# 
## Applying the parameter values -------
# get_fixed_params()
# apply_fixed_params()
# 
# measure_lag_pairs()
# build_possible_measures()
# 
# These ones are all for the desired states
# config_states()
# 
# resolve_summary()
# resolve_objective()
# 
# build_es_config()
# build_eval_settings()
# 
# build_deoptim_control()

# ---------------------------------------------------------------------------
# config_funs.R
# Turn input/optimisation_config.yml into the objects the optimisation uses.
# Needs: yaml, dplyr, tibble, DEoptim (and range_obj etc. from objective_funs
#        sourced before build_es_config() is called)
# ---------------------------------------------------------------------------

#' Title
#'
#' @param path path to the configuration yaml
#'
#' @returns list
#' @export
#'
#' @examples
load_config <- function(path) {
  if (!file.exists(path)) {
    stop("Config file not found: ", path)
  }
  
  cfg <- yaml::read_yaml(path)
  
  missing <- setdiff(c("run", "targets", "measures", "optimiser"), names(cfg))
  
  if (length(missing) > 0) {
    stop("Config is missing section(s): ", paste(missing, collapse = ", "))
  }
  
  cfg
}

# ---- small internals --------------------------------------------------------
.is_on <- function(x) is.null(x$enabled) || isTRUE(x$enabled)   # enabled unless set false
.num   <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
.chr   <- function(x) if (is.null(x)) NA_character_ else as.character(x)
.lag_name <- function(measure, lag) if (is.null(lag$name)) paste0(measure, "_lag") else lag$name
.rate_name <- function(measure, rate) if (is.null(rate$name)) paste0(measure, "_rate") else rate$name
# Value of a lag / rate block that is switched off but given a fixed `value` (NA if none)
.fixed_val <- function(blk) if (!is.null(blk) && !.is_on(blk) && !is.null(blk$value)) as.numeric(blk$value) else  NA_real_

parse_days <- function(x) {
  if (is.list(x) && all(c("start", "end") %in% names(x))) return(seq(x$start, x$end))
  if (is.numeric(x)) return(x)
  stop("eval_days must be {start: , end: } or a list of days.")
}

# ============================== MEASURES ======================================

# Parameters fixed by enabled measures (e.g. fManVeg needs cDayManVeg1). A parameter
# that is itself an enabled (optimised) measure takes precedence and is skipped. Used inside apply_fixed_params()
#' Title
#'
#' @param cfg optimisation configuration
#'
#' @returns
#' @export
#'
#' @examples
get_fixed_params <- function(cfg) {
  on <- purrr::keep(cfg$measures, .is_on)
  out <- list()
  for (nm in names(on)) {
    fp <- on[[nm]]$fixed_params
    for (p in names(fp)) {
      if (p %in% names(on)) {
        message("'", p, "' is optimised, so fixed_params of '", nm, "' is ignored for it.")
      } else {
        out[[p]] <- fp[[p]]
      }
    }
  }
  out
}

#' Title
#'
#' @param cfg optimisation configuration
#' @param lDATM the read in DATM file, as a list
#' @param set which DATM set to manipulate, Default sSet2
#'
#' @returns modified DATM list
#' @export
#'
#' @examples
apply_fixed_params <- function(cfg, lDATM, set = "sSet2") {
  fp <- get_fixed_params(cfg)
  for (p in names(fp)) {
    if (!p %in% rownames(lDATM$params)) {
      stop("Fixed parameter '", p, "' not found in lDATM$params.")
    }
    lDATM$params[p, set] <- fp[[p]]
    message("Fixed ", p, " = ", fp[[p]], " (", set, ")")
  }
  lDATM
}

# main parameter -> its lag and rate parameter (NA if not optimised), any fixed lag / rate
# value from the config, plus its current value. Includes enabled measures that have an
# optimised or fixed lag and/or rate.

#' Title
#'
#' @param cfg 
#'
#' @returns
#' @export
#'
#' @examples
measure_lag_pairs <- function(cfg) {
  on <- purrr::keep(cfg$measures, .is_on)
  purrr::imap(on, ~ {
    lg <- .x$lag
    rt <- .x$rate
    
    lg_on <- !is.null(lg) && .is_on(lg)
    rt_on <- !is.null(rt) && .is_on(rt)
    
    lg_fx <- .fixed_val(lg)
    rt_fx <- .fixed_val(rt)
    
    if (!lg_on && !rt_on && is.na(lg_fx) && is.na(rt_fx)) {
      return(NULL)
    }
    
    tibble::tibble(main = .y,
                   lag  = .lag_name(.y, lg))
  }) |> 
    list_rbind()
}




#' Title
#'
#' @param cfg optimisation config 
#'
#' @returns possible_measures dataframe
#' @export
#'
#' @examples

# Same data frame as the old possible_measures.csv + filter():
# parameter, upper_bound, lower_bound, current_val, note
# plus two columns describing each row for run_pathway():
#   role     'measure', 'lag' (when it starts, years) or 'rate' (years to reach the full value)
#   measure  the measure a lag / rate belongs to
build_possible_measures <- function(cfg) {
  on <- purrr::keep(cfg$measures, .is_on)
  
  # Apply checking for enabled measures
  if (length(on) == 0) {
    stop("No measures are enabled in the config.")
  }
  
  if (any(grepl("lag", names(on)))) {
    stop(
      "Measure names must not contain 'lag': ",
      paste(grep("lag", names(on), value = TRUE), collapse = ", ")
    )
  }
  
  if (any(grepl("_rate\\d*$", names(on)))) {
    stop(
      "Measure names must not end in '_rate': ",
      paste(grep("_rate\\d*$", names(on), value = TRUE), collapse = ", ")
    )
  }
  
  # Helper function
  mk <- function(parameter, x, current_val, note, role, measure) {
    tibble::tibble(parameter   = parameter,
                   upper_bound = .num(x$upper_bound),
                   lower_bound = .num(x$lower_bound),
                   current_val = current_val,
                   note        = note,
                   role        = role,
                   measure     = measure)
  }
  
  # identify the main variable
  main <- purrr::imap_dfr(on, function(x, nm) {
    mk(parameter   = nm,
       x           = x,
       current_val = .num(x$current_val),
       note        = .chr(x$note),
       role        = "measure",
       measure     = nm)
  })
  
  # lags of the variables
  lags <- purrr::imap_dfr(on, function(x, nm) {
    
    lg <- x$lag
    
    if (is.null(lg) || !.is_on(lg)) {
      return(NULL)
    }
    
    mk(parameter   = .lag_name(nm, lg),
       x           = lg,
       current_val = NA_real_,
       note        = lg$note %||% "must have main var",
       role        = "lag",
       measure     = nm)
  })
  
  # the rates of the variables
  rates <- purrr::imap_dfr(on, function(x, nm) {
    
    rt <- x$rate
    
    if (is.null(rt) || !.is_on(rt)) {
      return(NULL)
    }
    
    mk(parameter   = .rate_name(nm, rt),
       x           = rt,
       current_val = NA_real_,
       note = rt$note %||%
         "years to reach the full value once the lag has passed",
       role        = "rate",
       measure     = nm)
  })
  
  # all measures first, then all lags, then all rates
  out <- dplyr::bind_rows(main, lags, rates)
  
  # ---- validation ----
  # Checks requirements and the names for duplication etc. 
  if (anyNA(out$upper_bound) || anyNA(out$lower_bound)) {
    stop("Every enabled measure/lag/rate needs lower_bound and upper_bound: ", 
         paste(out$parameter[is.na(out$upper_bound) | is.na(out$lower_bound)], collapse = ", "))
  }
  
  if (any(out$lower_bound > out$upper_bound)) {
    stop("lower_bound > upper_bound for: ", paste(out$parameter[out$lower_bound > out$upper_bound], collapse = ", "))
  }
  
  if (any(out$role != "measure" & out$lower_bound < 0)) {
    stop("Lags and rates cannot be negative: ",
         paste(out$parameter[out$role != "measure" & out$lower_bound < 0], collapse = ", "))
  }
  
  if (anyDuplicated(out$parameter)) {
    stop("Duplicate parameter names: ", paste(out$parameter[duplicated(out$parameter)],collapse = ", "))
  }
  
  bad_lag <- out$parameter[out$role == "lag" & !grepl("lag", out$parameter)]
  
  bad_rate <- out$parameter[out$role == "rate" & !grepl("_rate", out$parameter)]
  
  if (length(bad_lag) > 0) {
    stop("Lag parameter names must contain 'lag': ",paste(bad_lag, collapse = ", "))
  }
  
  if (length(bad_rate) > 0) {
    stop("Rate parameter names must contain '_rate': ",paste(bad_rate, collapse = ", "))
  }
  
  for (nm in names(on)) {
    miss <- setdiff(unlist(on[[nm]]$requires),names(on))
    
    if (length(miss) > 0) {
      stop("Measure '", nm, "' requires ", paste(miss, collapse = ", ")," to be enabled as well.")
    }
  }
  
  # lags / rates that are NOT optimised but fixed at a value in the config
  # (e.g.  rate: {enabled: false, value: 5}). Kept as an attribute so the object
  # that is passed to run_pathway() as measure_map carries them.
  
  fixed <- purrr::imap_dfr(on, function(x, nm) {
    
    purrr::map_dfr(c("lag", "rate"), function(r) {
      
      blk <- x[[r]]
      
      if (is.null(blk) || is.null(blk$value)) {
        return(NULL)
      }
      
      if (.is_on(blk)) {
        return(NULL)
      }
      
      message("Using fixed ", r, " for '", nm,"': ", blk$value)
      
      tibble::tibble(measure = nm,
                     role = r,
                     value = .num(blk$value))
    })
  })
  
  if (nrow(fixed) > 0) {
    if (anyNA(fixed$value) || any(fixed$value < 0)) {
    
    bad <- is.na(fixed$value) | fixed$value < 0
    
    stop("Fixed lag / rate values must be numbers >= 0: ",
         paste(paste(fixed$measure, fixed$role)[bad],collapse = ", "))
  }
  
  attr(out, "fixed_timing") <- fixed
  }
  
  
  out
}

# ============================== TARGETS =======================================

# States used by enabled ESs/indicators (no function lookup needed - safe in plotting sessions)
#' Title
#'
#' @param cfg optimisation config
#'
#' @returns the unique states/indicators listed under the ES
#' @export
#'
#' @examples
config_states <- function(cfg) {
  cfg$targets$ecosystem_services |>
    purrr::keep(.is_on) |>
    purrr::map(
      ~ .x$indicators |>
        purrr::keep(.is_on) |>
        names()
    ) |>
    unlist(use.names = FALSE) |>
    unique()
  
}

#' Title
#'
#' @param name 
#'
#' @returns function
#' @export
#'
#' @examples
#' Convert a text string from YAML into an actual R function.
resolve_summary <- function(name) {
  
  if (is.null(name)) {
    stop("Every indicator needs a 'summary' (max, min, mean, median or pNN).")
  }
  
  if (grepl("^p[0-9]+$", name)) {
    q <- as.numeric(sub("^p", "", name)) / 100
    return(function(x) stats::quantile(x, probs = q, names = FALSE))
  }
  
  switch(name, max = max, min = min, mean = mean, median = stats::median,
         stop("Unknown summary '", name, "'. Use max, min, mean, median or pNN (e.g. p90)."))
}

#' Title
#'
#' @param name 
#'
#' @returns list list(f = below_obj, n = 1L)
#' @export
#'
#' @examples
#' Convert an objective name in the config into: a) the objective function (range, below, etc.); b)the number of target values required
resolve_objective <- function(name) {
  # looked up when called, so optim_functions.R must already be sourced
  objs <- list(range = list(f = range_obj, n = 2L), below = list(f = below_obj, n = 1L),
               above = list(f = above_obj, n = 1L), exact = list(f = exact_obj, n = 1L))
  
  if (is.null(name) || !name %in% names(objs)) {
    stop("Unknown objective '", name, "'. Use one of: ", paste(names(objs), collapse = ", "))
  }
  
  objs[[name]]
}

# Config list in the structure evaluate_es() expects (functions resolved from names)
#' Title
#'
#' @param cfg optimisation config
#'
#' @returns
#' @export
#'
#' @examples
#' Main constructor of the desired ESs. Converts a YAMLconfig into the structure expected by evaluate_pathway_es().
build_es_config <- function(cfg) {
  es_cfg <- purrr::keep(cfg$targets$ecosystem_services, .is_on)
  
  if (length(es_cfg) == 0) {
    stop("No ecosystem services enabled in the config.")
  }
  
  purrr::map(es_cfg, function(es) {
    inds <- purrr::keep(es$indicators, .is_on)
    
    if (length(inds) == 0) {
      stop("An enabled ecosystem service has no enabled indicators.")
    }
    
    ind_list <- purrr::imap(inds, function(ind, st) {
      obj <- resolve_objective(ind$objective)
      tg  <- as.numeric(ind$target)
      
      if (length(tg) != obj$n) {
        stop("Indicator '", st,"': objective '", ind$objective,
             "' needs ", obj$n, " target value(s).")
      }
      
      out <- list(target  = tg,
                  summary = resolve_summary(ind$summary),
                  obj     = obj$f)
      
      for (f in c("weight", "scale", "met_tol")) {
        if (!is.null(ind[[f]])) {
          out[[f]] <- ind[[f]]
        }
      }
      
      if (!is.null(ind$eval_days)) {
        out$eval_days <- parse_days(ind$eval_days)
      }
      
      out
    })
    
    out <- list(aggregate = es$aggregate %||% "mean",
                indicators = ind_list)
    
    if (!is.null(es$weight)) {
      out$weight <- es$weight
    }
    
    out
  })
}

# Arguments for evaluate_es(): eval_days, eval_year, nonneg, bound
#' Title
#'
#' @param cfg optimisation config
#'
#' @returns list, generic evaluation settings
#' @export
#'
#' @examples
build_eval_settings <- function(cfg) {
  e <- cfg$targets$evaluation
  list(eval_days = parse_days(e$eval_days),
       eval_year = if (is.null(e$eval_year)) "max" else e$eval_year,
       nonneg    = if (is.null(e$nonneg)) TRUE else isTRUE(e$nonneg),
       bound     = if (is.null(e$bound)) FALSE else isTRUE(e$bound))
}


# ============================== OPTIMISER =====================================

build_deoptim_control <- function(cfg, n_measures) {
  optim_config <- cfg$optimiser
  ctrl <- optim_config[setdiff(names(optim_config), "NP_per_measure")]
  bad  <- setdiff(names(ctrl), names(formals(DEoptim::DEoptim.control))) # formals identifies the arguments needed
  
  if (length(bad) > 0) {
    stop("Unknown DEoptim option(s) in config: ", paste(bad, collapse = ", "))
  }
  
  if (is.null(ctrl$NP)) {
    if (is.null(optim_config$NP_per_measure)) stop("Set optimiser$NP_per_measure (or NP).")
    ctrl$NP <- optim_config$NP_per_measure * n_measures
  }
  ctrl$parallelType <- "foreach"     # the script sets up the foreach cluster itself
  ctrl
}

