#--------------------------------------#
## Project: Pathway optimisation - Lake District
## Script purpose: functions for evaluating pathways and interpretting model outputs
## Date: 2026-11-12
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#
# Contents------------------------------
# es_states()
# validate_es_config()
# weight_audit()
# validate_eval_settings()
# evaluate_pathway_es()
# flatten_es_result()
# evaluate_pathway()   # legacy version, copied from optim_functions 12/11/2026


# ---------------------------------------------------------------------------
# Ecosystem-service (ES) level evaluation
# Replaces state-level evaluate_pathway(). The ES / indicator definitions (states, targets,
# summaries, weights, aggregation) are NOT set here: they come from
# input/optimisation_config.yml via build_es_config() in R/config_funs.R.
# An indicator = (ES, state) pair, so one state can serve several ESs.
# Weights are optional and re-normalised (equal by default) so the number of
# indicators in an ES never changes how much that ES counts.
# ---------------------------------------------------------------------------
#' Title
#'
#' @param es_config the config for ES, from build_es_config
#'
#' @returns
#' @export
#'
#' @examples
es_states <- function(es_config) {
  unique(unlist(lapply(es_config, function(e) names(e$indicators))))
}

# Fill defaults, validate, and normalise weights #
#' Title
#'
#' @param es_config the config for ES, from build_es_config
#'
#' @returns
#' @export
#'
#' @examples
#' validates names, validates required fields (target, summary, obj), normalises
#' ES and indicator weights, adds in missing defaults
validate_es_config <- function(es_config) {
  
  # Define helper function
  chk_names <- function(x) {
    if (is.null(names(x)) || any(names(x) == "") || anyDuplicated(names(x)) ||
        any(make.names(names(x)) != names(x)) || any(grepl("__", names(x), fixed = TRUE))) {
      stop("ES / state names must be unique, syntactic R names without '__': ",
           paste(names(x), collapse = ", "))
    }
  }
  
  norm_w <- function(w, what) {
    if (all(is.na(w))) w[] <- 1
    if (any(is.na(w))) stop("Give a weight to every ", what, " or to none.")
    if (any(w < 0) || sum(w) <= 0) stop("Weights for ", what, " must be non-negative and not all zero.")
    w / sum(w)
  }
  
  get_w <- function(x) vapply(x, function(i) if (is.null(i$weight)) NA_real_ else i$weight, numeric(1))
  #---------------
  
  chk_names(es_config)
  es_w <- norm_w(get_w(es_config), "ES")
  
  for (es in names(es_config)) {
    ind <- es_config[[es]]$indicators
    if (length(ind) == 0) stop("ES '", es, "' has no indicators.")
    chk_names(ind)
    w <- norm_w(get_w(ind), paste0("indicator in ", es))
    for (st in names(ind)) {
      i <- ind[[st]]
      if (!is.function(i$summary) || !is.function(i$obj) || is.null(i$target)) {
        stop("Indicator ", es, "/", st, " needs target, summary (function) and obj (function).")
      }
      ind[[st]]$weight  <- w[[st]]
      if (is.null(i$scale))   ind[[st]]$scale   <- 1
      if (is.null(i$met_tol)) ind[[st]]$met_tol <- 1e-9
    }
    es_config[[es]]$indicators <- ind
    es_config[[es]]$weight     <- es_w[[es]]
    if (is.null(es_config[[es]]$aggregate)) es_config[[es]]$aggregate <- "mean"
  }
  es_config
}

# Audit the weighting: how much does each indicator / state really count? 
#' Title
#'
#' @param es_config the config for ES, from build_es_config
#'
#' @returns data.frame of weights
#' @export
#'
#' @examples
weight_audit <- function(es_config) {
  d <- validate_es_config(es_config) |>
    purrr::imap_dfr(~ tibble::tibble(es = .y,
                                     state = names(.x$indicators), 
                                     es_weight = .x$weight,       
                                     ind_weight = purrr::map_dbl(.x$indicators, "weight"))) |>
    dplyr::mutate(eff_weight = es_weight * ind_weight)
  
  d
  # sum(d$eff_weight) == 1; aggregate(eff_weight ~ state, d, sum) = total pull of each state
}


# Validate the evaluation settings and fill in defaults. Unknown keys are an error,
# so a typo in targets.evaluation in the config is caught rather than silently ignored.
#' Title
#'
#' @param eval_settings list from build_eval_settings(config)
#'
#' @returns modified eval_settings list
#' @export
#'
#' @examples
validate_eval_settings <- function(eval_settings = list()) {
  defaults <- list(eval_days = 120:275, eval_year = "max", nonneg = TRUE, bound = FALSE)
  unknown <- setdiff(names(eval_settings), names(defaults))
  if (length(unknown) > 0) {
    stop("Unknown evaluation setting(s): ", paste(unknown, collapse = ", "),
         ". Allowed: ", paste(names(defaults), collapse = ", "))
  }
  purrr::list_modify(defaults, !!!eval_settings)
}


#' Evaluate a pathway at ES level
#'
#' @param PCLake_output output from run_pathway
#' @param es_config     see structure above
#' @param eval_settings named list (from build_eval_settings(cfg)); any missing element
#'   takes its default:
#'   \itemize{
#'     \item eval_days: days of year evaluated (an indicator can override), default 120:275
#'     \item eval_year: 'max' (last simulated year) or a year number, default 'max'
#'     \item nonneg: floor indicator errors at 0, so beating one target cannot buy back
#'       a miss elsewhere (and obj == 0 means "all targets met"), default TRUE
#'     \item bound: if TRUE, squash errors to [0,1) via e/(1+e) so no single indicator
#'       can dominate through sheer magnitude (e.g. a huge chlorophyll value). Keeps a
#'       gradient for DEoptim, unlike capping. Default FALSE
#'   }
#' @param return_detail FALSE -> scalar for DEoptim; TRUE -> list with total, es table,
#'                      indicator table, n_es_met, all_met
evaluate_pathway_es <- function(PCLake_output,
                                es_config,
                                eval_settings = list(),
                                return_detail = FALSE) {
  
  eval_settings <- validate_eval_settings(eval_settings)
  
  eval_days <- eval_settings$eval_days
  eval_year <- eval_settings$eval_year
  nonneg <- eval_settings$nonneg
  bound <- eval_settings$bound
  
  es_config <- validate_es_config(es_config)
  
  missing_states <- setdiff(es_states(es_config),
                            colnames(PCLake_output))
  
  if (length(missing_states) > 0) {
    stop("PCLake output doesn't contain: ",paste(missing_states, collapse = ", "))
  }
  
  # same year / day-of-year logic as evaluate_pathway()
  out <- PCLake_output
  
  out$year <- floor((out$time - 1) / 365) + 1
  
  out$doy <- lubridate::yday(lubridate::as_date(out$time - (out$year * 365) + 364,
                                                origin = "2025-01-01"))
  
  if (identical(eval_year, "max")) {
    eval_year <- max(out$year)
  }
  
  out <- dplyr::filter(out, year == eval_year)
  
  res <- purrr::imap(es_config, function(cfg, es) {
    
    # ---- indicator level: what happened, and was the target met? --------------
    d <- purrr::imap_dfr(cfg$indicators, function(ind, st) {
      
      days <- ind$eval_days %||% eval_days #"if NULL, use the right-hand side".
      x <- out[[st]][out$doy %in% days]
      tg <- ind$target
      
      value <- ind$summary(x) # calls the summary function as specified in the config (e.g. mean of the period)
      
      raw <- ind$obj(value, tg) # compares with the target and calls the eval function (range_obj, etc.); <= met_tol means target met
      met <- raw <= ind$met_tol
      
      err <- raw / ind$scale
      
      if (nonneg) {
        err <- max(0, err)
      }
      
      if (bound && err > 0) {
        err <- err / (1 + err)
      }
      
      # signed gap to the target in STATE units (0 = inside range)
      gap <- if (length(tg) == 2) {
        value - min(max(value, min(tg)), max(tg))
      } else {
        value - tg
      }
      
      tibble::tibble(es = es,
                     state = st,
                     summary_value = value,
                     target_lo = tg[1],
                     target_hi = if (length(tg) == 2) tg[2] else NA_real_,
                     gap = gap,
                     met = met,
                     status = ifelse(met, "met", ifelse(gap > 0, "too_low", "too_high")),
                     # how often does the state itself sit within target across the eval window?
                     prop_days_met = mean(purrr::map_lgl(x,~ ind$obj(.x, tg) <= ind$met_tol)),
                     error = err,
                     ind_weight = ind$weight)
    })
    
    # ---- ES level --------------------------------------------------------------
    score <- switch(cfg$aggregate,
                    mean = sum(d$error * d$ind_weight),
                    max = max(d$error),
                    stop("aggregate must be 'mean' or 'max'"))
    
    lim <- if (all(d$met)) {
      NA_character_
    } else {
      d$state[which.max(d$error)]
    }
    
    d <- dplyr::mutate(d,
                       es_weight = cfg$weight,
                       eff_weight = cfg$weight * ind_weight,
                       contribution = if (cfg$aggregate == "max") {
                         cfg$weight * ifelse(state == lim, error, 0)
                       } else { eff_weight * error})
    
    list(detail = d,
         es = tibble::tibble(es = es,
                             es_weight = cfg$weight,
                             es_score = score,
                             es_met = all(d$met), # strict: every indicator met
                             prop_indicators_met = sum(d$ind_weight * d$met), # partial credit
                             limiting_indicator = lim)) # biggest reason the ES fails
  })
  
  detail <- dplyr::bind_rows(purrr::map(res, "detail"))
  es_tab <- dplyr::bind_rows(purrr::map(res, "es"))
  total <- sum(es_tab$es_weight * es_tab$es_score)
  
  es_tab <- dplyr::mutate(es_tab, 
                          share_of_total = ifelse(total > 0,
                                                  es_weight * es_score / total,
                                                  0))
  
  detail <- dplyr::mutate(detail,
                          share_of_total = ifelse(total > 0,
                                                  contribution / total,
                                                  0))
  
  if (!return_detail) {
    return(total)
  }
  
  list(total = total,
       es = es_tab,
       detail = detail,
       n_es_met = sum(es_tab$es_met),
       all_met = all(es_tab$es_met))
}


# ---------------------------------------------------------------------------
# One-row wide table of everything above, for the per-member log
# columns: ind__<ES>__<state>__<field>, es__<ES>__<field>, n_es_met, all_met
# ---------------------------------------------------------------------------
flatten_es_result <- function(res) {
  d <- res$detail
  e <- res$es
  
  ind_fields <- c("summary_value", "gap", "status", "met", "prop_days_met",
                  "error", "eff_weight", "share_of_total")
  es_fields  <- c("es_score", "es_met", "prop_indicators_met", "limiting_indicator", "share_of_total")
  
  out <- list()
  
  for (i in seq_len(nrow(d))) for (f in ind_fields) {
    out[[paste("ind", d$es[i], d$state[i], f, sep = "_")]] <- d[[f]][i]
  }
  
  for (j in seq_len(nrow(e))) for (f in es_fields) {
    out[[paste("es", e$es[j], f, sep = "_")]] <- e[[f]][j]
  }
  
  out$n_es_met <- res$n_es_met
  out$all_met  <- res$all_met
  do.call(data.frame, c(out, list(check.names = FALSE, stringsAsFactors = FALSE)))
}

#' Evaluate a parameter pathway - states only (deprecated!!)
#'
#' @param PCLake_output output from run_pathway
#' @param future_states dataframe of variable and target value, with weights
#' @param eval_target how should the target be evaluated, list of functions to call, matched to the state names
#' @param eval_days which days in the last year should be evaluated
#'
#' @returns numeric value output from eval_target function call
#' @export
#'
#' @examples
evaluate_pathway <- function(PCLake_output, 
                             future_states,
                             eval_target = list(function(out,target){abs(out-target)/target}),
                             eval_days = 121:244,
                             eval_funs, # need to provide a named function
                             eval_year = 'max') {
  
  # For debugging ----------------- #
  # PCLake_output = model_output
  # future_states = desired_states
  # eval_target = list(oChlaEpi = range_obj)
  # eval_days = 121:244
  # #--------------------------------#
  
  if (sum(!names(future_states) %in% colnames(PCLake_output)) > 0) {
    stop("PCLake output doesn't contain the evaluation variables from the future states.")
  }
  
  if (sum(!names(eval_target) %in% colnames(PCLake_output)) > 0) {
    stop("PCLake output doesn't contain the evaluation variables from the eval_target. Either remove from the obj_function or amend future_states.")
  }
  
  if(sum(!names(future_states) %in% names(eval_target)) > 0){
    stop("The list of evaluation functions doesn't include all the desired states.")
  }
  
  if (sum(unlist(lapply(future_states, get('weights')))) != 1) {
    stop("Weights should add to 1.")
  }
  
  if(sum(!names(eval_target) %in% names(eval_funs)) > 0){
    stop("The list of evaluation functions doesn't match the evaluation.")
  }  
  
  # Extract model output and compare with desired state
  if(eval_year == 'max'){
    eval_year <- PCLake_output |> 
      mutate(year = floor((time-1)/365) + 1) |> 
      summarise(max = max(year)) |> 
      pull()
  } 
  
  model_output <- PCLake_output |>
    mutate(year = floor((time-1)/365) + 1,
           doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
    filter(year == eval_year, # filters to summer in the last year of the simulation
           doy %in% eval_days) |> 
    select(names(future_states)) |> 
    summarise(across(any_of(names(future_states)),~ eval_funs[[cur_column()]](.x))) |>
    pivot_longer(cols = any_of(names(future_states)),
                 names_to = 'variable',
                 values_to = 'output')
  
  # need a more complex process if the evaluation happens seperately (different function by objective)
  # Check for matching names across the eval_functions and future states
  if (length(names(eval_target)) > 0) { # if the functions are named check them
    
    if (sum(stringr::str_equal(names(eval_target), names(future_states))) != length(future_states)) {
      stop('the name(s) of your eval_function(s) dont match the future states') 
    } else {
      
      message('matching evaluation function by state')
      
      pathway_error <- model_output |> 
        mutate(diff= NA,
               diff_weighted = NA)
      
      for (i in 1:nrow(pathway_error)) {
        use_fun <- eval_target[[which(names(eval_target) == unlist(pathway_error['variable'][i,]))]]
        pathway_error$diff[i] <- use_fun(pathway_error$output[i], future_states[[which(names(future_states) == pathway_error$variable[i])]]$target)
        
        pathway_error$diff_weighted[i] <- pathway_error$diff[i] * future_states[[which(names(future_states) == pathway_error$variable[i])]]$weights
        
      }
      
      pathway_error <- pathway_error  |> 
        summarise(total_error = sum(diff_weighted)) ### IS THIS HOW YOU WOULD SUM THEM????
      
    }
    
  } else {
    if (length(eval_target) == 1) { # if there is only one function 
      
      if (length(future_states) > 1) { # but multiple states give a warning message
        message('Warning: using a single evaluation function')
      }
      
      pathway_error <- model_output |> 
        mutate(diff= NA,
               diff_weighted = NA)
      
      for (i in 1:nrow(model_output)) {
        pathway_error$diff[i] <- eval_target[[1]](pathway_error$output[i],
                                                  future_states[[which(names(future_states) == pathway_error$variable[i])]]$target)
        pathway_error$diff_weighted[i] <- pathway_error$diff[i] * future_states[[which(names(future_states) == pathway_error$variable[i])]]$weights
        
      }
      
      pathway_error <- pathway_error |> 
        summarise(total_error = sum(diff_weighted)) 
      
    } else {
      stop("If you are passing more than one function they need to be named") # cannot pass more than one unnamed function
    }
  }
  
  
  out_val <- pathway_error |> pull(total_error)
  
  return(out_val)
  
}
