#--------------------------------------#
## Project: Pathway optimisation - Lake District
## Script purpose: functions for reporting and post-processing
## Date: 2026-11-12
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#
# Contents -----------------------------
# report_pathway()
# log_to_long()
# summarise_achievement()

#' report_pathway outputs
#'
#' @param val_pars 
#' @param name_pars 
#' @param obj_val 
#' @param log_dir 
#' @param extra optional extra columns in the log (1-row data.frame)
#'
#' @returns
#' @export
#'
#' @examples
report_pathway  <- function(val_pars, 
                            name_pars, 
                            obj_val,
                            log_dir = log_dir, extra = NULL) {
  
  
  # Log the results of the objective_function
  fname <- file.path(log_dir, paste0("eval_", 
                                     format(Sys.time(), "%Y%m%d%H%M%OS6"), "_",
                                     Sys.getpid(), "_", paste(sample(letters, 4), collapse=""), ".csv"))
  report <- data.frame(matrix(val_pars,1,length(name_pars)))
  colnames(report) <- name_pars
  
  out <- cbind.data.frame(obj = obj_val, report)
  
  if (!is.null(extra)) {
    out <- cbind.data.frame(out, extra)
  }
  
  readr::write_csv(out, fname)
  
  # could also log the things that are being evaluated, parameter values
  
  
}

# ---------------------------------------------------------------------------#
# Post-processing: turn the wide log back into tidy tables
#   ind_long <- log_to_long(log_df, "ind")  # one row per member x ES x state
#   es_long  <- log_to_long(log_df, "es")   # one row per member x ES
# ---------------------------------------------------------------------------#
log_to_long <- function(log_df, level = c("ind", "es")) {
  level <- match.arg(level)
  cols  <- grep(paste0("^", level, "__"), names(log_df), value = TRUE)
  parts <- strsplit(cols, "__", fixed = TRUE)
  ids   <- intersect("iteration", names(log_df))
  
  long <- do.call(rbind, lapply(seq_along(cols), function(k) {
    p <- parts[[k]]
    data.frame(log_df[ids], member = seq_len(nrow(log_df)),
               es = p[2], state = if (level == "ind") p[3] else NA_character_,
               field = p[length(p)], value = as.character(log_df[[cols[k]]]),
               stringsAsFactors = FALSE)
  }))
  wide <- tidyr::pivot_wider(long, names_from = field, values_from = value)
  if (level == "es") wide$state <- NULL
  
  keep_chr <- c(ids, "member", "es", "state", "status", "limiting_indicator")
  dplyr::mutate(wide, dplyr::across(-dplyr::any_of(keep_chr), ~ utils::type.convert(.x, as.is = TRUE)))
}

# Population-level answer to "which targets are achieved, and why not?" --------
summarise_achievement <- function(log_df) {
  ind <- log_to_long(log_df, "ind")
  es  <- log_to_long(log_df, "es")
  
  list(
    es = es |>
      dplyr::group_by(es) |>
      dplyr::summarise(pct_members_met = 100 * mean(es_met),
                       mean_prop_indicators_met = mean(prop_indicators_met),
                       most_common_limiting = {
                         l <- limiting_indicator[!is.na(limiting_indicator)]
                         if (length(l)) names(sort(table(l), decreasing = TRUE))[1] else NA_character_
                       },
                       .groups = "drop"),
    indicators = ind |>
      dplyr::group_by(es, state) |>
      dplyr::summarise(pct_met = 100 * mean(met),
                       pct_too_low  = 100 * mean(status == "too_low"),
                       pct_too_high = 100 * mean(status == "too_high"),
                       median_gap = stats::median(gap),
                       median_prop_days_met = stats::median(prop_days_met),
                       .groups = "drop")
  )
}


