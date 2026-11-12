#--------------------------------------#
## Project: Pathway optimisation - Lake District
## Script purpose: functions for evaluation of pathways used in the objective functions
## Date: 2026-11-12
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

#' range_obj
#'
#' @param out value being evaluated
#' @param target vector of values n = 2, lower and upper
#'
#' @returns value minimising to zero
#' @export
#'
#' @examples
range_obj <- function(out, target) {
  if (length(target) !=2) {
    stop('Can only take 2 target values')
  } else {
    
    range_vals <- range(c(target,out))
    fun_val <- (min(target) - range_vals[1] + range_vals[2] - max(target))/mean(target)
    return(fun_val)
  }
}


below_obj <- function(out,target){ # below better
  if (length(target) !=1) {
    stop('Can only take 1 target value')
  } else {
    (out-target)/target
  }
  
}  

exact_obj <- function(out,target){#,  # target exact value
  if (length(target) !=1) {
    stop('Can only take 1 target value')
  } else { 
    abs(out-target)/target
  }
}

above_obj <- function(out,target){ # above better
  if (length(target) !=1) {
    stop('Can only take 1 target value')
  } else {
    (target-out)/target
  }     
}

convert_Pload <- function(before = 0.002, after) {
  
  change <- 100*(before - after)/before
  
  return(change)
}