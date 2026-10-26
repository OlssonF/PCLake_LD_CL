#' Run PCLake management pathway with modified model parameters
#'
#' @param val_pars values of the pathway parameters
#' @param name_pars names of the pathway parameters
#' @param initial_conditions dataframe of initial values (from prepInitials() or use default)
#'
#' @returns a dataframe of PCLake output as defined in the DATM file
#' @export
#'
#' @examples

run_pathway <- function(val_pars, name_pars, current_val, initial_conditions = NULL) {
  
  # For debugging ----------------- #
  # val_pars <- possible_measures$lower_bound
  # name_pars <- possible_measures$parameter
  # current_val <- possible_measures$current_val
  # initial_conditions <- equilibrium_states
  #--------------------------------#
  
  # Setting parameter values -----------#
  lDATM_SETTINGS_obj <- lDATM_SETTINGS
  lag_pars_lag <- grep('_lag', name_pars) # these are model parameters with lags - take these out
  lag_pars <- which(name_pars %in% gsub('_lag', '', name_pars[lag_pars_lag]))
  
  # it wasn't working if there were no lagged pars
  if (length(lag_pars) == 0) {
    none_lag_pars <- name_pars
  } else {
    none_lag_pars <- name_pars[-c(lag_pars, lag_pars_lag)]
  }
  
  # check if these non-lagged parameter values are listed in the forcings and change there if it is, if not just change the parameter value
  for (p in none_lag_pars) {
    if (p %in% names(lDATM_SETTINGS_obj$forcings$sSet2)) {
      lDATM_SETTINGS_obj$forcings$sSet2[[p]]$value <- val_pars[which(name_pars ==p)]
    } else {
      lDATM_SETTINGS_obj$params[name_pars[-c(lag_pars, lag_pars_lag)], "sSet2"] <- val_pars[-c(lag_pars, lag_pars_lag)]
    }
  }
  
  
  # set the lagged variables via the forcings in the DATM file
  # extract lagged variable 
  # NOTE: this  only works if there is already something in the forcings and you are just modifying it
  # otherwise you have to recompile the model
  
  unique_lags <- unique(gsub("^\\d+|\\d+$", "", name_pars[lag_pars_lag])) # could have multiple lag values e.g. mPLoadEpi_lag2, so identify the unqiue unnumbered measures
  unique_names <- unique(gsub("^\\d+|\\d+$", "", name_pars[lag_pars])) 
  
  for (i in 1:length(unique_lags)) {
    
    # extract the names
    name1 <- name_pars[str_detect(name_pars, unique_lags[i])]
    name2 <- name_pars[str_detect(name_pars, paste0(unique_names[i], '[^_]', '|',
                                                    unique_names[i], '$'))]
    
    # what are the possible values (current values plus those specified in the pathway parameters)
    values_use <- c(unique(current_val[which(name_pars %in% name2)]), 
                    val_pars[which(name_pars == name2)])
    
    # what are the time points that things change? (the start t = 0, plus any lag values to implement measures)
    lags_use <- c(0, val_pars[which(name_pars %in% name1)])[order(c(0, val_pars[which(name_pars %in% name1)]))]
    
    df <- data.frame(l = floor(lags_use*365), # round to an integer
                     value = values_use[order(c(0, val_pars[which(name_pars %in% name1)]))]) # make sure the lag is matched with its associated parameter
    
    
    forcing_df <- data.frame(time = 0:(lDATM_SETTINGS_obj$run_settings["dReady", "Set2"]*365)) |> 
      full_join(df, by = join_by(time == l)) |> 
      fill(value, .direction = 'down')
    
    if (nrow(lDATM_SETTINGS_obj$forcings$sSet2[[unique_names[i]]]) != nrow(forcing_df)) {
      
      lDATM_SETTINGS_obj$forcings$sSet2[[unique_names[i]]] <- NULL # reset the df
    }
    
    lDATM_SETTINGS_obj$forcings$sSet2[[unique_names[i]]] <- forcing_df
    
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


#' Evaluate a parameter pathway
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


#' report_pathway outputs
#'
#' @param val_pars 
#' @param name_pars 
#' @param obj_val 
#' @param log_dir 
#'
#' @returns
#' @export
#'
#' @examples
report_pathway  <- function(val_pars, 
                            name_pars, 
                            obj_val,
                            log_dir = log_dir) {
  
  
  # Log the results of the objective_function
  fname <- file.path(log_dir, paste0("eval_", 
                                     format(Sys.time(), "%Y%m%d%H%M%OS6"), "_",
                                     Sys.getpid(), "_", paste(sample(letters, 4), collapse=""), ".csv"))
  report <- data.frame(matrix(val_pars,1,length(name_pars)))
  colnames(report) <- name_pars
  
  write_csv(cbind.data.frame(obj = obj_val,report), fname)
  
  # could also log the things that are being evaluated, parameter values
  
  
}

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


# function(out, lower, upper) {
#   range_vals <- range(c(lower,upper,out))
#   fun_val <- (lower - range_vals[1] + range_vals[2] - b)/mean(c(lower,upper))
#   return(fun_val)
# }


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
