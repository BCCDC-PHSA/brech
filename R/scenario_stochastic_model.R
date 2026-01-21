
###################################
##                               ##
##  Model scenarios section      ##
##                               ##
###################################

#' Create a set of scenarios that can
#' be incorporated into [projection_stochastic_model()]
#' @param sm Object of class \code{stochastic_model}.
#' @param parameters A data frame or tibble of prior parameter values. Each row
#' should correspond to one simulation, and columns should match the parameter
#' names expected by `update_parameters()`. If `NULL` then no parameters will
#' be updated an are all inherited from `sm`
#' @param states A data frame where each row represents an initial state. If
#' `NULL` then this is generated from the initial states in `sm`
#' @seealso [abc_stochastic_model()]
#' @export
scenario_stochastic_model <- function(sm,parameters=NULL,states=NULL){
  if(is.null(parameters) && is.null(states)){
    stop("Need to define at least one of `parameters` or `states`")
  }
  if(!is.null(parameters) && !is.null(states)){
    if(nrow(parameters) != nrow(states)){
      stop("`parameters` and `states` need the same number of rows")
    }
  }
  nreps <- nrow(if (!is.null(parameters)) parameters else states)
  if(is.null(states)){
    states <- as.data.frame(sm$initial_states)
    colnames(states) <- names(sm$initial_states)
    states <- states[rep(row.names(states), times = nreps), ]
  }
  if(is.null(parameters)){
    parameters <- as.data.frame(sm$params)
    colnames(parameters) <- names(sm$params)
    parameters <- parameters[rep(row.names(parameters), times = nreps), ]
  }

  structure(list("model" = sm, "parameters" = parameters, "states" = states),
            class = "scenario_stochastic_model")
}
