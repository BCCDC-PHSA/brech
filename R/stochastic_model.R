##################################################
##################################################
##                                              ##
##        stochastic_model functions            ##
##                                              ##
##################################################
##################################################

#' convert reactions list to transitions and rates
#' @param reactions list
#' @return list of rates and transitions
get_reactions_and_rates <- function(reactions){
  # General function to get vector of rates
  rates <- function(x, p, t) {
    sapply(reactions, function(rxn) rxn$rate(x, p, t))
  }
  transitions <- lapply(reactions, `[[`, "transition")
  return(list("rates"=rates,"transitions"=transitions))
}

#' constructor for `stochastic_model` class
#' @noRd
new_stochastic_model <- function(reactions,sim_scenario){
  rr <- get_reactions_and_rates(reactions)
  rr <- c(rr,sim_scenario)
  structure(rr, class="stochastic_model")
}

#' validator for stochastic model
#' @noRd
validate_reactions <- function(reactions) {
  if (!is.list(reactions)) stop("`reactions` must be a list.")

  for (name in names(reactions)) {
    reaction <- reactions[[name]]

    # Check if each reaction is a list
    if (!is.list(reaction)) stop(sprintf("Reaction '%s' must be a list.", name))

    # Check that 'transition' exists and is a named numeric vector
    if (!("transition" %in% names(reaction))) {
      stop(sprintf("Reaction '%s' must have a 'transition' element.", name))
    }
    transition <- reaction$transition
    if (!is.numeric(transition) || is.null(names(transition))) {
      stop(sprintf("Transition in reaction '%s' must be a named numeric vector.", name))
    }

    # Check that 'rate' exists and is a function
    if (!("rate" %in% names(reaction))) {
      stop(sprintf("Reaction '%s' must have a 'rate' function.", name))
    }
    rate_fn <- reaction$rate
    if (!is.function(rate_fn)) {
      stop(sprintf("Rate in reaction '%s' must be a function.", name))
    }

    # Check that the rate function has at least 3 arguments: x, p, t
    rate_formals <- names(formals(rate_fn))
    if (length(rate_formals) < 3 || !all(c("x", "p", "t") %in% rate_formals[1:3])) {
      stop(sprintf("Rate function in reaction '%s' must have arguments (x, p, t).", name))
    }
  }
  return(TRUE)
}

#' @noRd
validate_sim_scenario <- function(sim_scenario){
  for(name in c("initial_states","params")){
    if(!(name %in% names(sim_scenario))){
      stop("sim_scenario must have '%s' list", name)
    }
  }
  # check that `params` is a list
  if(typeof(sim_scenario$params)!= "list"){
    stop("`params` must be a list")
  }

  return(TRUE)
}

#' Create a stochastic model
#' @param reactions list of state transitions. Each reaction is itself a list
#' that requires a `transition` vector and a `rate` vector
#' @param sim_scenario list that defines the scenario including `params`,
#' `initial_states` and `sim_args` that must at least include simulation time
#' `T`
#' @name stochastic_model
#' @examples
#' reactions <- list(
#'   infection = list(
#'     transition = c("I" = +1),
#'     rate = function(x,p,t){p$beta})
#' )
#' example_scenario <- list(
#'   params = list(beta = 0.1),
#'   initial_states = list(I = 0),
#'   sim_args = list(T = 10)
#' )
#' stochastic_model(reactions,example_scenario)
#'
#' @export
stochastic_model <- function(reactions,sim_scenario){
  validate_sim_scenario(sim_scenario)
  validate_reactions(reactions)
  new_stochastic_model(reactions,sim_scenario)
}

#' Print method for stochastic_model objects
#'
#' This function provides a formatted and colorized summary of a stochastic model,
#' including initial values, parameters, and a list of reactions with their transitions
#' and rate function arguments.
#'
#' @param x An object of class \code{stochastic_model}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the \code{x} object.
#' @export
#' @method print stochastic_model
print.stochastic_model <- function(x, ...) {
  print_list <- function(x){
    for(i in seq_along(x)){
      name <- names(x)[i]
      cat(" ", crayon::bold(name),":", signif(x[[name]], digits = 3))
      if(i %% 3 == 0){
        cat("\n")
      }
    }
  }

  cat(crayon::bold$blue("Stochastic Model"), "\n")
  cat(crayon::blue(strrep("=", 18)), "\n")

  # Initial values
  cat("\n", crayon::green$bold("Initial Values:"), "\n", sep = "")
  print_list(x$initial_states)

  # Parameters
  cat("\n", crayon::green$bold("Parameters:"), "\n", sep = "")
  print_list(x$params)

  # sim arguments
  cat("\n", crayon::green$bold("Simulation arguments:"), "\n", sep = "")
  print_list(x$sim_args)

  # Reactions
  cat("\n", crayon::green$bold("Reactions:"), "\n", sep = "")
  transitions <- x$transitions
  for (name in names(transitions)) {
    transition <- transitions[[name]]
    trans_str <- paste(
      paste0(names(transition),
             ifelse(transition > 0, " +", " "),
             transition),
      collapse = ", "
    )
    cat(crayon::yellow$bold(sprintf(" - %s:\n", name)))
    cat(crayon::silver(sprintf("   Transition: [%s]\n", trans_str)))
  }

  invisible(x)
}



#' run stochastic model using adaptive tau set algorithm
#' @param x An object of class \code{stochastic_model}.
#' @param tf maximum time
#' @export
#' @rdname stochastic_model
run_sim <- function(x,tf=NULL){
  init_values <- unlist(x$initial_states)
  if(is.null(tf)){
    tf <- x$sim_args$T
  }
  adaptivetau::ssa.adaptivetau(init_values, x$transitions, x$rates,
                               x$params, tf=tf)
}

#' Interpolate State Trajectories to Whole Days
#'
#' Converts a time-series data frame or matrix with continuous time points
#' and one or more state columns into a tibble with daily interpolated values.
#'
#' @param x A matrix or data frame with a numeric column named "time"
#'   and one or more state columns with arbitrary names. output of `run`
#'
#' @return A tibble with one row per whole day and interpolated values for each state.
#'
#' @examples
#' mat <- cbind(
#'   time = c(0.1, 0.5, 1.2, 2.4, 3.8),
#'   S = c(1000, 900, 800, 600, 500),
#'   I = c(1, 10, 50, 100, 90)
#' )
#' interpolate_run_by_day(mat)
#'
#' @export
interpolate_run_by_day <- function(x) {

  x <- as.data.frame(x)
  stopifnot("time" %in% names(x))

  time_range <- floor(min(x$time)) : ceiling(max(x$time))
  states <- setdiff(names(x), "time")

  interpolated <- purrr::map(states, function(state_name) {
    out <- data.frame(stats::approx(x$time, x[[state_name]],
                             xout = time_range, method = "linear")$y)
    names(out) <- state_name
    out
  }) |>
    purrr::list_cbind()

  dplyr::tibble(time = time_range) |>
    dplyr::bind_cols(interpolated)
}



