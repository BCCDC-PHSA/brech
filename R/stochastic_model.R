###################################
##                               ##
## Parameter and state creation  ##
##                               ##
###################################

#' load model parameters
#' @param file file path to yaml
#' @export
load_model_params <- function(file=NULL){
  if(is.null(file)){
    file <- system.file("extdata", "params.yaml", package = "brech")
  }
  conditions <- yaml::read_yaml(file)
  params <- conditions$params
  # preserves vector if using age groups
  N <- Reduce('+',conditions$initial_states)

  params$theta <- 1 / params$latent_period
  params$gamma <- 1 / params$infectious_period
  params$case_rate <- 1 / params$case_delay

  # if beta isn't defined then base on R0 and population size
  if(!("beta" %in% names(params))){
    params$beta <- params$r_0 * params$gamma / N
  }
  params$N <- N
  conditions$params <- params
  # init states
  # set initial cases to 0 by convention
  # multiply by another state to retain correct shape
  first_state <- names(conditions$initial_states)[1]
  conditions$initial_states$C <- 0*conditions$initial_states[[first_state]]
  conditions$initial_states$D <- 0*conditions$initial_states[[first_state]]
  return(conditions)
}

#' change vaccination coverage in initial conditions
#' @param sim_scenario output of `load_model_params`
#' @param vaccine_rate numeric proportion vaccinated. Can be a numeric
#'   vector if using age groups. If named then vaccinates according to each
#' age group with format "S_\{age_group\}"
#' @export
vaccinate_initial_conditions <- function(sim_scenario,vaccine_rate){
  age_groups <- names(vaccine_rate)
  initial_states <- sim_scenario$initial_states

  # if not named
  if(is.null(age_groups)){
    S <- (initial_states$S + initial_states$V)
    initial_states$V <- round(S*vaccine_rate)
    initial_states$S <- round(S*(1 - vaccine_rate))
    # add back in discrepancy due to rounding
    initial_states$S <- initial_states$S  + (S - initial_states$V - initial_states$S)

  }else{
    for(age_group in age_groups){
      age_vaccine_rate <- vaccine_rate[[age_group]]
      S <- initial_states[[paste0("S_",age_group)]]
      V <- initial_states[[paste0("V_",age_group)]]
      S <- N <- S + V
      # vaccinate
      V <- round(S*age_vaccine_rate)
      S <- N - V
      # update initial states
      initial_states[[paste0("S_",age_group)]] <- S
      initial_states[[paste0("V_",age_group)]] <- V
    }
  }
  sim_scenario$initial_states <- initial_states
  return(sim_scenario)
}

#' Convert a vector of each state into a labeled state using
#' `age_groups` labels
#' @param sim_scenario output of `load_model_params`
#' @param age_groups vector of age group labels
#' @examples
#' sim_scenario <- list(
#'   "initial_states" = list(
#'     "S" = c(1,2,3),
#'     "I" = c(0,0,1)
#'   )
#' )
#' age_groups <- c("0-5","5-19","19+")
#' create_age_initial_conditions(sim_scenario,age_groups)
#' @export
create_age_initial_conditions <- function(sim_scenario,age_groups){
  initial_states <- sim_scenario$initial_states
  updated_initial_states <- list()
  for(state in names(initial_states)){
    for(i in seq_along(age_groups)){
      age_group <- age_groups[i]
      updated_initial_states[[paste0(state,"_",age_group)]] <- initial_states[[state]][i]
    }
  }
  sim_scenario$initial_states <- updated_initial_states
  return(sim_scenario)
}

#' get reaction list from number of age groups m and
#' contact_matrix
#' @param contact_matrix matrix m x m
#' @export
build_age_reactions <- function(contact_matrix){
  m <- ncol(contact_matrix)
  age_groups <- colnames(contact_matrix)
  if(m != nrow(contact_matrix)){
    stop("contact_matrix needs the same number of rows and columns.")
  }

  reactions <- list()
  for (i in seq_along(age_groups)){
        age_group <- age_groups[i]
        ######################################
        #  detectable infection transitions  #
        ######################################
        transition <- function() {
          trans <- c()
          trans[paste0("S_",age_group)] <- -1
          trans[paste0("E_",age_group)] <- +1
          trans[paste0("D_",age_group)] <- +1
          trans[paste0("INC_",age_group)] <- +1
          trans
        }
        detectable_infection_rate <-
          function(i,age_group){
            force(i)
            force(age_group)
            function(x, p, t) {
              I <- x[paste0("I_", age_groups)]
              S <- x[paste0("S_", age_group)] # just need the ith S component
              lambda <- p$case_ascertainment * p$beta * S * sum(contact_matrix[,i] * (I / p$N))
              lambda  # return i-th force of infection term
            }
          }
        reactions[[paste0("detectable_infection_",age_group)]] <- list(
          "transition"=transition(),
          "rate" = detectable_infection_rate(i,age_group)
        )
        ########################################
        #  undetectable infection transitions  #
        ########################################
        transition <- function() {
          trans <- c()
          trans[paste0("S_",age_group)] <- -1
          trans[paste0("E_",age_group)] <- + 1
          trans[paste0("INC_",age_group)] <- +1
          trans
        }
        undetectable_infection_rate <-
          function(i,age_group){
            force(i)
            force(age_group)
            function(x, p, t) {
              I <- x[paste0("I_", age_groups)]
              S <- x[paste0("S_", age_group)] # just need the ith S component
              lambda <- (1 - p$case_ascertainment) * p$beta * S * sum(contact_matrix[,i] * (I / p$N))
              lambda  # return i-th force of infection term
            }
          }
        reactions[[paste0("undetectable_infection_",age_group)]] <- list(
          "transition"=transition(),
          "rate" = undetectable_infection_rate(i,age_group)
        )
        ############################
        #  incubation transitions  #
        ############################
        transition <- function() {
          trans <- c()
          trans[paste0("E_",age_group)] <- -1
          trans[paste0("I_",age_group)] <- + 1
          trans
        }
        incubation_rate <-
        function(age_group){
          force(age_group)
          function(x, p, t) {
            E <- x[paste0("E_", age_group)] # ith component
            p$theta * E
          }
        }
        reactions[[paste0("incubation_",age_group)]] <- list(
          "transition"=transition(),
          "rate" = incubation_rate(age_group)
        )
        ##############################
        #    recovery transitions    #
        ##############################
        transition <- function() {
          trans <- c()
          trans[paste0("I_",age_group)] <- -1
          trans[paste0("R_",age_group)] <- + 1
          trans
        }
        recovery_rate <-
        function(age_group){
          force(age_group)
          function(x, p, t) {
            I <- x[paste0("I_", age_group)] # ith component
            p$gamma * I
          }
        }
        reactions[[paste0("recovery_",age_group)]] <- list(
          "transition"=transition(),
          "rate" = recovery_rate(age_group)
        )
        ##############################
        #        case detection      #
        ##############################
        transition <- function() {
          trans <- c()
          trans[paste0("C_",age_group)] <- + 1
          trans[paste0("D_",age_group)] <- - 1
          trans
        }
        case_detection_rate <-
        function(age_group){
          force(age_group)
          function(x, p, t) {
            D <- x[paste0("D_", age_group)] # ith component
            p$case_rate * D
          }
        }
        reactions[[paste0("case_detection_",age_group)]] <- list(
          "transition"=transition(),
          "rate" = case_detection_rate(age_group)
        )
        ##############################
        #              PEP           #
        ##############################
        transition <- function() {
          trans <- c()
          trans[paste0("E_",age_group)] <- - 1
          trans[paste0("R_",age_group)] <- + 1
          trans
        }
        pep_administered_rate <-
          function(age_group){
            force(age_group)
            function(x, p, t) {
              E <- x[paste0("E_", age_group)] # ith component
              p$pep_rate * E
            }
          }
        reactions[[paste0("pep_administered_",age_group)]] <- list(
          "transition"=transition(),
          "rate" = pep_administered_rate(age_group)
        )
        ###############################
        #        hospitalization      #
        ###############################
        transition <- function() {
          trans <- c()
          trans[paste0("H_",age_group)] <- + 1
          trans[paste0("I_",age_group)] <- - 1
          trans
        }
        hospitalization_rate <-
          function(i, age_group){
            force(i)
            force(age_group)
            function(x, p, t) {
              I <- x[paste0("I_", age_group)] # ith component
              prop_hosp <- p$hospitalization_proportion[i]
                p$gamma * (prop_hosp/(1-prop_hosp)) * I
            }
          }
        reactions[[paste0("hospitalization_",age_group)]] <- list(
          "transition"=transition(),
          "rate" = hospitalization_rate(i, age_group)
        )

    }

  return(reactions)

}

#' Using the next generation matrix approach, calculates beta based on a fixed
#' R0
#' @param sm Object of class \code{stochastic_model}
#' @param contact_matrix An m x m contact matrix
#' @param population_sizes An m vector of population size
#' @param R0 numeric reproduction number
#' @return numeric
#' @export
calculate_age_structured_beta <- function(sm, contact_matrix, population_sizes, R0){
  # beta calculation based on R0
  recovery_rate <- sm$params$gamma
  adjusted_matrix <- contact_matrix * outer(population_sizes, 1/population_sizes)
  spectral_radius <- max(Mod(eigen(adjusted_matrix)$values))
  beta <- recovery_rate*R0/spectral_radius

  return(beta)
}

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

#####
##### Model fitting section
#####

#' update parameters for a stochastic_model
#' @param x An object of class \code{stochastic_model}.
#' @param ... parameters to update
#' @param params parameters to update passed in as a list
#' @return An object of class \code{stochastic_model}.
#' @export
#' @rdname stochastic_model
#' @examples
#' \dontrun{
#' sm <- update_parameters(sm,beta=0.2)
#' }
update_parameters <- function(x,...,params=list()){
  args <- utils::modifyList(params, list(...))
  for(param in names(args)){
    x$params[[param]] <- args[[param]]
  }
  return(x)
}

#' update initial state for a stochastic_model
#' @param x An object of class \code{stochastic_model}.
#' @param state data.frame
#' @rdname stochastic_model
#' @export
update_state <- function(x,state){
  x$initial_states <- unlist(state)
  return(x)
}


#' Compute Daily Incident Cases from Cumulative Simulation Output
#'
#' This function takes a data frame or matrix containing irregularly spaced simulation
#' time points and cumulative case counts (`C`), and returns the number of new cases
#' per day using linear interpolation.
#'
#' @param r A data frame or matrix with columns `"time"` and `"C"`, representing
#'   time (in continuous units) and cumulative number of cases, respectively.
#'
#' @return A tibble with two columns:
#'   \describe{
#'     \item{day}{Integer day values from the floor of the minimum time to the ceiling of the maximum time}
#'     \item{daily_incidence}{Number of new cases per day, computed as the difference in interpolated cumulative cases}
#'   }
#'
#' @details
#' Since simulation outputs may record cumulative cases at non-integer time points,
#' this function performs linear interpolation to estimate cumulative cases at each
#' whole-number day. It then computes the number of incident cases as the first difference
#' of the interpolated values. The first day is assigned 0 new cases by default.
#'
#' @examples
#' sim_data <- as.matrix(dplyr::tibble(
#'   time = c(0, 0.1, 0.9, 1.4, 2.3, 3.8),
#'   C = c(0, 1, 5, 7, 10, 15)
#' ))
#' get_daily_cases(sim_data)
#' @export
get_daily_cases <- function(r){
  C_interp <- NULL
  day <- daily_incidence <- NULL
  times <- r[,"time"]
  cases <- rowSums(r[,startsWith(colnames(r),"C"),drop=FALSE])
  dplyr::tibble(day = floor(min(times)) : ceiling(max(times))) |>
    dplyr::mutate(
      C_interp = stats::approx(x = times, y = cases, xout = day, method = "linear")$y
    ) |>
    dplyr::mutate(
      daily_incidence = c(0, diff(C_interp))  # First day has no prior to diff from
    ) |>
    dplyr::select(day,daily_incidence)
}


#' Plot Prior and Posterior Distributions with Optional True Value
#'
#' This function generates a density plot comparing the prior and posterior distributions
#' of a parameter. An optional vertical dashed line can be added to indicate the true parameter value.
#'
#' @param prior_values A numeric vector of samples from the prior distribution.
#' @param posterior_values A numeric vector of samples from the posterior distribution.
#' @param true_value Optional numeric value representing the true parameter value to be shown as a vertical dashed line.
#'
#' @return A \code{ggplot2} object representing the density plot comparing the prior and posterior distributions.
#'
#' @details The function uses `geom_density()` to visualize both prior and posterior samples, and `geom_vline()` to optionally display the true parameter value. It uses `theme_minimal()` and places the legend at the bottom for clarity.
#'
#' @examples
#' prior <- rnorm(1000, mean = 0, sd = 2)
#' posterior <- rnorm(1000, mean = 1, sd = 0.5)
#' true_value <- 1.2
#' plot_prior_posterior(prior, posterior, true_value)
#'
#' @export
plot_prior_posterior <- function(prior_values,posterior_values,true_value = NULL){
  value <- type <- NULL
  g <- dplyr::bind_rows(
    dplyr::tibble(value = prior_values, type = "Prior"),
    dplyr::tibble(value = posterior_values, type = "Posterior")
  ) |>
    ggplot2::ggplot(ggplot2::aes(x = value, fill = type, color = type)) +
    ggplot2::geom_density(alpha = 0.4, linewidth = 1)

  if(!is.null(true_value)){
    g <- g +
      ggplot2::geom_vline(xintercept = true_value,
                          linetype = "dashed", color = "black", linewidth = 1)
  }
  g <- g +
    ggplot2::labs(
      x = "Parameter Value",
      y = "Density",
      fill = "Distribution",
      color = "Distribution"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  return(g)
}

#' Approximate Bayesian Computation for Stochastic Model
#'
#' Runs an Approximate Bayesian Computation (ABC) procedure using a stochastic model and a set of prior parameter draws.
#'
#' @param sm An object of class \code{stochastic_model}.
#' @param priors A data frame or tibble of prior parameter values. Each row should correspond to one simulation, and columns should match the parameter names expected by `update_parameters()`.
#' @param target A named vector or data frame of observed summary statistics to match against simulated summary statistics.
#' @param stat_func A function that computes summary statistics from a simulation result. Should accept a model object and return a named numeric vector or tibble row.
#' @param states Optional. A data frame or tibble of initial states. Each row corresponds to a complete model state e.g. for an SIR model each row has three columns
#' @param tol A numeric tolerance (between 0 and 1) for the ABC rejection algorithm. Defaults to 0.1.
#' @param method A character string indicating the ABC method to use. Defaults to `"rejection"`. Passed directly to [abc::abc()].
#' @param seed An integer specifying the random seed for reproducibility. Defaults to 123.
#'
#' @return An object of class `"abc_stochastic_model"` containing:
#' \describe{
#'   \item{model}{The original stochastic model.}
#'   \item{fit}{An object returned by [abc::abc()] representing the ABC posterior.}
#'   \item{priors}{The input prior values used for simulations.}
#'   \item{posteriors}{The output posterior values used for simulations.}
#'   \item{prior_predictive}{The summary statistics defined in `stat_func`
#'   drawn from the prior distribution.}
#'   \item{state}{Final state of model.}
#' }
#'
#' @details
#' This function performs ABC by:
#' \enumerate{
#'   \item Sampling parameter combinations from `priors`
#'   \item Updating the `stochastic_model` with each sampled parameter set
#'   \item Running the model and computing summary statistics via `stat_func`
#'   \item Comparing simulated statistics to `target` using the `abc::abc()` function
#' }
#'
#' Parallelization is handled using `furrr::future_map()`, allowing reproducible simulation across parameter sets using a fixed seed.
#'
#' @examples
#' \dontrun{
#' priors <- expand.grid(beta = seq(0.1, 0.5, length.out = 10))
#' target <- c(mean_infected = 50)
#' abc_result <- abc_stochastic_model(
#'   stochastic_model = my_model,
#'   priors = priors,
#'   target = target,
#'   stat_func = function(sim) data.frame(mean_infected = mean(sim$I)),
#'   tol = 0.2
#' )
#' }
#'
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr pmap map reduce
#' @importFrom dplyr bind_rows
#' @importFrom abc abc
#' @export
abc_stochastic_model <- function(sm,priors,target,stat_func,
                                 states=NULL,
                                 tol=0.1,method="rejection",seed=123){

  prior_list <- priors |> purrr::pmap(~ list(...))

  # Create a list of state updates (or NULLs if states is NULL)
  if (!is.null(states)) {
    state_list <- states |> purrr::pmap(~ list(...))
  } else {
    state_list <- rep(list(NULL), length(prior_list))
  }


  progress_update <- progressr::progressor(steps=length(prior_list))

  lists_of_results <-
    furrr::future_map2(
      prior_list,
      state_list,
      function(prior,state){
        m <- update_parameters(sm, params = prior)
        if(!is.null(state)){
          m <- update_state(m, state)
        }
        m <- run_sim(m)
        final_state <- m[nrow(m), -1]
        stat <- stat_func(m)
        progress_update()

        list("stat" = stat,
             "final_state" = final_state)
    },
    .options = furrr::furrr_options(seed = seed))

    stat <- purrr::map(lists_of_results, "stat") |>
      purrr::reduce(dplyr::bind_rows)
    final_state <- purrr::map(lists_of_results, "final_state") |>
      purrr::reduce(dplyr::bind_rows)
    # simulation <- purrr::map(lists_of_results, "simulation") |>
    #   purrr::map(as.data.frame) |>
    #   purrr::imap(~dplyr::mutate(.x, index = .y)) |>
    #   dplyr::bind_rows()



  fit <- abc::abc(target=target,param=priors,sumstat=stat,
                  tol=tol,method=method)

  # get posterior values
  final_state <- final_state[fit$region,]
  posteriors <- priors[fit$region,]
  indexes <- seq_len(length(fit$region))[fit$region]
  # simulation <- simulation |> dplyr::filter(index %in% indexes)

  structure(list("model" = sm, "fit" = fit,
                 "priors" = priors, "posteriors" = posteriors,
                 "prior_predictive" = stat,
                 "state" = final_state),
                 #"projection" = simulation),
            class="abc_stochastic_model")
}

#' type of plot either "prior_posterior" or "simulations"
#' @param x An object of class \code{abc_stochastic_model}.
#' @param ... list of additional plot arguments.
#' `param` - if `type` is "prior_posterior" then parameter to plot
#' `state` - if `type` is "simulations" then state to plot
#' @export
plot.abc_stochastic_model <- function(x, ...){

  asm <- x
  args <- list(...)
  param <- args$param
  state <- args$state

  if(!is.null(param)){
    plot_prior_posterior(dplyr::pull(asm$priors,param),
                         dplyr::pull(asm$posteriors,param),
                         true_value = asm$model$params[[param]])
  }else if(!is.null(state)){
    plot_projections(asm$simulation,state)
  }else{
    stop(paste0("plot type not recognized"))
  }
}


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


###################################
##                               ##
##   Model projections section   ##
##                               ##
###################################

#' Project a Stochastic Model Forward in Time
#'
#' Runs forward projections from the posterior draws of a fitted
#' \code{abc_stochastic_model} for a fixed number of time steps (default 30),
#' using region-specific states and parameter sets. Each simulation is seeded
#' for reproducibility and run in parallel.
#'
#' @param asm An object of class \code{abc_stochastic_model} or \code{scenario_stochastic_model}.
#' @param project_time Number of time steps to project forward. Defaults to
#'  \code{30}.
#' @param seed An integer used to seed the parallel simulations for
#'   reproducibility. Defaults to \code{123}.
#'
#' @return An object of class \code{projection_stochastic_model}, which is a list
#'   with the following elements:
#'   \item{model}{The stochastic model used for simulation (of class \code{stochastic_model}).}
#'   \item{projection}{A data frame with the projected simulations concatenated across
#'   posterior draws and days, including any parameter updates appended as columns.}
#'
#' @details
#' This function takes a fitted \code{abc_stochastic_model}, extracts the
#' posterior state and parameter sets, then runs simulations in parallel using
#' \pkg{furrr}. The model is updated with each posterior draw, run forward for a
#' fixed time frame, interpolated to daily resolution, and combined
#' into a single projection data frame.
#'
#' @importFrom furrr future_map2 furrr_options
#' @importFrom purrr pmap reduce
#' @importFrom dplyr bind_cols bind_rows
#' @name projection_stochastic_model
#' @export
projection_stochastic_model <- function(asm, project_time = 30,
                                        seed = 123){
  plist <- function(df){
    df |> purrr::pmap(~ list(...))
  }
  if("abc_stochastic_model" %in% class(asm)){
    state_list <- plist(asm$state)
    project_params <- plist(asm$posteriors)
  }else if("scenario_stochastic_model" %in% class(asm)){
    state_list <- plist(asm$states)
    project_params <- plist(asm$parameters)
  }else{
    stop(paste0("Object of class ",class(asm)," not known."))
  }

  progress_update <- progressr::progressor(steps=length(project_params))

  projections <- furrr::future_map2(
    .x = project_params,
    .y = state_list,
    .f = function(params,state){
      progress_update()

      asm$model |>
        update_state(state) |>
        update_parameters(params=params) |>
        run_sim(tf=project_time) |>
        interpolate_run_by_day() |>
        dplyr::bind_cols(as.data.frame(params))
    },
    .options = furrr::furrr_options(seed = seed)) |>
    purrr::reduce(dplyr::bind_rows)

  structure(list("model"=asm$model,"projection"=projections),
            class="projection_stochastic_model")
}

#' Add date to projection using a `start_date`
#' @rdname projection_stochastic_model
#' @param psm An object of class \code{projection_stochastic_model}
#' @param start_date a date in yyyy-mm-dd format
#' @return An object of class \code{projection_stochastic_model}
#' @export
add_projection_date <- function(psm,start_date){
  time <- NULL
  start_date <- lubridate::ymd(start_date)
  psm$projection <- psm$projection |>
    dplyr::mutate(date = start_date + lubridate::days(time))
  return(psm)
}

#' plot projections
#' @description
#' If state doesn't exist e.g. "I" but has states "I_\{age\}" then
#' will sum over all states first before plotting
#'
#' @rdname projection_stochastic_model
#' @param x An object of class \code{projection_stochastic_model}
#' @param ... other arguments including `state` state to plot, `"I"` by default
#' and `type` either show individual trajectories "samples" or summarize "quantiles",
#' `"quantiles"` by default
#' @export
plot.projection_stochastic_model <- function(
    x,
    ...
  ){

  psm <- x
  args <- list(...)
  state <- args$state
  type <- args$type

  # default args
  if(is.null(state)){
    state <- "I"
  }
  if(is.null(type)){
    type <- "quantiles"
  }

  if(!(state %in% names(psm$projection))){
    projection <- collapse_states(psm)
  }else{
    projection <- psm$projection
  }

  if(type == "quantiles"){
    plot_projections(projection,state)
  }else if(type == "samples"){
    plot_projection_samples(projection,state)
  }else{
    stop("`type` should be 'quantiles' or 'samples'")
  }

}

#' @rdname projection_stochastic_model
#'
plot_projections <- function(projection,state){
  time <- date <- q0.025 <- q0.25 <- q0.5 <- q0.75 <- q0.975 <- NULL
  .data <- NULL
  probs <- c(0.025,0.25, 0.5, 0.75,0.975)
  time_col <- if("date" %in% colnames(projection)) "date" else "time"
  plot_data <- create_projection_quantiles(projection, probs = probs)
  plot_data |>
    tidyr::pivot_wider(id_cols=tidyselect::all_of(time_col),
                       names_from="quantile",
                       values_from=tidyselect::all_of(state),
                       names_prefix = "q") |>
    ggplot2::ggplot(ggplot2::aes(x=.data[[time_col]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=q0.025,ymax=q0.975),
                         fill="deepskyblue", alpha = 0.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=q0.25,ymax=q0.75),
                         fill="deepskyblue", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y=q0.5), color="black") +
    ggplot2::theme_minimal() +
    ggplot2::ylab(state)
}

#' Plot projection samples as opposed to quantiles.
#' Color by max outcome per sample.
#' @rdname projection_stochastic_model
#'
plot_projection_samples <- function(projection,state){
  time <- date <- q0.025 <- q0.5 <- q0.75 <- q0.975 <- NULL
  group <- .data <- max_state <- NULL
  probs <- c(0.5)
  time_col <- if("date" %in% colnames(projection)) "date" else "time"
  median_data <- create_projection_quantiles(projection, probs = probs) |>
    tidyr::pivot_wider(id_cols=tidyselect::all_of(time_col),
                       names_from="quantile",
                       values_from=tidyselect::all_of(state),
                       names_prefix = "q")

  get_projection_dataframe(projection) |>
  dplyr::mutate(group = cumsum(time == 0)) |>
  dplyr::group_by(group) |>
  dplyr::mutate(max_state = max(.data[[state]])) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    max_state = dplyr::ntile(max_state, 4),
    max_state = paste0("q", max_state)
    ) |>
  ggplot2::ggplot(ggplot2::aes(x=.data[[time_col]])) +
    ggplot2::geom_line(ggplot2::aes(y=.data[[state]], color=max_state,
                                    group = group),
                       alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(x=.data[[time_col]], y=q0.5),
                       color = "black", size=1.5, data = median_data) +
    ggplot2::scale_color_viridis_d(option = "viridis") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ylab(state)
}

#' @rdname projection_stochastic_model
#' @param projection projection grouped using [projection_quantiles_by_age_group]
#' @param state state to plot
#' @examples
#' \dontrun{
#'   projection <- psm |> projection_quantiles_by_age_group("C")
#'   plot_projections_by_age_group(projection,"C")
#' }
#' @export
plot_projections_by_age_group <- function(projection,state){
  age_group <- q0.025 <- q0.25 <- q0.5 <- q0.75 <- q0.975 <- NULL
  .data <- NULL
  time_col <- if("date" %in% colnames(projection)) "date" else "time"
  projection |>
    tidyr::pivot_wider(id_cols=c(tidyselect::all_of(time_col),"age_group"),
                       names_from="quantile",
                       values_from="values",
                       names_prefix = "q") |>
    ggplot2::ggplot(ggplot2::aes(x=.data[[time_col]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=q0.025,ymax=q0.975),
                         alpha = 0.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=q0.25,ymax=q0.75),
                         alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y=q0.5)) +
    ggplot2::facet_wrap(ggplot2::vars(age_group)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(y=state)
}

#' @rdname projection_stochastic_model
#' @param psm object of type \code{projection_stochastic_model} or tibble
#' @param probs vector of quantiles to plot
#' @param by_cols vector of columns to group by
#' @export
create_projection_quantiles <- function(
    psm,
    probs=c(0.025,0.25, 0.5, 0.75,0.975),
    by_cols = c("time","date")){
  projection <- get_projection_dataframe(psm)

  project_time <- length(unique(projection$time))

  projection |>
    dplyr::reframe(
      dplyr::across(
        .cols = tidyselect::everything(),
        .fns = ~ quantile(.x, probs = probs, na.rm = TRUE)
      ),
      .by = tidyselect::any_of(by_cols)
    ) |>
    dplyr::mutate(quantile=rep(probs,project_time))
}

#' @rdname projection_stochastic_model
#' @param psm object of type \code{projection_stochastic_model} or tibble
#' @param state model state e.g. "C"
#' @export
projection_quantiles_by_age_group <- function(psm,state,
  probs=c(0.025,0.25, 0.5, 0.75,0.975)
  ){
  by_cols <- c("date", "time", "age_group")
  projection <- get_projection_dataframe(psm)
  projection <- projection |>
    create_age_group_column(state) |>
    dplyr::reframe(
      quantile = probs,
      values = purrr::map_dbl(probs, ~ quantile(.data[[state]], probs = .x, na.rm = TRUE)),
      .by = tidyselect::any_of(by_cols)
    )


  return(projection)

}

#' projection tibble by age group
#' @rdname projection_stochastic_model
#' @param psm  \code{projection_stochastic_model}
#' @param state a model state e.g. "C" for cases
#' @return tibble
#' @export
projection_by_age_group <- function(psm,state){
  projection <- get_projection_dataframe(psm)
  projection |>
    create_age_group_column(state)
}

#' create age group column for \code{projection_stochastic_model$projection}.
#' @rdname projection_stochastic_model
#' @param psm a tibble from \code{projection_stochastic_model} or \code{projection_stochastic_model}
#' @param state a model state e.g. "C" for cases
#' @param index_cols columns to preserve defaults to \code{c("time","date")}
#' @export
create_age_group_column <- function(psm, state,
                                    index_cols = c("date","time")){
  projection <- get_projection_dataframe(psm)
  projection |>
    dplyr::select(c(
      tidyselect::starts_with(paste0(state,"_")),
      tidyselect::any_of(index_cols)
    )) |>
    tidyr::pivot_longer(
      cols = tidyselect::starts_with(paste0(state,"_")),
      names_to = "age_group",
      names_prefix = paste0(state,"_"),
      values_to = state
    )
}

#' collapse states of a \code{projection_stochastic_model} for plotting
#' @rdname projection_stochastic_model
#' @param psm Object of class \code{projection_stochastic_model}
#' @return tibble
#' @export
collapse_states <- function(psm){
  projection <- get_projection_dataframe(psm)

  # Exclude time column
  other_cols <- setdiff(names(projection), c("time", "date"))

  # Extract MATCH part (prefix before the underscore)
  match_parts <- stringr::str_extract(other_cols, "^[^_]+")

  # Group columns by prefix
  grouped <- split(other_cols, match_parts)

  # Build the output
  for (prefix in names(grouped)) {
    projection[[prefix]] <- projection |>
      dplyr::select(tidyselect::all_of(grouped[[prefix]])) |>
      rowSums(na.rm=TRUE)
  }

  return(projection)
}

#' Take difference of columns
#' @rdname projection_stochastic_model
#' @description
#' take the difference of column(s) that match on state `state`. This is
#' useful for any cumulative columns e.g. `INC` and `C` by converted them into
#' the incidence and case incidence respectively
#' @param psm Object of class \code{projection_stochastic_model}
#' @param state character string e.g. "C"
#' @return tibble
#' @export
difference_of_states <- function(psm,state){
  time <- sim_id <- NULL
  projection <- get_projection_dataframe(psm)
  projection |>
    dplyr::mutate(sim_id = cumsum(time == 0)) |>
    dplyr::group_by(sim_id) |>
    dplyr::mutate(
      dplyr::across(tidyselect::starts_with(paste0(state, "_")),
      ~ c(0, diff(.)))
    ) |>
    dplyr::ungroup()
}

#' Convert cumulative state variables back to zero at a reference time
#'
#' @description
#' For cumulative state variables such as cases (e.g., `C_[0,1)`), this function
#' resets their values to 0 at a specified reference time (default is 0). This is
#' useful when wanting to compute cumulative changes from a particular time point
#' in a simulation or projection.
#'
#' @param psm An object of class `projection_stochastic_model`, or a compatible
#' data structure as accepted by `get_projection_dataframe()`.
#' @param state A character string specifying the prefix of the cumulative state variable
#' to reset (e.g., `"C"` to match columns like `"C_[0,1)"`, `"C_[1,5)"`, etc.).
#' @param reset_time A numeric value indicating the time at which to reset the cumulative
#' values back to zero. Defaults to 0.
#'
#' @return A `tibble` with the same structure as the original projection data frame
#' @examples
#' \dontrun{
#' # Reset cumulative cases to 0 at time 10
#' reset_state(psm = fitted_model, state = "C", reset_time = 10)
#'
#' # Reset cumulative hospitalizations
#' reset_state(psm = fitted_model, state = "H", reset_time = 5)
#' }
#' @export
reset_state <- function(psm, state, reset_time = 0){
  time <- sim_id <- NULL
  projection <- get_projection_dataframe(psm)
  projection |>
    dplyr::mutate(sim_id = cumsum(time == 0)) |>
    dplyr::group_by(sim_id) |>
    dplyr::mutate(
      dplyr::across(tidyselect::starts_with(paste0(state, "_")),
                    ~ (. - .[time == reset_time]))
    ) |>
    dplyr::ungroup()
}

#' get projection tibble
#' Internal function used to either allow a projection_stochastic_model object
#' or the underlying projection tibble to be included in functions like
#' `collapse_states`. Useful for piping multiple functions that modify the
#' projections without updating the underlying object
#' @noRd
get_projection_dataframe <- function(psm) {
  cls <- class(psm)

  if ("projection_stochastic_model" %in% cls) {
    return(psm$projection)
  } else if ("tbl_df" %in% cls || "data.frame" %in% cls) {
    return(psm)
  } else {
    stop("Unsupported object class: ", paste(cls, collapse = ", "))
  }
}
