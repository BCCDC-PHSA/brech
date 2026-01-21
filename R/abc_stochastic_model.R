
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
