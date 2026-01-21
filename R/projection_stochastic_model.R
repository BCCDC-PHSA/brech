
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
