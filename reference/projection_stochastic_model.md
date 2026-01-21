# Project a Stochastic Model Forward in Time

Runs forward projections from the posterior draws of a fitted
`abc_stochastic_model` for a fixed number of time steps (default 30),
using region-specific states and parameter sets. Each simulation is
seeded for reproducibility and run in parallel.

If state doesn't exist e.g. "I" but has states "I\_{age}" then will sum
over all states first before plotting

take the difference of column(s) that match on state `state`. This is
useful for any cumulative columns e.g. `INC` and `C` by converted them
into the incidence and case incidence respectively

## Usage

``` r
projection_stochastic_model(asm, project_time = 30, seed = 123)

add_projection_date(psm, start_date)

# S3 method for class 'projection_stochastic_model'
plot(x, ...)

plot_projections(projection, state)

plot_projection_samples(projection, state)

plot_projections_by_age_group(projection, state)

create_projection_quantiles(
  psm,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  by_cols = c("time", "date")
)

projection_quantiles_by_age_group(
  psm,
  state,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)

projection_by_age_group(psm, state)

create_age_group_column(psm, state, index_cols = c("date", "time"))

collapse_states(psm)

difference_of_states(psm, state)
```

## Arguments

- asm:

  An object of class `abc_stochastic_model` or
  `scenario_stochastic_model`.

- project_time:

  Number of time steps to project forward. Defaults to `30`.

- seed:

  An integer used to seed the parallel simulations for reproducibility.
  Defaults to `123`.

- psm:

  Object of class `projection_stochastic_model`

- start_date:

  a date in yyyy-mm-dd format

- x:

  An object of class `projection_stochastic_model`

- ...:

  other arguments including `state` state to plot, `"I"` by default and
  `type` either show individual trajectories "samples" or summarize
  "quantiles", `"quantiles"` by default

- projection:

  projection grouped using projection_quantiles_by_age_group

- state:

  character string e.g. "C"

- probs:

  vector of quantiles to plot

- by_cols:

  vector of columns to group by

- index_cols:

  columns to preserve defaults to `c("time","date")`

## Value

An object of class `projection_stochastic_model`, which is a list with
the following elements:

- model:

  The stochastic model used for simulation (of class
  `stochastic_model`).

- projection:

  A data frame with the projected simulations concatenated across
  posterior draws and days, including any parameter updates appended as
  columns.

An object of class `projection_stochastic_model`

tibble

tibble

tibble

## Details

This function takes a fitted `abc_stochastic_model`, extracts the
posterior state and parameter sets, then runs simulations in parallel
using furrr. The model is updated with each posterior draw, run forward
for a fixed time frame, interpolated to daily resolution, and combined
into a single projection data frame.

## Examples

``` r
if (FALSE) { # \dontrun{
  projection <- psm |> projection_quantiles_by_age_group("C")
  plot_projections_by_age_group(projection,"C")
} # }
```
