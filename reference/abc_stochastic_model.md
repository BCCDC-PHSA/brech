# Approximate Bayesian Computation for Stochastic Model

Runs an Approximate Bayesian Computation (ABC) procedure using a
stochastic model and a set of prior parameter draws.

## Usage

``` r
abc_stochastic_model(
  sm,
  priors,
  target,
  stat_func,
  states = NULL,
  tol = 0.1,
  method = "rejection",
  seed = 123
)

# S3 method for class 'abc_stochastic_model'
plot(x, ...)
```

## Arguments

- sm:

  An object of class `stochastic_model`.

- priors:

  A data frame or `tibble()` of prior parameter values. Each row should
  correspond to one simulation, and columns should match the parameter
  names expected by
  [`update_parameters()`](https://bccdc-phsa.github.io/brech/reference/stochastic_model.md).

- target:

  A named vector or data frame of observed summary statistics to match
  against simulated summary statistics.

- stat_func:

  A function that computes summary statistics from a simulation result.
  Should accept a model object and return a named numeric vector or
  `tibble()` row.

- states:

  Optional. A data frame or `tibble()` of initial states. Each row
  corresponds to a complete model state e.g. for an SIR model each row
  has three columns

- tol:

  A numeric tolerance (between 0 and 1) for the ABC rejection algorithm.
  Defaults to 0.1.

- method:

  A character string indicating the ABC method to use. Defaults to
  `"rejection"`. Passed directly to
  [`abc::abc()`](https://rdrr.io/pkg/abc/man/abc.html).

- seed:

  An integer specifying the random seed for reproducibility. Defaults to
  123.

- x:

  An object of class `abc_stochastic_model`.

- ...:

  list of additional plot arguments. `param` - if `type` is
  "prior_posterior" then parameter to plot `state` - if `type` is
  "simulations" then state to plot

## Value

An object of class `"abc_stochastic_model"` containing:

- model:

  The original stochastic model.

- fit:

  An object returned by
  [`abc::abc()`](https://rdrr.io/pkg/abc/man/abc.html) representing the
  ABC posterior.

- priors:

  The input prior values used for simulations.

- posteriors:

  The output posterior values used for simulations.

- prior_predictive:

  The summary statistics defined in `stat_func` drawn from the prior
  distribution.

- state:

  Final state of model.

A `ggplot` object

## Details

This function performs ABC by:

1.  Sampling parameter combinations from `priors`

2.  Updating the `stochastic_model` with each sampled parameter set

3.  Running the model and computing summary statistics via `stat_func`

4.  Comparing simulated statistics to `target` using the
    [`abc::abc()`](https://rdrr.io/pkg/abc/man/abc.html) function

Parallelization is handled using
[`furrr::future_map()`](https://furrr.futureverse.org/reference/future_map.html),
allowing reproducible simulation across parameter sets using a fixed
seed.

## Methods (by generic)

- `plot(abc_stochastic_model)`: Plot of either prior and posterior or
  projections

## Examples

``` r
reactions <- list(
  infection = list(
  transition = c("I" = +1),
  rate = function(x,p,t){ p$beta })
)
example_scenario <- list(
  params = list(beta = 0.1),
  initial_states = list(I = 0),
  sim_args = list(T = 10)
)
sm <- stochastic_model(reactions,example_scenario)
priors <- expand.grid(beta = seq(0.1, 0.5, length.out = 10))
target <- c(total_infected = 50)
abc_result <- abc_stochastic_model(
  sm,
  priors = priors,
  target = target,
  stat_func = function(sim) data.frame(total_infected = sim[nrow(sim), "I"]),
  tol = 0.2
)
```
