# Create a stochastic model

This function provides a formatted and colorized summary of a stochastic
model, including initial values, parameters, and a list of reactions
with their transitions and rate function arguments.

## Usage

``` r
stochastic_model(reactions, sim_scenario)

# S3 method for class 'stochastic_model'
print(x, ...)

run_sim(x, tf = NULL)

update_parameters(x, ..., params = list())

update_state(x, state)

interpolate_run_by_day(x)
```

## Arguments

- reactions:

  list of state transitions. Each reaction is itself a list that
  requires a `transition` vector and a `rate` vector

- sim_scenario:

  list that defines the scenario including `params`, `initial_states`
  and `sim_args` that must at least include simulation time `T`

- x:

  A matrix or data frame with a numeric column named "time" and one or
  more state columns with arbitrary names. output of `run`

- ...:

  parameters to update

- tf:

  maximum time

- params:

  parameters to update passed in as a list

- state:

  data.frame

## Value

Invisibly returns the `x` object.

An object of class `stochastic_model`.

An object of class `stochastic_model`.

A `tibble()` with one row per whole day and interpolated values for each
state.

## Methods (by generic)

- `print(stochastic_model)`: Print method for stochastic_model objects

## Functions

- `run_sim()`: run stochastic model using adaptive tau set algorithm

- `update_parameters()`: update parameters for a stochastic_model

- `update_state()`: update initial state for a stochastic_model

- `interpolate_run_by_day()`: Interpolate State Trajectories to Whole
  Days

  Converts a time-series data frame or matrix with continuous time
  points and one or more state columns into a `tibble()` with daily
  interpolated values.

## Examples

``` r
reactions <- list(
  infection = list(
    transition = c("I" = +1),
    rate = function(x,p,t){p$beta})
)
example_scenario <- list(
  params = list(beta = 0.1),
  initial_states = list(I = 0),
  sim_args = list(T = 10)
)
stochastic_model(reactions,example_scenario)
#> Stochastic Model 
#> ================== 
#> 
#> Initial Values:
#>   I : 0
#> Parameters:
#>   beta : 0.1
#> Simulation arguments:
#>   T : 10
#> Reactions:
#>  - infection:
#>    Transition: [I +1]
#> 

reactions <- list(
  infection = list(
    transition = c("I" = +1),
    rate = function(x,p,t){p$beta})
)
example_scenario <- list(
  params = list(beta = 0.1),
  initial_states = list(I = 0),
  sim_args = list(T = 10)
)
sm <- stochastic_model(reactions,example_scenario)
sm <- update_parameters(sm,beta=0.2)

mat <- cbind(
  time = c(0.1, 0.5, 1.2, 2.4, 3.8),
  S = c(1000, 900, 800, 600, 500),
  I = c(1, 10, 50, 100, 90)
)
interpolate_run_by_day(mat)
#> # A tibble: 5 × 3
#>    time     S     I
#>   <int> <dbl> <dbl>
#> 1     0   NA   NA  
#> 2     1  829.  38.6
#> 3     2  667.  83.3
#> 4     3  557.  95.7
#> 5     4   NA   NA  
```
