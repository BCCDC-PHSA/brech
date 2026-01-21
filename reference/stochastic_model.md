# update parameters for a stochastic_model

update parameters for a stochastic_model

update initial state for a stochastic_model

Create a stochastic model

run stochastic model using adaptive tau set algorithm

## Usage

``` r
update_parameters(x, ..., params = list())

update_state(x, state)

stochastic_model(reactions, sim_scenario)

run_sim(x, tf = NULL)
```

## Arguments

- x:

  An object of class `stochastic_model`.

- ...:

  parameters to update

- params:

  parameters to update passed in as a list

- state:

  data.frame

- reactions:

  list of state transitions. Each reaction is itself a list that
  requires a `transition` vector and a `rate` vector

- sim_scenario:

  list that defines the scenario including `params`, `initial_states`
  and `sim_args` that must at least include simulation time `T`

- tf:

  maximum time

## Value

An object of class `stochastic_model`.

## Examples

``` r
if (FALSE) { # \dontrun{
sm <- update_parameters(sm,beta=0.2)
} # }
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
```
