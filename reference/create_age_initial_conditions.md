# Convert a vector of each state into a labeled state using `age_groups` labels

Convert a vector of each state into a labeled state using `age_groups`
labels

## Usage

``` r
create_age_initial_conditions(sim_scenario, age_groups)
```

## Arguments

- sim_scenario:

  output of `load_model_params`

- age_groups:

  vector of age group labels

## Examples

``` r
sim_scenario <- list(
  "initial_states" = list(
    "S" = c(1,2,3),
    "I" = c(0,0,1)
  )
)
age_groups <- c("0-5","5-19","19+")
create_age_initial_conditions(sim_scenario,age_groups)
#> $initial_states
#> $initial_states$`S_0-5`
#> [1] 1
#> 
#> $initial_states$`S_5-19`
#> [1] 2
#> 
#> $initial_states$`S_19+`
#> [1] 3
#> 
#> $initial_states$`I_0-5`
#> [1] 0
#> 
#> $initial_states$`I_5-19`
#> [1] 0
#> 
#> $initial_states$`I_19+`
#> [1] 1
#> 
#> 
```
