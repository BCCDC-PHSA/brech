# change vaccination coverage in initial conditions

change vaccination coverage in initial conditions

## Usage

``` r
vaccinate_initial_conditions(sim_scenario, vaccine_rate)
```

## Arguments

- sim_scenario:

  output of `load_model_params`

- vaccine_rate:

  numeric proportion vaccinated. Can be a numeric vector if using age
  groups. If named then vaccinates according to each age group with
  format "S\_{age_group}"
