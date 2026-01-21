# Create a set of scenarios that can be incorporated into [`projection_stochastic_model()`](https://bccdc-phsa.github.io/brech/reference/projection_stochastic_model.md)

Create a set of scenarios that can be incorporated into
[`projection_stochastic_model()`](https://bccdc-phsa.github.io/brech/reference/projection_stochastic_model.md)

## Usage

``` r
scenario_stochastic_model(sm, parameters = NULL, states = NULL)
```

## Arguments

- sm:

  Object of class `stochastic_model`.

- parameters:

  A data frame or tibble of prior parameter values. Each row should
  correspond to one simulation, and columns should match the parameter
  names expected by
  [`update_parameters()`](https://bccdc-phsa.github.io/brech/reference/stochastic_model.md).
  If `NULL` then no parameters will be updated an are all inherited from
  `sm`

- states:

  A data frame where each row represents an initial state. If `NULL`
  then this is generated from the initial states in `sm`

## See also

[`abc_stochastic_model()`](https://bccdc-phsa.github.io/brech/reference/abc_stochastic_model.md)
