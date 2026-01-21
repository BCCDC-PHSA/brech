# Using the next generation matrix approach, calculates beta based on a fixed R0

Using the next generation matrix approach, calculates beta based on a
fixed R0

## Usage

``` r
calculate_age_structured_beta(sm, contact_matrix, population_sizes, R0)
```

## Arguments

- sm:

  Object of class `stochastic_model`

- contact_matrix:

  An m x m contact matrix

- population_sizes:

  An m vector of population size

- R0:

  numeric reproduction number

## Value

numeric
