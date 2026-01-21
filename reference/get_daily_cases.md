# Compute Daily Incident Cases from Cumulative Simulation Output

This function takes a data frame or matrix containing irregularly spaced
simulation time points and cumulative case counts (`C`), and returns the
number of new cases per day using linear interpolation.

## Usage

``` r
get_daily_cases(r)
```

## Arguments

- r:

  A data frame or matrix with columns `"time"` and `"C"`, representing
  time (in continuous units) and cumulative number of cases,
  respectively.

## Value

A tibble with two columns:

- day:

  Integer day values from the floor of the minimum time to the ceiling
  of the maximum time

- daily_incidence:

  Number of new cases per day, computed as the difference in
  interpolated cumulative cases

## Details

Since simulation outputs may record cumulative cases at non-integer time
points, this function performs linear interpolation to estimate
cumulative cases at each whole-number day. It then computes the number
of incident cases as the first difference of the interpolated values.
The first day is assigned 0 new cases by default.

## Examples

``` r
sim_data <- as.matrix(dplyr::tibble(
  time = c(0, 0.1, 0.9, 1.4, 2.3, 3.8),
  C = c(0, 1, 5, 7, 10, 15)
))
get_daily_cases(sim_data)
#> # A tibble: 5 × 2
#>     day daily_incidence
#>   <int>           <dbl>
#> 1     0            0   
#> 2     1            5.4 
#> 3     2            3.6 
#> 4     3            3.33
#> 5     4           NA   
```
