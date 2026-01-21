# Interpolate State Trajectories to Whole Days

Converts a time-series data frame or matrix with continuous time points
and one or more state columns into a tibble with daily interpolated
values.

## Usage

``` r
interpolate_run_by_day(x)
```

## Arguments

- x:

  A matrix or data frame with a numeric column named "time" and one or
  more state columns with arbitrary names. output of `run`

## Value

A tibble with one row per whole day and interpolated values for each
state.

## Examples

``` r
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
