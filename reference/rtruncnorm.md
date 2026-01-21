# Sample from a Truncated Normal Distribution

Generates random samples from a normal distribution truncated between
specified lower and upper bounds.

## Usage

``` r
rtruncnorm(n, mean, sd, low, high)
```

## Arguments

- n:

  Integer. Number of samples to draw.

- mean:

  Numeric. Mean of the normal distribution.

- sd:

  Numeric. Standard deviation of the normal distribution.

- low:

  Numeric. Lower truncation bound.

- high:

  Numeric. Upper truncation bound.

## Value

A numeric vector of length `n` containing random draws from the
truncated normal distribution.

## Details

The function uses inverse transform sampling by first computing the
cumulative probabilities corresponding to the lower and upper bounds,
then sampling uniformly from this range and transforming back via the
normal quantile function.

## Examples

``` r
set.seed(123)
rtruncnorm(n = 5, mean = 0, sd = 1, low = -1, high = 1)
#> [1] -0.3719060  0.5152845 -0.1563984  0.7110777  0.8441328
```
