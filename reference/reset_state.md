# Convert cumulative state variables back to zero at a reference time

For cumulative state variables such as cases (e.g., `C_[0,1)`), this
function resets their values to 0 at a specified reference time (default
is 0). This is useful when wanting to compute cumulative changes from a
particular time point in a simulation or projection.

## Usage

``` r
reset_state(psm, state, reset_time = 0)
```

## Arguments

- psm:

  An object of class `projection_stochastic_model`, or a compatible data
  structure as accepted by `get_projection_dataframe()`.

- state:

  A character string specifying the prefix of the cumulative state
  variable to reset (e.g., `"C"` to match columns like `"C_[0,1)"`,
  `"C_[1,5)"`, etc.).

- reset_time:

  A numeric value indicating the time at which to reset the cumulative
  values back to zero. Defaults to 0.

## Value

A `tibble` with the same structure as the original projection data frame

## Examples

``` r
if (FALSE) { # \dontrun{
# Reset cumulative cases to 0 at time 10
reset_state(psm = fitted_model, state = "C", reset_time = 10)

# Reset cumulative hospitalizations
reset_state(psm = fitted_model, state = "H", reset_time = 5)
} # }
```
