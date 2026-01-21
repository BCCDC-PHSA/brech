# Plot Prior and Posterior Distributions with Optional True Value

This function generates a density plot comparing the prior and posterior
distributions of a parameter. An optional vertical dashed line can be
added to indicate the true parameter value.

## Usage

``` r
plot_prior_posterior(prior_values, posterior_values, true_value = NULL)
```

## Arguments

- prior_values:

  A numeric vector of samples from the prior distribution.

- posterior_values:

  A numeric vector of samples from the posterior distribution.

- true_value:

  Optional numeric value representing the true parameter value to be
  shown as a vertical dashed line.

## Value

A `ggplot2` object representing the density plot comparing the prior and
posterior distributions.

## Details

The function uses
[`geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html)
to visualize both prior and posterior samples, and
[`geom_vline()`](https://ggplot2.tidyverse.org/reference/geom_abline.html)
to optionally display the true parameter value. It uses
[`theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
and places the legend at the bottom for clarity.

## Examples

``` r
prior <- rnorm(1000, mean = 0, sd = 2)
posterior <- rnorm(1000, mean = 1, sd = 0.5)
true_value <- 1.2
plot_prior_posterior(prior, posterior, true_value)

```
