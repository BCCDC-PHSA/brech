#' Sample from a Truncated Normal Distribution
#'
#' Generates random samples from a normal distribution truncated between
#' specified lower and upper bounds.
#'
#' @param n Integer. Number of samples to draw.
#' @param mu Numeric. Mean of the normal distribution.
#' @param sigma Numeric. Standard deviation of the normal distribution.
#' @param low Numeric. Lower truncation bound.
#' @param high Numeric. Upper truncation bound.
#'
#' @return A numeric vector of length \code{n} containing random draws from the truncated normal distribution.
#'
#' @details The function uses inverse transform sampling by first computing the
#' cumulative probabilities corresponding to the lower and upper bounds, then
#' sampling uniformly from this range and transforming back via the normal quantile function.
#'
#' @examples
#' set.seed(123)
#' rtruncnorm(n = 5, mean = 0, sd = 1, low = -1, high = 1)
#'
#' @export
rtruncnorm <- function(n, mean, sd, low, high) {

  p_low <- stats::pnorm(low, mean, sd)
  p_high <- stats::pnorm(high, mean, sd)


  stats::qnorm(stats::runif(n, p_low, p_high), mean, sd)
}
