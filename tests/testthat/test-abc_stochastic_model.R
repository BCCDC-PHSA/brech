library(testthat)
library(dplyr)
library(ggplot2)

test_that("update_parameters updates model parameters", {
  model <- list(params = list(beta = 0.1))
  class(model) <- "stochastic_model"

  updated <- update_parameters(model, beta = 0.2, gamma = 0.3)
  expect_equal(updated$params$beta, 0.2)
  expect_equal(updated$params$gamma, 0.3)

  updated2 <- update_parameters(model, params = list(beta = 0.5))
  expect_equal(updated2$params$beta, 0.5)
})

test_that("update_state replaces initial states", {
  model <- list(initial_states = list(S = 100, I = 1))
  class(model) <- "stochastic_model"

  new_state <- data.frame(S = 80, I = 20)
  updated <- update_state(model, new_state)

  expect_named(updated$initial_states, c("S", "I"))
  expect_equal(updated$initial_states[["S"]], 80)
  expect_equal(updated$initial_states[["I"]], 20)
})

test_that("get_daily_cases returns daily incidence correctly", {
  sim_data <- as.matrix(tibble(
    time = seq(from=0,to=10,length.out=100),
    C =    seq(from=0,to=100,length.out=100)
  ))
  out <- get_daily_cases(sim_data)
  expect_s3_class(out, "tbl_df")
  expect_equal(names(out), c("day", "daily_incidence"))
  expect_equal(nrow(out), 11)  # days 0 to 10
})

test_that("plot_prior_posterior creates plot object", {
  prior <- rnorm(1000, mean = 0, sd = 1)
  posterior <- rnorm(1000, mean = 1, sd = 0.5)

  g <- plot_prior_posterior(prior, posterior, true_value = 1.2)
  expect_s3_class(g, "ggplot")
})

test_that("abc_stochastic_model works on small input", {
  skip_if_not_installed("abc")
  skip_if_not_installed("furrr")
  skip_if_not_installed("adaptivetau")
  skip_if_not_installed("progressr")
  set.seed(123)

  reactions <- list(
    infection = list(
      transition = c("I" = +1),
      rate = function(x, p, t) p$beta
    )
  )
  sm <- stochastic_model(reactions, list(
    initial_states = list(I = 0),
    params = list(beta = 0.1),
    sim_args = list(T = 1)
  ))

  priors <- tibble(beta = pmax(0,rnorm(100,mean=0.1,sd=0.05)))
  target <- 1
  stat_func <- function(sim) sim[nrow(sim), "I"]

  result <- abc_stochastic_model(sm, priors, target, stat_func, tol = 1)

  expect_s3_class(result, "abc_stochastic_model")
  expect_named(result, c("model", "fit", "priors", "posteriors", "prior_predictive", "state"))
})

test_that("plot.abc_stochastic_model generates prior/posterior plot", {
  prior <- rnorm(1000, 0, 1)
  posterior <- rnorm(1000, 1, 0.5)

  fake_model <- structure(list(
    params = list(beta = 0.5)
  ), class = "stochastic_model")

  abc_model <- structure(list(
    model = fake_model,
    priors = tibble(beta = prior),
    posteriors = tibble(beta = posterior)
  ), class = "abc_stochastic_model")

  g <- plot(abc_model, param = "beta")
  expect_s3_class(g, "ggplot")
})
