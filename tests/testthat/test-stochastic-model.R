library(testthat)
library(purrr)
library(dplyr)

test_that("get_reactions_and_rates returns correct structure", {
  reactions <- list(
    infection = list(
      transition = c("I" = +1),
      rate = function(x, p, t) p$beta
    )
  )
  rr <- get_reactions_and_rates(reactions)
  expect_named(rr, c("rates", "transitions"))
  expect_type(rr$rates, "closure")
  expect_type(rr$transitions, "list")
  expect_equal(rr$transitions[[1]], c(I = +1))
})

test_that("validate_reactions catches invalid input", {
  expect_error(validate_reactions("not_a_list"), "must be a list")

  invalid_rxn <- list(foo = list())
  expect_error(validate_reactions(invalid_rxn), "must have a 'transition'")

  invalid_rxn <- list(foo = list(transition = c(I = 1, S = -1)))
  expect_error(validate_reactions(invalid_rxn), "must have a 'rate'")

  invalid_rxn <- list(foo = list(
    transition = c(x = 1),
    rate = 5
  ))
  expect_error(validate_reactions(invalid_rxn), "must be a function")

  invalid_rxn <- list(foo = list(
    transition = c(x = 1),
    rate = function(a, b) {}
  ))
  expect_error(validate_reactions(invalid_rxn), "must have arguments \\(x, p, t\\)")
})

test_that("validate_reactions passes valid input", {
  valid_rxn <- list(
    infection = list(
      transition = c("I" = +1),
      rate = function(x, p, t) p$beta
    )
  )
  expect_true(validate_reactions(valid_rxn))
})

test_that("validate_sim_scenario requires initial_states and params", {
  bad <- list(params = list())
  expect_error(validate_sim_scenario(bad), "sim_scenario must have")

  bad2 <- list(initial_states = list())
  expect_error(validate_sim_scenario(bad2), "sim_scenario must have")

  good <- list(initial_states = list(), params = list())
  expect_true(validate_sim_scenario(good))
})

test_that("stochastic_model returns classed object with expected structure", {
  reactions <- list(
    infection = list(
      transition = c("I" = +1),
      rate = function(x, p, t) p$beta
    )
  )
  sim <- list(
    initial_states = list(I = 0),
    params = list(beta = 0.2),
    sim_args = list(T = 5)
  )

  sm <- stochastic_model(reactions, sim)
  expect_s3_class(sm, "stochastic_model")
  expect_true("rates" %in% names(sm))
  expect_true("transitions" %in% names(sm))
  expect_equal(sm$params$beta, 0.2)
})

test_that("print.stochastic_model does not error", {
  reactions <- list(
    infection = list(
      transition = c("I" = +1),
      rate = function(x, p, t) p$beta
    )
  )
  sim <- list(
    initial_states = list(I = 0),
    params = list(beta = 0.2),
    sim_args = list(T = 5)
  )

  sm <- stochastic_model(reactions, sim)
  expect_invisible(print(sm))
})

test_that("run returns a data frame with time column", {

  reactions <- list(
    infection = list(
      transition = c("I" = +1),
      rate = function(x, p, t) p$beta
    )
  )
  sim <- list(
    initial_states = list(I = 0),
    params = list(beta = 0.1),
    sim_args = list(T = 2)
  )

  sm <- stochastic_model(reactions, sim)
  sim_out <- run_sim(sm)
  expect_true(is.matrix(sim_out))
  expect_true(all(c("time","I") %in% colnames(sim_out)))
})

test_that("interpolate_run_by_day interpolates correctly", {
  x <- cbind(
    time = c(0.5, 1.5, 2.5, 3.5),
    S = c(100, 90, 70, 50),
    I = c(0, 10, 30, 50)
  )

  out <- interpolate_run_by_day(x)
  expect_s3_class(out, "tbl_df")
  expect_true(all(c("time", "S", "I") %in% names(out)))
  expect_equal(out$time, 0:4)
  expect_equal(nrow(out), 5)
})

test_that("Adding params as a numeric vector returns an error", {

  reactions <- list(
    infection = list(
      transition = c("I" = +1),
      rate = function(x, p, t) p$beta
    )
  )
  sim <- list(
    initial_states = list(I = 0),
    params = c(beta = 0.1),
    sim_args = list(T = 2)
  )

  expect_error(stochastic_model(reactions, sim), "`params` must be a list")
})
