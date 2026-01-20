library(testthat)
library(dplyr)

test_that("scenario_stochastic_model returns correct structure", {
  sm <- list(
    initial_states = list(S = 100, I = 1),
    params = list(beta = 0.2, gamma = 0.1)
  )
  class(sm) <- "stochastic_model"

  parameters <- tibble(beta = c(0.1, 0.2), gamma = c(0.1, 0.2))
  states <- tibble(S = c(100, 80), I = c(1, 5))

  scen <- scenario_stochastic_model(sm, parameters, states)

  expect_s3_class(scen, "scenario_stochastic_model")
  expect_equal(nrow(scen$parameters), 2)
  expect_equal(nrow(scen$states), 2)
  expect_equal(scen$model, sm)
})

test_that("scenario_stochastic_model replicates states when only parameters are provided", {
  sm <- list(
    initial_states = list(S = 50, I = 2),
    params = list(beta = 0.1)
  )
  class(sm) <- "stochastic_model"

  parameters <- tibble(beta = c(0.1, 0.2, 0.3))

  scen <- scenario_stochastic_model(sm, parameters = parameters)

  expect_equal(nrow(scen$parameters), 3)
  expect_equal(nrow(scen$states), 3)
  expect_equal(unique(scen$states$S), 50)
  expect_equal(unique(scen$states$I), 2)
})

test_that("scenario_stochastic_model replicates parameters when only states are provided", {
  sm <- list(
    initial_states = list(S = 60, I = 4),
    params = list(beta = 0.5, gamma = 0.1)
  )
  class(sm) <- "stochastic_model"

  states <- tibble(S = c(60, 70, 80), I = c(4, 5, 6))

  scen <- scenario_stochastic_model(sm, states = states)

  expect_equal(nrow(scen$states), 3)
  expect_equal(nrow(scen$parameters), 3)
  expect_equal(unique(scen$parameters[,"beta"]), 0.5)
})

test_that("scenario_stochastic_model errors if neither parameters nor states are provided", {
  sm <- list(
    initial_states = list(S = 100, I = 1),
    params = list(beta = 0.2)
  )
  class(sm) <- "stochastic_model"

  expect_error(scenario_stochastic_model(sm),
               "Need to define at least one of `parameters` or `states`")
})

test_that("scenario_stochastic_model errors if parameters and states have mismatched rows", {
  sm <- list(
    initial_states = list(S = 100, I = 1),
    params = list(beta = 0.2)
  )
  class(sm) <- "stochastic_model"

  parameters <- tibble(beta = c(0.1, 0.2))
  states <- tibble(S = c(100), I = c(1))

  expect_error(
    scenario_stochastic_model(sm, parameters, states),
    "`parameters` and `states` need the same number of rows"
  )
})
