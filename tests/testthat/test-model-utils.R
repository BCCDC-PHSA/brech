library(testthat)
library(yaml)

test_that("load_model_params loads and transforms parameters correctly", {
  temp_yaml <- tempfile(fileext = ".yaml")
  writeLines("
params:
  latent_period: 2
  infectious_period: 4
  case_delay: 3
  r_0: 2
initial_states:
  S: [1000]
  V: [0]
", con = temp_yaml)

  result <- load_model_params(temp_yaml)

  expect_type(result, "list")
  expect_true("params" %in% names(result))
  expect_equal(result$params$theta, 0.5)
  expect_equal(result$params$gamma, 0.25)
  expect_equal(result$params$case_rate, 1/3)
  expect_true("C" %in% names(result$initial_states))
  expect_equal(result$initial_states$C, 0)
})

test_that("vaccinate_initial_conditions adjusts ungrouped states correctly", {
  sim <- list(initial_states = list(S = 80, V = 20))
  updated <- vaccinate_initial_conditions(sim, 0.75)
  expect_equal(updated$initial_states$V, 75)
  expect_equal(updated$initial_states$S, 25)
})

test_that("vaccinate_initial_conditions adjusts grouped states correctly", {
  sim <- list(initial_states = list(S_0 = 80, V_0 = 20, S_1 = 100, V_1 = 0))
  updated <- vaccinate_initial_conditions(sim, c("0" = 0.5, "1" = 0.25))
  expect_equal(updated$initial_states$V_0, 50)
  expect_equal(updated$initial_states$S_0, 50)
  expect_equal(updated$initial_states$V_1, 25)
  expect_equal(updated$initial_states$S_1, 75)
})

test_that("create_age_initial_conditions expands states by age group", {
  sim <- list(initial_states = list(S = c(10, 20), I = c(1, 2)))
  result <- create_age_initial_conditions(sim, c("[0-5)", "[5-19)"))
  expect_equal(result$initial_states[["S_[0-5)"]], 10)
  expect_equal(result$initial_states[["I_[5-19)"]], 2)
})

test_that("build_age_reactions returns expected transitions and rates", {
  cm <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
  colnames(cm) <- c("0-4", "5-19")
  reactions <- build_age_reactions(cm)

  expect_type(reactions, "list")
  expect_true("detectable_infection_0-4" %in% names(reactions))
  expect_true("recovery_5-19" %in% names(reactions))

  det_trans <- reactions[["detectable_infection_0-4"]]$transition
  expect_equal(det_trans[["S_0-4"]], -1)
  expect_equal(det_trans[["E_0-4"]], 1)

  # Test one rate function works with dummy data
  x <- c("I_0-4" = 10, "I_5-19" = 5, "S_0-4" = 90)
  p <- list(beta = 0.02, N = 2000, case_ascertainment = 1)
  t <- 0
  rate_fn <- reactions[["detectable_infection_0-4"]]$rate
  expect_type(rate_fn(x, p, t), "double")
})

test_that("calculate_age_structured_beta returns correct beta", {
  contact_matrix <- matrix(c(10, 5, 5, 20), ncol = 2)
  population_sizes <- c(1000, 2000)
  sm <- list(params = list(gamma = 0.25))
  R0 <- 2.5

  beta <- calculate_age_structured_beta(sm, contact_matrix, population_sizes, R0)
  expect_type(beta, "double")
  expect_gt(beta, 0)
})

