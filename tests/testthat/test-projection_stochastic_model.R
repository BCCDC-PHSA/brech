library(dplyr)

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

abc_fit <- abc_stochastic_model(sm, priors, target, stat_func, tol = 0.5)

psm <- projection_stochastic_model(abc_fit, project_time = 5, seed = 42)

test_that("projection_stochastic_model works correctly with abc_stochastic_model", {

  expect_s3_class(psm, "projection_stochastic_model")
  expect_true("projection" %in% names(psm))
  expect_true("model" %in% names(psm))
  expect_gt(nrow(psm$projection), 0)
})

test_that("add_projection_date works correctly", {
  # Reuse the previous projection
  start_date <- "2020-01-01"
  psm <- add_projection_date(psm, start_date)
  expect_true("date" %in% names(psm$projection))
  expect_equal(min(psm$projection$date), as.Date(start_date))
})

test_that("collapse_states collapses prefixed states correctly", {
  proj <- list("projection" = tibble::tibble(
    time = 0:2,
    I_1 = c(1, 2, 3),
    I_2 = c(2, 3, 4)
  ))
  class(proj) <- c("projection_stochastic_model", class(proj))
  result <- collapse_states(proj)
  expect_true("I" %in% names(result))
  expect_equal(result$I, c(3, 5, 7))
})

test_that("create_age_group_column works with proper format", {
  proj <- tibble::tibble(
    time = 0:1,
    C_0_4 = c(1, 2),
    C_5_17 = c(3, 4)
  )
  res <- create_age_group_column(proj, state = "C")
  expect_true("age_group" %in% names(res))
  expect_true("C" %in% names(res))
})

test_that("difference_of_states returns differences by group", {
  proj <- list("projection" = tibble::tibble(
    time = c(0, 1, 2, 0, 1, 2),
    C_0_4 = c(0,1, 2, 0, 1, 2),
    C_5_17 = c(0, 2, 4,  0, 2, 4)
  ))
  class(proj) <- c("projection_stochastic_model", class(proj))
  diff <- difference_of_states(proj, "C")
  expect_equal(diff$C_0_4, c(0, 1, 1, 0, 1, 1))
  expect_equal(diff$C_5_17, c(0, 2, 2, 0, 2, 2))
})

