# Temporary test scaffolding for Graf score splitting
library(ranger)
library(survival)
context("ranger_surv_graf")

# Extremely rudimentary "it does not explode" test
test_that("Graf splitting is callable", {

  expect_failure(expect_error(
    ranger(
      Surv(time, status) ~ .,
      data = veteran,
      splitrule = "brier",
      num.trees = 10
    )
  ))
})
