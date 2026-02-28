## Tests for competing risks survival forests

library(ranger)
library(survival)
context("ranger_competing_risks")

## Create test data with competing risks
set.seed(42)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
time <- rexp(n, rate = exp(0.3 * x1))
event_type <- sample(0:2, n, replace = TRUE, prob = c(0.3, 0.4, 0.3))
dat_cr <- data.frame(time = time, status = event_type, x1 = x1, x2 = x2)

## Standard single-event survival data for backward compatibility
dat_surv <- dat_cr
dat_surv$status <- ifelse(dat_cr$status > 0, 1, 0)

test_that("competing risks forest grows without error", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20)
  expect_is(rf, "ranger")
  expect_equal(rf$treetype, "Survival")
})

test_that("correct number of event types detected", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20)
  expect_equal(rf$num.event.types, 2)
})

test_that("CHF output is list of matrices for competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20)
  expect_is(rf$chf, "list")
  expect_equal(length(rf$chf), 2)
  expect_is(rf$chf[[1]], "matrix")
  expect_is(rf$chf[[2]], "matrix")
  expect_equal(nrow(rf$chf[[1]]), n)
  expect_equal(ncol(rf$chf[[1]]), length(rf$unique.death.times))
})

test_that("CIF output is computed and valid", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20)
  expect_is(rf$cif, "list")
  expect_equal(length(rf$cif), 2)
  expect_equal(nrow(rf$cif[[1]]), n)

  # CIF values should be between 0 and 1
  for (e in 1:2) {
    expect_true(all(rf$cif[[e]] >= -1e-10))
    expect_true(all(rf$cif[[e]] <= 1 + 1e-10))
  }

  # Sum of CIFs should be <= 1 at all timepoints
  cif_sum <- rf$cif[[1]] + rf$cif[[2]]
  expect_true(all(cif_sum <= 1 + 1e-10))
})

test_that("predictions work for new data with competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20)
  pred <- predict(rf, data = dat_cr[1:10, ])
  expect_is(pred$chf, "list")
  expect_equal(length(pred$chf), 2)
  expect_equal(nrow(pred$chf[[1]]), 10)
  expect_is(pred$cif, "list")
  expect_equal(length(pred$cif), 2)
})

test_that("backward compatibility: single-event survival identical behavior", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_surv, num.trees = 20, seed = 123)
  expect_is(rf$chf, "matrix")
  expect_is(rf$survival, "matrix")
  expect_null(rf$cif)
  expect_equal(rf$num.event.types, 1)
})

test_that("formula interface with Surv() works for competing risks", {
  rf <- ranger(Surv(time, status) ~ x1 + x2, data = dat_cr, num.trees = 20)
  expect_equal(rf$num.event.types, 2)
})

test_that("alternative interface works for competing risks", {
  rf <- ranger(dependent.variable.name = "time", status.variable.name = "status",
               data = dat_cr, num.trees = 20)
  expect_equal(rf$num.event.types, 2)
  expect_is(rf$chf, "list")
})

test_that("OOB error is computed for competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 50)
  expect_true(!is.null(rf$prediction.error))
  expect_true(rf$prediction.error >= 0 && rf$prediction.error <= 1)
})

test_that("extratrees splitrule works with competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20,
               splitrule = "extratrees")
  expect_is(rf, "ranger")
  expect_equal(rf$num.event.types, 2)
})

test_that("maxstat splitrule works with competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20,
               splitrule = "maxstat")
  expect_is(rf, "ranger")
  expect_equal(rf$num.event.types, 2)
})

test_that("C-index splitrule errors for competing risks", {
  expect_error(
    ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20, splitrule = "C"),
    "C-index.*splitting not supported for competing risks"
  )
})

test_that("variable importance works with competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 50,
               importance = "permutation")
  expect_true(!is.null(rf$variable.importance))
  expect_equal(length(rf$variable.importance), 2)
})

test_that("impurity importance works with competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 50,
               importance = "impurity")
  expect_true(!is.null(rf$variable.importance))
  expect_equal(length(rf$variable.importance), 2)
})

test_that("invalid status values give error", {
  dat_bad <- dat_cr
  dat_bad$status[1] <- -1
  expect_error(
    ranger(dependent.variable.name = "time", status.variable.name = "status",
           data = dat_bad, num.trees = 20),
    "consecutive integers"
  )
})

test_that("three event types work", {
  dat3 <- dat_cr
  dat3$status <- sample(0:3, n, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
  rf <- ranger(Surv(time, status) ~ ., data = dat3, num.trees = 20)
  expect_equal(rf$num.event.types, 3)
  expect_equal(length(rf$chf), 3)
  expect_equal(length(rf$cif), 3)
})

test_that("forest object stores num.event.types", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20)
  expect_equal(rf$forest$num.event.types, 2)
})

test_that("print works for competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 20)
  expect_output(print(rf), "Number of event types")
})

test_that("CHF values are non-decreasing for competing risks", {
  rf <- ranger(Surv(time, status) ~ ., data = dat_cr, num.trees = 50)
  for (e in 1:2) {
    for (i in 1:min(10, nrow(rf$chf[[e]]))) {
      diffs <- diff(rf$chf[[e]][i, ])
      expect_true(all(diffs >= -1e-10))
    }
  }
})
