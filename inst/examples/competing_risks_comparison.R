## Competing risks: ranger vs randomForestSRC comparison
##
## Both packages are fit with cause-specific log-rank splitting
## (splitrule = "logrank") and otherwise similar settings.
## The resulting CIF curves should look similar for a large number of trees.

library(ranger)
library(randomForestSRC)
library(survival)

## ---------------------------------------------------------------------------
## 1. Simulate competing risks data
## ---------------------------------------------------------------------------
set.seed(42)
n <- 800

x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rbinom(n, 1, 0.5)

## Latent event times for two competing causes
t1 <- rexp(n, rate = exp(0.5 * x1 + 0.2 * x3))   # cause 1
t2 <- rexp(n, rate = exp(0.3 * x2 - 0.1 * x3))   # cause 2
tc <- rexp(n, rate = 0.3)                          # censoring

## Observed time and status
time   <- pmin(t1, t2, tc)
status <- ifelse(time == tc, 0L,
          ifelse(time == t1, 1L, 2L))

dat <- data.frame(time = time, status = status, x1 = x1, x2 = x2, x3 = x3)

cat("Event distribution:\n")
print(table(status = dat$status))

## ---------------------------------------------------------------------------
## 2. Fit forests with matched settings
## ---------------------------------------------------------------------------

## Common settings
ntree    <- 1000
nodesize <- 15
mtry     <- 2

## ranger: uses composite cause-specific log-rank by default
## replace = TRUE is the default (bootstrap with replacement)
t_ranger <- system.time(
  rf_ranger <- ranger(
    Surv(time, status) ~ .,
    data      = dat,
    num.trees = ntree,
    min.node.size = nodesize,
    mtry      = mtry,
    replace   = TRUE,
    seed      = 1
  )
)

## randomForestSRC: use splitrule = "logrank" for cause-specific log-rank
## (the default "logrankCR" uses Gray's test, which differs)
## Match ranger defaults: bootstrap w/ replacement, all split points, all timepoints
t_rfsrc <- system.time(
  rf_rfsrc <- rfsrc(
    Surv(time, status) ~ .,
    data      = dat,
    ntree     = ntree,
    nodesize  = nodesize,
    mtry      = mtry,
    splitrule = "logrank",
    nsplit    = 0,
    samptype  = "swr",
    ntime     = 0,
    seed      = -1
  )
)

cat("\nTiming:\n")
cat("  ranger:          ", round(t_ranger["elapsed"], 2), "s\n")
cat("  randomForestSRC: ", round(t_rfsrc["elapsed"], 2), "s\n")

## ---------------------------------------------------------------------------
## 3. Extract CIF predictions (OOB)
## ---------------------------------------------------------------------------

## ranger: cif is a list of K matrices (n x num_timepoints)
cif_ranger_1 <- rf_ranger$cif[[1]]
cif_ranger_2 <- rf_ranger$cif[[2]]
times_ranger <- rf_ranger$unique.death.times

## randomForestSRC: cif.oob is a 3d array (n x num_timepoints x K)
cif_rfsrc_1 <- rf_rfsrc$cif.oob[, , 1]
cif_rfsrc_2 <- rf_rfsrc$cif.oob[, , 2]
times_rfsrc <- rf_rfsrc$time.interest

## ---------------------------------------------------------------------------
## 4. Compare mean CIF curves across all subjects
## ---------------------------------------------------------------------------

mean_cif_ranger_1 <- colMeans(cif_ranger_1)
mean_cif_ranger_2 <- colMeans(cif_ranger_2)

mean_cif_rfsrc_1 <- colMeans(cif_rfsrc_1)
mean_cif_rfsrc_2 <- colMeans(cif_rfsrc_2)

## ---------------------------------------------------------------------------
## 5. Plot mean CIF curves
## ---------------------------------------------------------------------------

pdf("competing_risks_comparison.pdf", width = 12, height = 10)

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

## Cause 1
plot(times_ranger, mean_cif_ranger_1, type = "s", col = "steelblue", lwd = 2,
     xlab = "Time", ylab = "Mean CIF", main = "Cause 1 -Mean CIF",
     ylim = c(0, max(mean_cif_ranger_1, mean_cif_rfsrc_1) * 1.1))
lines(times_rfsrc, mean_cif_rfsrc_1, type = "s", col = "firebrick", lwd = 2, lty = 2)
legend("bottomright", legend = c("ranger", "randomForestSRC"),
       col = c("steelblue", "firebrick"), lwd = 2, lty = c(1, 2), bty = "n")

## Cause 2
plot(times_ranger, mean_cif_ranger_2, type = "s", col = "steelblue", lwd = 2,
     xlab = "Time", ylab = "Mean CIF", main = "Cause 2 -Mean CIF",
     ylim = c(0, max(mean_cif_ranger_2, mean_cif_rfsrc_2) * 1.1))
lines(times_rfsrc, mean_cif_rfsrc_2, type = "s", col = "firebrick", lwd = 2, lty = 2)
legend("bottomright", legend = c("ranger", "randomForestSRC"),
       col = c("steelblue", "firebrick"), lwd = 2, lty = c(1, 2), bty = "n")

## ---------------------------------------------------------------------------
## 6. Compare CIF for individual subjects (high/low risk for cause 1)
## ---------------------------------------------------------------------------

## Pick a high-risk (large x1) and low-risk (small x1) subject
i_high <- which.max(dat$x1)
i_low  <- which.min(dat$x1)

## High-risk subject
plot(times_ranger, cif_ranger_1[i_high, ], type = "s", col = "steelblue", lwd = 2,
     xlab = "Time", ylab = "CIF (Cause 1)",
     main = sprintf("High-risk subject (x1 = %.1f)", dat$x1[i_high]),
     ylim = c(0, 1))
lines(times_rfsrc, cif_rfsrc_1[i_high, ], type = "s", col = "firebrick", lwd = 2, lty = 2)
legend("bottomright", legend = c("ranger", "randomForestSRC"),
       col = c("steelblue", "firebrick"), lwd = 2, lty = c(1, 2), bty = "n")

## Low-risk subject
plot(times_ranger, cif_ranger_1[i_low, ], type = "s", col = "steelblue", lwd = 2,
     xlab = "Time", ylab = "CIF (Cause 1)",
     main = sprintf("Low-risk subject (x1 = %.1f)", dat$x1[i_low]),
     ylim = c(0, 1))
lines(times_rfsrc, cif_rfsrc_1[i_low, ], type = "s", col = "firebrick", lwd = 2, lty = 2)
legend("bottomright", legend = c("ranger", "randomForestSRC"),
       col = c("steelblue", "firebrick"), lwd = 2, lty = c(1, 2), bty = "n")

## ---------------------------------------------------------------------------
## 7. Prediction on new data
## ---------------------------------------------------------------------------

newdat <- data.frame(
  x1 = c(-1, 0, 1),
  x2 = c( 0, 0, 0),
  x3 = c( 0, 1, 1)
)

pred_ranger <- predict(rf_ranger, data = newdat)
pred_rfsrc  <- predict(rf_rfsrc, newdata = newdat)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (i in 1:3) {
  plot(pred_ranger$unique.death.times, pred_ranger$cif[[1]][i, ],
       type = "s", col = "steelblue", lwd = 2,
       xlab = "Time", ylab = "CIF (Cause 1)",
       main = sprintf("x1=%.0f, x2=%.0f, x3=%.0f", newdat$x1[i], newdat$x2[i], newdat$x3[i]),
       ylim = c(0, 1))
  lines(pred_rfsrc$time.interest, pred_rfsrc$cif[i, , 1],
        type = "s", col = "firebrick", lwd = 2, lty = 2)
  legend("bottomright", legend = c("ranger", "rfsrc"),
         col = c("steelblue", "firebrick"), lwd = 2, lty = c(1, 2), bty = "n")
}

dev.off()

## ---------------------------------------------------------------------------
## 8. Numeric comparison
## ---------------------------------------------------------------------------

## Interpolate rfsrc CIF onto ranger's time grid for direct comparison
interp <- function(x_from, y_from, x_to) {
  stepfun(x_from, c(0, y_from))(x_to)
}

cif1_rfsrc_interp <- t(sapply(1:n, function(i) {
  interp(times_rfsrc, cif_rfsrc_1[i, ], times_ranger)
}))

## Per-subject mean absolute difference
mad_per_subject <- rowMeans(abs(cif_ranger_1 - cif1_rfsrc_interp))

cat("\nCause-1 CIF comparison (ranger vs rfsrc):\n")
cat("  Per-subject mean absolute difference:\n")
cat("    median: ", round(median(mad_per_subject), 4), "\n")
cat("    95th:   ", round(quantile(mad_per_subject, 0.95), 4), "\n")
cat("    max:    ", round(max(mad_per_subject), 4), "\n")
cat("  Correlation of subject-level CIF at final timepoint:\n")
cat("   ", round(cor(cif_ranger_1[, ncol(cif_ranger_1)],
                     cif1_rfsrc_interp[, ncol(cif1_rfsrc_interp)]), 4), "\n")

cat("\nPlots saved to competing_risks_comparison.pdf\n")
