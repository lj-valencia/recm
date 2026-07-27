# Unit tests for predict.recmfit() in R/predict.R

test_that("in-sample prediction matches fitted values", {
  set.seed(2026)
  n <- 120
  dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
  dat <- data.frame(ystar = cumsum(dys))
  dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))

  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
                       var_lags = 2, restarts = 1, quiet = TRUE)

  preds <- predict(fit, type = "diff")
  expect_equal(length(preds), nrow(dat))
  # Non-NA predictions on index match object$fitted exactly
  expect_equal(preds[fit$index], fit$fitted)
})

test_that("predict level prediction equals lag(y) + diff prediction", {
  set.seed(2026)
  n <- 120
  dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
  dat <- data.frame(ystar = cumsum(dys))
  dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))

  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
                       var_lags = 2, restarts = 1, quiet = TRUE)

  pred_diff <- predict(fit, type = "diff")
  pred_level <- predict(fit, type = "level")

  lag_y <- c(NA_real_, dat$y[-nrow(dat)])
  expect_equal(pred_level, lag_y + pred_diff)
})

test_that("out-of-sample predictions work on newdata", {
  set.seed(2026)
  n <- 150
  dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
  dat <- data.frame(ystar = cumsum(dys))
  dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))

  train_dat <- dat[1:120, ]
  test_dat  <- dat[115:150, ]

  fit <- recm_estimate("y", "ystar", beta = 0.995, data = train_dat, m = 2,
                       var_lags = 2, restarts = 1, quiet = TRUE)

  preds_out <- predict(fit, newdata = test_dat, type = "diff")
  expect_equal(length(preds_out), nrow(test_dat))
  expect_true(any(is.finite(preds_out)))
})

test_that("confidence and prediction intervals without boot_object return correct structure", {
  set.seed(2026)
  n <- 120
  dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
  dat <- data.frame(ystar = cumsum(dys))
  dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))

  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
                       var_lags = 2, restarts = 1, quiet = TRUE)

  ci <- predict(fit, interval = "confidence", quiet = TRUE)
  pi <- predict(fit, interval = "prediction", quiet = TRUE)

  expect_s3_class(ci, "data.frame")
  expect_named(ci, c("fit", "se.fit", "lwr", "upr"))
  expect_named(pi, c("fit", "se.fit", "lwr", "upr"))

  ok <- fit$index
  expect_true(all(ci$lwr[ok] <= ci$fit[ok]))
  expect_true(all(ci$upr[ok] >= ci$fit[ok]))

  # Prediction intervals should be wider than confidence intervals
  expect_true(mean(pi$upr[ok] - pi$lwr[ok]) > mean(ci$upr[ok] - ci$lwr[ok]))
})

test_that("bootstrap integrated intervals work with recm_boot()", {
  set.seed(2026)
  n <- 100
  dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
  dat <- data.frame(ystar = cumsum(dys))
  dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))

  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
                       var_lags = 2, restarts = 0, quiet = TRUE)

  bs <- recm_boot(fit, R = 15, seed = 2026)
  ci_boot <- predict(fit, interval = "confidence", boot_object = bs)

  expect_s3_class(ci_boot, "data.frame")
  expect_named(ci_boot, c("fit", "se.fit", "lwr", "upr"))

  ok <- fit$index
  expect_true(all(ci_boot$lwr[ok] <= ci_boot$fit[ok]))
  expect_true(all(ci_boot$upr[ok] >= ci_boot$fit[ok]))
})

test_that("predict errors on missing columns in newdata", {
  set.seed(2026)
  n <- 100
  dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
  dat <- data.frame(ystar = cumsum(dys))
  dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))

  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
                       var_lags = 2, restarts = 0, quiet = TRUE)

  bad_data <- data.frame(wrong_col = 1:10)
  expect_error(predict(fit, newdata = bad_data), "column 'y'")
})


# ---- the construction, checked against an independent route ----
#
# predict() reaches the forward sum through the expectations mechanism. The
# tests below write the construction out by hand instead -- the fitted
# alpha, the FITTED H, states built from the new observations, lagged one
# period -- so that the two routes have to agree. The second test is what
# stops the first from being vacuous: it shows that refitting the VAR on
# `newdata` would give a visibly different answer, so reproducing the
# by-hand value pins WHICH H was used and not merely that some H was.

test_that("prediction on newdata is the decision rule written out by hand", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac[1:120, ],
                       m = 2, var_lags = 2, restarts = 1, quiet = TRUE)
  nd <- dgp_pac[121:178, ]
  Tn <- nrow(nd)

  S <- .var_states(cbind(d.ystar = c(NA, diff(nd$ystar))), fit$var$p)
  Slag <- rbind(NA, S[-Tn, , drop = FALSE])
  hv <- .hvec(fit$alpha, fit$beta, fit$var$H, 2L)
  by_hand <- fit$a[1] * .lagv(nd$y - nd$ystar, 1) +
    fit$a[2] * .lagv(c(NA, diff(nd$y)), 1) +
    drop(Slag %*% hv)

  expect_equal(predict(fit, newdata = nd), by_hand)
  expect_true(sum(is.finite(predict(fit, newdata = nd))) > 40)

  # And in levels, which only adds y_{t-1} from the data.
  expect_equal(predict(fit, newdata = nd, type = "level"),
               .lagv(nd$y, 1) + by_hand)
})

test_that("predict does not re-estimate the auxiliary VAR on newdata", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac[1:120, ],
                       m = 2, var_lags = 2, restarts = 1, quiet = TRUE)
  nd <- dgp_pac[121:178, ]
  Tn <- nrow(nd)

  Xe <- cbind(d.ystar = c(NA, diff(nd$ystar)))
  Slag <- rbind(NA, .var_states(Xe, fit$var$p)[-Tn, , drop = FALSE])
  H_refit <- .var_companion(Xe, fit$var$p)$H

  # A VAR refitted on the 58 held-out observations is a different VAR ...
  expect_false(isTRUE(all.equal(H_refit, fit$var$H)))
  Z_refit <- drop(Slag %*% .hvec(fit$alpha, fit$beta, H_refit, 2L))
  Z_fitted <- drop(Slag %*% .hvec(fit$alpha, fit$beta, fit$var$H, 2L))
  # ... and the forward sum it implies is visibly different, so the previous
  # test distinguishes the two rather than passing either way.
  expect_false(isTRUE(all.equal(Z_refit, Z_fitted)))
})

test_that("predicting on the estimation data reproduces the in-sample call", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, restarts = 1, quiet = TRUE)
  expect_equal(predict(fit, newdata = dgp_pac), predict(fit))
})


# ---- the linear block is unpacked positionally ----
#
# `par` is c(theta, a_f?, delta), and predict() reads the tail of it by
# position, the same way recm_estimate() reads `r$lin`. Reproducing fitted()
# with BOTH optional pieces present is what pins that layout: with only one
# of them, several wrong orderings still agree.

test_that("fitted values survive free_forward plus an extra regressor", {
  fit <- recm_estimate("y", "ystar", vars = "w", beta = 0.995, data = dgp_pac,
                       m = 3, var_lags = 2, free_forward = TRUE,
                       restarts = 1, quiet = TRUE)
  expect_equal(predict(fit)[fit$index], fitted(fit))
  expect_equal(predict(fit, newdata = dgp_pac), predict(fit))
})

test_that("fitted values are reproduced with an expectations variable", {
  fit <- recm_estimate("y", "ystar", expectations = "income", beta = 0.995,
                       data = dgp_pac, m = 2, var_lags = 2, restarts = 1,
                       quiet = TRUE)
  expect_equal(fit$var$k, 2L)
  expect_equal(predict(fit)[fit$index], fitted(fit))
  expect_equal(predict(fit, newdata = dgp_pac), predict(fit))
})

test_that("the growth column is applied and is required on newdata", {
  dat <- dgp_pac
  dat$g <- .DGP_G_Y
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
                       var_lags = 2, growth = "g", restarts = 1, quiet = TRUE)
  expect_true(fit$growth)
  expect_equal(predict(fit)[fit$index], fitted(fit))
  # The correction is (1 - sum a_i - sum d_i) * g_t, so dropping the column
  # must fail loudly rather than predict the uncorrected path.
  expect_error(predict(fit, newdata = dat[, c("y", "ystar")]), "trend growth")
})


# ---- intervals ----

test_that("predict warns that conditional intervals are too narrow", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, restarts = 0, quiet = TRUE)
  expect_warning(predict(fit, interval = "confidence"),
                 "condition on the auxiliary VAR")
  expect_silent(predict(fit, interval = "confidence", quiet = TRUE))
  expect_error(predict(fit, interval = "confidence", level = 1), "level")
  expect_error(predict(fit, interval = "confidence", level = c(0.9, 0.95)),
               "level")
  expect_error(predict(fit, interval = "confidence", boot_object = list()),
               "recm_boot")
})

test_that("bootstrap intervals are reproducible across calls", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac[1:100, ],
                       m = 2, var_lags = 2, restarts = 0, quiet = TRUE)
  bs <- recm_boot(fit, R = 10, seed = 4)
  # The prediction interval is formed analytically, not by drawing residual
  # noise: a predict() that returned different numbers on each call would be
  # unusable for anything reported.
  expect_equal(predict(fit, interval = "prediction", boot_object = bs),
               predict(fit, interval = "prediction", boot_object = bs))
  ci <- predict(fit, interval = "confidence", boot_object = bs)
  pi <- predict(fit, interval = "prediction", boot_object = bs)
  ok <- fit$index
  expect_true(mean(pi$upr[ok] - pi$lwr[ok]) > mean(ci$upr[ok] - ci$lwr[ok]))
  wrong <- structure(list(replicates = cbind(1, 2, 3)), class = "recm_boot")
  expect_error(predict(fit, interval = "confidence", boot_object = wrong),
               "not from the same fit")
})
