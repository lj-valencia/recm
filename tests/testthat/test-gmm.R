# Covers the endogeneity correction, instrument dating, and Hansen J.
# See docs/04-testing.md.
#
# Note on instrument dating: Z_t = h'z_{t-1} and z_{t-1} holds d(ystar) at
# lags 0..p-1, so if ystar is measured with error nu the composite error
# carries nu_{t-1} through nu_{t-p-1}. Instruments built from y or ystar
# must therefore be dated t-(p+2) or earlier. `iv_lag` dates the AUTOMATIC
# set; instruments the user supplies are still the user's responsibility,
# because recognising the provenance of an arbitrary instrument expression
# is not something the package can do reliably.
#
# iv_lag defaults to 1, which is NOT the dating rule. That is deliberate and
# is tested below: the rule-respecting set has valid moments but does not
# identify theta. See the Instrument dating section of ?recm_estimate and
# roadmap R-10.

test_that("GMM runs and reports its diagnostics", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, method = "gmm", restarts = 1,
                       quiet = TRUE)
  expect_equal(fit$method, "gmm")
  expect_true(is.finite(fit$extras$J))
  expect_true(is.finite(fit$extras$hansen_p))
  expect_true(is.finite(fit$extras$first_stage_F))
  expect_gt(fit$extras$n_instruments, 0L)
  expect_true(all(is.finite(fit$residuals)))
})

test_that("the GMM criterion is what the profile closure evaluates", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, method = "gmm", restarts = 1,
                       quiet = TRUE)
  expect_true(is.finite(fit$objective(fit$theta)))
  # The optimum should not be beaten by a nearby point by much.
  expect_lte(fit$objective(fit$theta),
             fit$objective(fit$theta + 0.5) + 1e-8)
})

test_that("GMM is consistent on a known PAC DGP", {
  skip_on_cran()
  # NOT tested here: that GMM removes the bias measurement error in ystar
  # induces in NLS. It does not, at either dating. See the two tests below.
  big <- simulate_dgp_pac(n = 4000)
  fit <- recm_estimate("y", "ystar", beta = attr(big, "beta"), data = big,
                       m = 2, var_lags = 2, method = "gmm", restarts = 3,
                       quiet = TRUE)
  expect_lt(abs(fit$a[1] - attr(big, "a0_true")), 0.03)
})

test_that("iv_lag dates only the y/ystar-derived instrument columns", {
  # The dating rule applies to instruments built from y or ystar. Other
  # state variables are exogenous with respect to nu and must stay at t-1,
  # where they are far stronger. Column names carry the TOTAL lag.
  nm <- function(L) {
    fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                         expectations = ~ income, var_lags = 2,
                         method = "gmm", iv_lag = L, restarts = 0,
                         quiet = TRUE)
    colnames(get("Zi", envir = environment(fit$objective)))
  }
  expect_equal(nm(1), c("const", "ecm.l1", "dy.l1", "d.ystar.l1",
                        "income.l1", "d.ystar.l2", "income.l2"))
  # ystar-derived columns move to t-4 and t-5; income stays at t-1, t-2.
  expect_equal(nm(4), c("const", "ecm.l4", "dy.l4", "d.ystar.l4",
                        "income.l1", "d.ystar.l5", "income.l2"))
})

test_that("the intercept survives the zero-variance instrument filter", {
  # It is the one column that is meant to be constant. Dropping it deletes
  # the E[e_t] = 0 moment, which is what used to happen.
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, method = "gmm", restarts = 0,
                       quiet = TRUE)
  Zi <- get("Zi", envir = environment(fit$objective))
  expect_equal(colnames(Zi)[1], "const")
  expect_true(all(Zi[, 1] == 1))
})

test_that("the default instrument set has valid moments", {
  skip_on_cran()
  # At iv_lag = var_lags + 2 the residual must be orthogonal to every
  # instrument. This is the property the whole dating rule exists to buy.
  big <- simulate_dgp_pac(n = 3000, seed = 8101)
  cm <- list(y = "y", ystar = "ystar", beta = attr(big, "beta"), data = big,
             m = 2, expectations = ~ income, var_lags = 2, restarts = 1,
             quiet = TRUE)
  f_nls <- do.call(recm_estimate, c(cm, list(method = "nls")))
  f_4 <- do.call(recm_estimate, c(cm, list(method = "gmm")))
  expect_equal(f_4$extras$iv_lag, 4L)

  e <- residuals(f_nls)
  Zi <- get("Zi", envir = environment(f_4$objective))
  n <- min(nrow(Zi), length(e))
  ps <- apply(tail(Zi, n)[, -1, drop = FALSE], 2,
              function(z) stats::cor.test(z, tail(e, n))$p.value)
  expect_true(all(ps > 0.05))
})

test_that("correct dating corrects measurement-error bias; t-1 does not", {
  skip_on_cran()
  # The reason iv_lag defaults to the rule. With ystar badly mismeasured,
  # NLS and the t-1 instrument set are both inconsistent -- at T = 3000
  # their bias is about +0.087 and does not shrink -- while the correctly
  # dated set lands within 0.02. Bounds are loose because the point is the
  # order of magnitude, not the third digit.
  big <- simulate_dgp_pac(n = 3000, seed = 9601, measurement_error = TRUE,
                          sd_nu = 0.02)
  a0 <- attr(big, "a0_true")
  cm <- list(y = "y", ystar = "ystar", beta = attr(big, "beta"), data = big,
             m = 2, expectations = ~ income, var_lags = 2, restarts = 1,
             quiet = TRUE)
  f_nls <- do.call(recm_estimate, c(cm, list(method = "nls")))
  f_1 <- do.call(recm_estimate, c(cm, list(method = "gmm", iv_lag = 1)))
  f_4 <- do.call(recm_estimate, c(cm, list(method = "gmm")))

  expect_gt(abs(f_nls$a[1] - a0), 0.05)
  expect_gt(abs(f_1$a[1] - a0), 0.05)   # t-1 buys nothing over NLS
  expect_lt(abs(f_4$a[1] - a0), 0.03)
})

test_that("dating below the rule warns", {
  expect_warning(
    recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                  var_lags = 2, method = "gmm", iv_lag = 1, restarts = 0),
    "inside the window"
  )
})

test_that("first-stage F falls as the expectations VAR lengthens", {
  skip_on_cran()
  # A longer expectations VAR buys forecast accuracy and pays in instrument
  # strength. The reference dataset went 13.1 -> 6.0 -> 3.3 at p = 1, 2, 4;
  # this fixture has a far more predictable ecm so the levels are much
  # higher, but the direction is the documented one.
  ff <- vapply(c(1, 2, 4), function(p) {
    recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                  var_lags = p, method = "gmm", restarts = 0,
                  quiet = TRUE)$extras$first_stage_F
  }, numeric(1))
  expect_true(all(diff(ff) < 0))
})

test_that("a strong first stage does not warn", {
  # The warning fires below F = 10. On this fixture F is above 100, so
  # silence is the correct behaviour -- this pins the branch so a future
  # change cannot start warning unconditionally.
  expect_no_warning(
    recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                  var_lags = 2, method = "gmm", restarts = 0)
  )
})

test_that("collinear instrument columns are dropped by QR pivot", {
  fit <- recm_estimate("y", "ystar", vars = ~ w, beta = 0.995,
                       data = dgp_pac, m = 2, var_lags = 2, method = "gmm",
                       instruments = ~ w + I(2 * w), restarts = 0,
                       quiet = TRUE)
  expect_lt(fit$extras$n_instruments, fit$extras$n_instruments_supplied)
})

test_that("Hansen J has the right degrees of freedom", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, method = "gmm", restarts = 0,
                       quiet = TRUE)
  expect_equal(fit$extras$Jdf, fit$extras$n_instruments - fit$npar)
  expect_gte(fit$extras$J, 0)
})
