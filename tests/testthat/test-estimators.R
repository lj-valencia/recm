# The estimator extension point (docs/01-architecture.md).
#
# A branch in recm_estimate() must set seven things: th, r, par_all, V, fit,
# extras and objective. Only the first, fourth and sixth were documented,
# and the object assembly reads all seven. These tests assert the contract
# in its positive form -- every implemented estimator produces a complete
# recmfit -- so a new branch that sets only what docs/01 used to name fails
# here rather than producing a half-built object.

test_that("every implemented estimator satisfies the branch contract", {
  for (mth in c("nls", "gmm")) {
    fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                         var_lags = 2, method = mth, restarts = 0,
                         quiet = TRUE)
    info <- mth

    # th and par_all
    expect_true(is.numeric(fit$theta) && length(fit$theta) == 2L,
                info = info)
    expect_equal(length(fit$par), fit$npar, info = info)
    expect_false(any(is.na(names(fit$par))), info = info)

    # V, conformable and named
    expect_equal(dim(fit$vcov), c(fit$npar, fit$npar), info = info)
    expect_equal(rownames(fit$vcov), names(fit$par), info = info)

    # r -- the residual object feeds every one of these
    expect_length(fit$residuals, fit$n)
    expect_length(fit$fitted, fit$n)
    expect_true(all(is.finite(fit$alpha)), info = info)
    expect_true(all(fit$k > 0), info = info)

    # fit -- only $convergence is read, and it must have been
    expect_true(is.logical(fit$converged) && !is.na(fit$converged),
                info = info)

    # extras -- a list, possibly empty
    expect_true(is.list(fit$extras), info = info)

    # objective -- a closure, finite at the optimum
    expect_true(is.function(fit$objective), info = info)
    expect_true(is.finite(fit$objective(fit$theta)), info = info)

    # and the fit must actually be a fit
    expect_lt(fit$a[1], 0)                     # INVARIANT 2
  }
})

test_that("only .METHODS_IV estimators build instruments", {
  # The instrument block used to gate on a bare `method == "gmm"`, an
  # invisible second edit site; a new estimator reached it with Ziv = NULL
  # and died inside apply() with "dim(X) must have a positive length".
  expect_true("gmm" %in% .METHODS_IV)
  expect_false("nls" %in% .METHODS_IV)

  f_nls <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                         var_lags = 2, method = "nls", restarts = 0,
                         quiet = TRUE)
  f_gmm <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                         var_lags = 2, method = "gmm", restarts = 0,
                         quiet = TRUE)
  # nls carries no instrument diagnostics; gmm carries the full set.
  expect_null(f_nls$extras$iv_lag)
  expect_null(f_nls$extras$n_instruments)
  expect_false(is.null(f_gmm$extras$iv_lag))
  expect_false(is.null(f_gmm$extras$n_instruments))
})

test_that("an unknown estimator is rejected", {
  expect_error(
    recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                  var_lags = 2, method = "cue", restarts = 0, quiet = TRUE),
    "should be one of"
  )
})

test_that("the summary reports instrument diagnostics only when they exist", {
  # print.summary.recmfit tests each extras field for NULL rather than
  # branching on method, which is what lets a new estimator supply its own
  # diagnostics without touching methods.R. Pin that.
  s_nls <- utils::capture.output(
    print(summary(recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac,
                                m = 2, var_lags = 2, restarts = 0,
                                quiet = TRUE)))
  )
  s_gmm <- utils::capture.output(
    print(summary(recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac,
                                m = 2, var_lags = 2, method = "gmm",
                                restarts = 0, quiet = TRUE)))
  )
  expect_false(any(grepl("Hansen J", s_nls)))
  expect_false(any(grepl("Instruments:", s_nls)))
  expect_true(any(grepl("Hansen J", s_gmm)))
  expect_true(any(grepl("Instruments:", s_gmm)))
  # The generated-regressor caveat is unconditional -- both, always.
  expect_true(any(grepl("generated regressor", s_nls, ignore.case = TRUE)))
  expect_true(any(grepl("generated regressor", s_gmm, ignore.case = TRUE)))
})
