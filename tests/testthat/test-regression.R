# Frozen numbers from a fixed fixture and seed. Catches unintended drift.
# See docs/04-testing.md, "Regression tests".
#
# When one of these fails, that is information, not an inconvenience:
#
#   1. work out WHICH invariant or numerical route changed
#   2. if the mathematics is now more correct, update the frozen value AND
#      add a NEWS.md entry explaining the change
#   3. if you cannot explain the change, the change is a bug
#
# Do not batch-update frozen values to make the suite green.

test_that("the factorisation of the reference cost vector is stable", {
  # Route: .lq_alpha (Riccati, gain-based convergence at .TOL_RICCATI),
  # m = 2, beta = 0.98, k = c(1, 27.6486, 3.2239).
  lq <- .lq_alpha(.ref_k, .ref_beta)
  expect_equal(lq$alpha, c(-0.9300004263, 0.08000186671),
               tolerance = 1e-8)

  # The spectral route must land in the same place. Its normalisation
  # constant is a separate frozen quantity because INVARIANT 3 is stated
  # in terms of it.
  fs <- factorise_spectral(.ref_k, .ref_beta)
  expect_equal(fs$cc, 40.29780968, tolerance = 1e-8)
})

test_that("the reference scalars are stable", {
  s <- .scalars(.lq_alpha(.ref_k, .ref_beta)$alpha, .ref_beta)
  expect_equal(s$A1, 0.1500014404, tolerance = 1e-9)
  expect_equal(s$Ab, 0.165433375, tolerance = 1e-9)
  expect_equal(s$sum_d, 0.837051537, tolerance = 1e-9)
  expect_equal(s$sum_h, 0.1500014404, tolerance = 1e-9)
  expect_equal(s$rhoG, 0.8174029702, tolerance = 1e-9)
  # sum_h = A(1)A(beta) iota'(I-G)^{-1} iota happens to equal A1 here.
  expect_equal(.halflife(.alpha_to_a(.lq_alpha(.ref_k, .ref_beta)$alpha))
               $halflife, 3.943326466, tolerance = 1e-7)
})

test_that("the reference NLS fit is stable", {
  # dgp_pac at its default seed, m = 2, beta = 0.995, var_lags = 2,
  # restarts = 3. Nelder-Mead, so the parameter tolerance is looser than
  # the objective tolerance by design -- two very different theta can give
  # nearly identical criteria.
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, restarts = 3, quiet = TRUE)
  expect_equal(fit$n, 175L)
  expect_equal(unname(fit$a), c(-0.1172746121, 0.1980108212),
               tolerance = 1e-5)
  expect_equal(fit$a_f, 0.7811280592, tolerance = 1e-5)
  expect_equal(fit$objective(fit$theta), 0.0008094822876, tolerance = 1e-8)
  expect_equal(sd(fit$residuals), 0.001965315746, tolerance = 1e-8)
})

test_that("the reference GMM fit is stable", {
  # These values MOVED in 0.0.0.9000, deliberately, for two reasons. The
  # intercept column of the instrument matrix used to be deleted by the
  # zero-variance filter, silently dropping the E[e_t] = 0 moment; it is now
  # retained. And iv_lag now defaults to var_lags + 2 rather than dating the
  # automatic instruments at t-1. Previously frozen at
  # a = c(-0.09958854565, 0.3249709416), F = 137.1755111, 4 instruments.
  # Not a tolerance change -- a specification change.
  #
  # The sample is 172, not the 175 of the NLS fit: the longer instrument
  # lags cost three observations.
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, method = "gmm", restarts = 1,
                       quiet = TRUE)
  expect_equal(fit$extras$iv_lag, 4L)
  expect_equal(fit$n, 172L)
  expect_equal(unname(fit$a), c(-0.0653096843, 0.3439159914),
               tolerance = 1e-4)
  expect_equal(fit$extras$n_instruments, 5L)
  expect_equal(fit$extras$first_stage_F, 54.74009946, tolerance = 1e-5)
})

test_that("the t-1 dated GMM fit is stable", {
  # The opt-out path, frozen so it cannot drift unnoticed. Instruments are
  # stronger here (F 137 against 55) and the sample is three longer -- that
  # is the efficiency the default gives up in exchange for validity when
  # ystar is measured with error.
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, method = "gmm", iv_lag = 1,
                       restarts = 1, quiet = TRUE)
  expect_equal(fit$n, 175L)
  expect_equal(unname(fit$a), c(-0.1067493569, 0.2007726722),
               tolerance = 1e-4)
  expect_equal(fit$extras$first_stage_F, 137.1755111, tolerance = 1e-5)
})

test_that("the fixture itself has not drifted", {
  # If this fails, every other frozen number in this file is meaningless.
  expect_equal(nrow(dgp_pac), 178L)
  expect_equal(mean(diff(dgp_pac$y)), 0.008835444247, tolerance = 1e-10)
  expect_equal(mean(diff(dgp_pac$ystar)), 0.008891236821, tolerance = 1e-10)
  expect_equal(attr(dgp_pac, "a0_true"), -0.1550013654, tolerance = 1e-9)
  expect_equal(attr(dgp_pac, "a1_true"), 0.07936951194, tolerance = 1e-9)
})
