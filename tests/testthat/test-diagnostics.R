# recm_profile() and recm_boot(). Correctness of recm_boot() beyond "it
# runs" is roadmap item R-1 and is explicitly NOT claimed here; see
# docs/04-testing.md, "What is not tested yet".

fit_diag <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac,
                          m = 2, var_lags = 2, restarts = 1, quiet = TRUE)

test_that("recm_profile returns a grid bracketing the estimate", {
  pr <- recm_profile(fit_diag, which = 1, ngrid = 9)
  expect_s3_class(pr, "recm_profile")
  expect_length(pr$grid, 9L)
  expect_length(pr$value, 9L)
  expect_gte(pr$estimate, min(pr$grid))
  expect_lte(pr$estimate, max(pr$grid))
  # The profiled criterion is minimised at or near the estimate.
  expect_lte(min(pr$value), fit_diag$objective(fit_diag$theta) + 1e-6)
})

test_that("recm_profile reports flatness rather than hiding it", {
  pr <- recm_profile(fit_diag, which = 1, ngrid = 9)
  expect_true(is.finite(pr$rel_range))
  out <- paste(capture.output(print(pr)), collapse = " ")
  expect_match(out, "relative range")
  if (pr$rel_range < 0.02) expect_match(out, "weakly identified")
})

test_that("recm_profile handles a one-parameter free block via Brent", {
  pr <- recm_profile(fit_diag, which = 2, ngrid = 5)
  expect_true(all(is.finite(pr$value)))
})

test_that("recm_boot runs and returns a covariance of the right shape", {
  skip_on_cran()
  # R is deliberately tiny: this asserts that the machinery runs and is
  # shaped correctly, NOT that the bootstrap is calibrated. Calibration is
  # roadmap item R-1.
  bs <- recm_boot(fit_diag, R = 6, seed = 42)
  expect_s3_class(bs, "recm_boot")
  expect_equal(dim(bs$cov), c(fit_diag$npar, fit_diag$npar))
  expect_length(bs$se, fit_diag$npar)
  expect_gte(bs$R_ok, 2L)
  expect_equal(bs$R_ok + bs$R_failed, 6L)
  expect_true(all(is.finite(bs$se)))
})

test_that("recm_boot refuses a growth-corrected fit rather than guessing", {
  d <- dgp_pac
  d$g <- 0.0104
  fg <- recm_estimate("y", "ystar", beta = 0.995, data = d, m = 2,
                      var_lags = 2, growth = "g", restarts = 0,
                      quiet = TRUE)
  expect_error(recm_boot(fg, R = 2), "growth")
})

test_that("plot methods draw without error", {
  skip_if_not(capabilities("png"), "no png device")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  pr <- recm_profile(fit_diag, which = 1, ngrid = 5)
  expect_silent(invisible(plot(pr)))
})
