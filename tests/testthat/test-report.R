# Covers both reporting formats, INVARIANT 7, and the equation()
# round-trip. See docs/04-testing.md.

fit_report <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac,
                            m = 3, var_lags = 2, restarts = 1,
                            quiet = TRUE)

test_that("summary prints both reporting formats", {
  out <- capture.output(print(summary(fit_report)))
  expect_true(any(grepl("FORMAT 1: EXPLICIT", out)))
  expect_true(any(grepl("FORMAT 2: COMPRESSED", out)))
  expect_true(any(grepl("Fit and diagnostics", out)))
})

test_that("the generated-regressor caveat is always printed", {
  # Do NOT add a switch to suppress this. It is always true, and it is
  # listed under "explicitly out of scope" in docs/05-roadmap.md.
  out <- paste(capture.output(print(summary(fit_report))), collapse = " ")
  expect_match(out, "generated regressor")
  expect_match(out, "recm_boot")
})

test_that("normalised lead weights sum to one", {       # INVARIANT 7
  # f_i = d_i / sum_d uses the CLOSED-FORM total, so a truncated vector
  # reaches 1 only up to a tail decaying like rho(G)^i. At horizon 2000
  # that tail is below machine epsilon.
  lw <- lead_weights(fit_report, horizon = 2000, normalised = TRUE)
  expect_equal(sum(lw$weight), 1, tolerance = 1e-12)
})

test_that("equation() round-trips the stored coefficients", {
  eq <- capture.output(res <- equation(fit_report, format = "both"))
  expect_equal(res$a, fit_report$a, tolerance = 1e-12)
  expect_equal(res$f, res$d / fit_report$scalars$sum_d, tolerance = 1e-12)
  # a_f = sum(d_i) under the Euler restriction.
  expect_equal(res$a_f, fit_report$scalars$sum_d, tolerance = 1e-12)
})

test_that("the two formats describe the same equation", {
  # Explicit fixes the forward coefficient at 1 and carries d_i; compressed
  # carries a_f and f_i. a_f * f_i must equal d_i.
  invisible(capture.output(res <- equation(fit_report)))
  expect_equal(res$a_f * res$f, res$d, tolerance = 1e-12)
})

test_that("max_lags truncates and reports the suppressed count", {
  big <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 8,
                       var_lags = 2, restarts = 0, quiet = TRUE)
  out <- paste(capture.output(equation(big, max_lags = 3)), collapse = " ")
  expect_match(out, "further lags suppressed")
})

test_that("lead weight standard errors are not identically zero", {
  # .halflife() interpolates the 0.5 crossing because an integer step has
  # zero numerical derivative and the delta-method SE came back as exactly
  # 0 before the fix. Guard the same failure mode here.
  lw <- lead_weights(fit_report, horizon = 12)
  expect_true(all(is.finite(lw$se)))
  expect_true(all(lw$se[-1] > 0))
  expect_equal(nrow(lw), 13L)
  expect_true(all(lw$lo <= lw$weight & lw$weight <= lw$hi))
})

test_that("lead_weights exposes the Jacobian for R-3", {
  lw <- lead_weights(fit_report, horizon = 6)
  J <- attr(lw, "jacobian")
  expect_equal(dim(J), c(7L, length(fit_report$theta)))
})

test_that("print and the standard extractors work", {
  expect_output(print(fit_report), "REC equation")
  expect_length(coef(fit_report), fit_report$npar)
  expect_equal(dim(vcov(fit_report)),
               c(fit_report$npar, fit_report$npar))
  expect_length(residuals(fit_report), fit_report$n)
  expect_equal(fitted(fit_report) + residuals(fit_report),
               fit_report$dep, tolerance = 1e-12)
  expect_equal(dim(confint(fit_report)), c(fit_report$npar, 2L))
  expect_true(is.finite(as.numeric(logLik(fit_report))))
})

test_that("plot draws without error", {
  skip_if_not(capabilities("png"), "no png device")
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  expect_silent(invisible(plot(fit_report, horizon = 8)))
})
