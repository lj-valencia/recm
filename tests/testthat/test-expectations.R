test_that("the closed-form h matches a 500-term explicit sum", {
  phi <- make_companion(0.5, c(0.6, -0.2))
  ev <- matrix(c(0, 1, 0), 3L, 1L)
  set.seed(11)
  zlag <- cbind(1, matrix(stats::rnorm(80), 40L, 2L))

  for (a in list(0.35, c(0.4, 0.25), c(0.5, 0.28, 0.06))) {
    for (beta in c(0.95, 0.99, 1)) {
      alg <- recm:::pac_algebra(a, beta)
      skip_if_not(recm:::is_admissible(alg, max(Mod(eigen(phi)$values))))
      expect_equal(
        recm:::z_var(alg, phi, ev, zlag),
        recm:::z_var_truncated(alg, phi, ev, zlag, 500L),
        tolerance = 1e-10
      )
    }
  }
})

test_that("the forward term collapses to sum(d_i) * g on a constant path", {
  # A target growing at a constant rate g must load sum(d_i) * g, whichever
  # branch computes it. This is what ties "mce" to "var".
  g <- 0.5
  phi <- make_companion(g * (1 - 0.6), 0.6)
  ev <- matrix(c(0, 1), 2L, 1L)
  zlag <- matrix(c(1, g), 1L, 2L)

  for (a in list(0.35, c(0.4, 0.25), c(0.45, 0.2, 0.05))) {
    beta <- 0.97
    alg <- recm:::pac_algebra(a, beta)
    expect_equal(recm:::z_var(alg, phi, ev, zlag)[1], alg$d_sum * g,
                 tolerance = 1e-10)
    expect_equal(recm:::z_mce(alg, rep(g, 3000L), 1L), alg$d_sum * g,
                 tolerance = 1e-9)
  }
})

test_that("mce and var agree on a deterministic path once shared", {
  # With no future shocks both branches compute the same rational expectation.
  # The var branch conditions on t-1, so it is compared against the mce path
  # of the same date.
  beta <- 0.98
  a <- c(0.4, 0.25)
  alg <- recm:::pac_algebra(a, beta)
  ar <- 0.7
  const <- 0
  phi <- make_companion(const, ar)
  ev <- matrix(c(0, 1), 2L, 1L)

  n <- 60L
  dystar <- numeric(n)
  dystar[10L] <- 1
  for (t in 11:n) {
    dystar[t] <- ar * dystar[t - 1L]
  }
  # The pad must continue the autoregressive decay, not jump to zero: the mce
  # branch reads the whole path, so a discontinuity there is a real forecast
  # error rather than a truncation.
  path <- c(dystar, dystar[n] * ar^seq_len(3000L))

  z_m <- recm:::z_mce(alg, path, n)
  zlag <- cbind(1, c(0, dystar[-n]))
  z_v <- recm:::z_var(alg, phi, ev, zlag)

  # From the period after the shock the two information sets coincide. The
  # comparison is on absolute error: the mce branch truncates at the end of
  # the padded path, so late values are tiny and their relative error is not.
  expect_lt(max(abs(z_m[12:n] - z_v[12:n])), 1e-9)
})

test_that("the mce pad is long enough for the terminal condition", {
  alg <- recm:::pac_algebra(0.2, 0.99)
  g <- 0.4
  short <- recm:::z_mce(alg, rep(g, 20L), 1L)
  long <- recm:::z_mce(alg, rep(g, 4000L), 1L)
  expect_lt(abs(long - alg$d_sum * g), 1e-9)
  expect_gt(abs(short - alg$d_sum * g), 1e-9)
  expect_gte(recm:::mce_pad_length(alg), 200L)
})
