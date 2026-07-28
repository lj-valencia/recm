test_that("growth neutrality is automatic at a unit discount factor", {
  # At beta = 1 the operator A(beta F) A(L) is symmetric in L, so a path
  # linear in t passes through exactly and R(a) vanishes for every a and m.
  for (a in list(0.35, c(0.4, 0.25), c(0.5, 0.28, 0.06),
                 c(0.45, 0.2, 0.1, 0.03))) {
    expect_equal(recm:::gn_gap(a, 1), 0, tolerance = 1e-10)
    expect_lt(max(abs(recm:::gn_grad(a, 1))), 1e-6)
  }
})

test_that("growth neutrality binds when the discount factor is below one", {
  for (a in list(c(0.4, 0.25), c(0.5, 0.28, 0.06))) {
    gap <- recm:::gn_gap(a, 0.95)
    expect_gt(abs(gap), 1e-3)
    expect_gt(max(abs(recm:::gn_grad(a, 0.95))), 1e-3)
  }
})

test_that("at m = 1 and beta < 1 only a0 = 1 satisfies the restriction", {
  root <- stats::uniroot(
    function(a0) recm:::gn_gap(a0, 0.95),
    c(0.05, 1.5), tol = 1e-12
  )$root
  expect_equal(root, 1, tolerance = 1e-6)
})

test_that("the gap equals the permanent level gap of a simulated path", {
  # Simulate the equation deterministically on a balanced growth path and
  # confirm the steady-state deviation is -R(a) g / a0.
  a <- c(0.4, 0.25)
  beta <- 0.95
  g <- 0.5
  alg <- recm:::pac_algebra(a, beta)
  n <- 400L
  ystar <- 100 + g * seq_len(n)
  y <- numeric(n)
  y[1:2] <- ystar[1:2] - 1.5
  dy <- numeric(n)
  for (t in 3:n) {
    dy[t] <- a[1L] * (ystar[t - 1L] - y[t - 1L]) + a[2L] * dy[t - 1L] +
      alg$d_sum * g
    y[t] <- y[t - 1L] + dy[t]
  }
  expect_equal(y[n] - ystar[n], recm:::gn_implied_gap(a, beta, g),
               tolerance = 1e-6)
})

test_that("the restriction locus is non-degenerate at m = 2", {
  # For each a1 there is an a0 solving R(a) = 0, so the restriction traces a
  # curve in coefficient space rather than pinning down a point.
  beta <- 0.95
  a0 <- vapply(
    seq(0.05, 0.35, by = 0.05),
    function(a1) {
      stats::uniroot(
        function(v) recm:::gn_gap(c(v, a1), beta),
        c(1e-4, 0.999), tol = 1e-12
      )$root
    },
    numeric(1)
  )
  expect_true(all(is.finite(a0)))
  expect_true(all(a0 > 0 & a0 < 1))
  expect_gt(stats::sd(a0), 1e-4)
})
