# Covers INVARIANT 4, the VAR companion layout, and INVARIANT 8.
# See docs/04-testing.md.

test_that("closed-form Z matches brute-force summation", {   # INVARIANT 4
  set.seed(404)
  beta <- .ref_beta
  alpha <- c(-0.93, 0.08)

  X <- cbind(dystar = rnorm(200, 0.0104, 0.004),
             income = rnorm(200, 0.0124, 0.008))
  vc <- .var_companion(X, p = 2)
  H <- vc$H
  z <- vc$states[nrow(vc$states), ]

  h <- .hvec(alpha, beta, H, sel = 2L)
  d <- .dweights(alpha, beta, 600)$d
  Hi <- H
  bf <- 0
  for (i in 0:600) {
    bf <- bf + d[i + 1] * (Hi %*% z)[2]
    Hi <- Hi %*% H
  }
  expect_equal(drop(h %*% z), bf, tolerance = 1e-10)
})

test_that("the state vector layout is 1, X_t, X_{t-1}, ...", {
  set.seed(405)
  X <- cbind(dystar = rnorm(120), income = rnorm(120))
  p <- 2L
  k <- 2L
  vc <- .var_companion(X, p)

  expect_equal(vc$n_z, 1L + k * p)
  expect_equal(ncol(vc$states), 1L + k * p)
  # The leading constant carries the VAR intercept: no demeaning anywhere.
  expect_equal(vc$H[1, 1], 1)
  expect_equal(vc$H[1, -1], rep(0, vc$n_z - 1L))
  # states is stored row-wise, one row per period, and the first p-1 rows
  # are NA because the lag block is not yet filled.
  expect_true(all(vc$states[p:nrow(X), 1] == 1))
  for (tt in c(p, 50L, nrow(X))) {
    expect_equal(unname(vc$states[tt, 2:(k + 1)]), unname(X[tt, ]))
    expect_equal(unname(vc$states[tt, (k + 2):(1 + 2 * k)]),
                 unname(X[tt - 1, ]))
  }
})

test_that("d(ystar) sits at position 2 of the state, as sel = 2L assumes", {
  set.seed(408)
  X <- cbind(dystar = rnorm(120), income = rnorm(120))
  vc <- .var_companion(X, p = 2)
  expect_equal(colnames(vc$states)[2], "dystar.l0")
  expect_equal(unname(vc$states[60, 2]), unname(X[60, 1]))
})

test_that("rho(H) excluding the constant block is reported separately", {
  set.seed(406)
  X <- cbind(dystar = rnorm(120), income = rnorm(120))
  vc <- .var_companion(X, p = 2)
  # The raw spectral radius is always 1 because of the constant block, so
  # rho_dyn is what should be displayed.
  expect_equal(vc$rho, 1, tolerance = 1e-8)
  expect_lt(vc$rho_dyn, 1)
})

test_that("the forward sum refuses to diverge", {       # INVARIANT 8
  set.seed(407)
  X <- cbind(dystar = rnorm(120), income = rnorm(120))
  H <- .var_companion(X, p = 2)$H
  # An alpha whose G has spectral radius at or above 1 must yield NULL,
  # never a silently truncated or divergent sum. G's eigenvalues here are
  # +/- sqrt(1.5), so rho(G) = 1.2247.
  bad <- c(0, -1.5)
  expect_gt(.scalars(bad, 1)$rhoG, 1)
  expect_null(.hvec(bad, 1, H, sel = 2L))
})
