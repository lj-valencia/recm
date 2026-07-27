# Covers INVARIANTS 3, 5, 6 and the conditioning warning above m = 5.
# See docs/04-testing.md.

test_that("k_0 equals A(1)A(beta) under the cc = 1 normalisation", {
  # INVARIANT 3, in the form the mathematics actually supports.
  #
  # Matching leading coefficients in z^m P(q(z)) = cc * z^m A(beta/z) A(z)
  # gives cc = k_m * prod(zeta_i); evaluating at z = 1, where q(1) = 0,
  # gives k_0 = cc * A(1)A(beta). The identity k_0 = A(1)A(beta) therefore
  # holds only when cc = 1, which is NOT the k_0 = 1 normalisation the
  # package estimates in. For the reference vector cc = 40.29781, so
  # k[1] = 1 while A(1)A(beta) = 0.02481524.
  #
  # `k_norm` is the input rescaled to cc = 1. See docs/02-math-spec.md
  # section 3 and the note in R/factorise.R.
  k <- .ref_k
  beta <- .ref_beta
  fs <- factorise_spectral(k, beta)
  s <- .scalars(fs$alpha, beta)
  expect_equal(fs$k_norm[1], s$A1 * s$Ab, tolerance = 1e-10)
  expect_equal(k[1] / fs$cc, s$A1 * s$Ab, tolerance = 1e-10)
})

test_that("h_i equal the long division of A(1)A(beta)/A(beta F)", {  # INV 5
  k <- .ref_k
  beta <- .ref_beta
  al <- factorise_spectral(k, beta)$alpha
  s <- .scalars(al, beta)
  w <- .dweights(al, beta, horizon = 50)

  # Recursive psi_i from phi = alpha * beta^(1:m): the power series
  # coefficients of 1 / A(beta F).
  m <- length(al)
  phi <- al * beta^seq_len(m)
  psi <- numeric(51)
  psi[1] <- 1
  for (i in seq_len(50)) {
    j <- seq_len(min(i, m))
    psi[i + 1] <- -sum(phi[j] * psi[i + 1 - j])
  }
  expect_equal(w$h, s$A1 * s$Ab * psi, tolerance = 1e-12)
})

test_that("Riccati and spectral agree for small m", {   # INVARIANT 6
  set.seed(606)
  for (m in 2:4) {
    k <- c(1, sort(runif(m, .5, 30), decreasing = TRUE))
    expect_equal(.lq_alpha(k, 0.98)$alpha,
                 factorise_spectral(k, 0.98)$alpha, tolerance = 1e-6)
  }
})

test_that("the spectral route is well conditioned at small m", {
  fs <- factorise_spectral(.ref_k, .ref_beta)
  # Reciprocity |zeta_in * zeta_out / beta - 1| stays near machine epsilon
  # here and degrades to about 1e-3 by m = 20.
  expect_lt(fs$reciprocity, 1e-10)
  expect_length(fs$roots, 2L)
  expect_true(all(Mod(fs$roots) > 1))
})

test_that("the spectral route warns above m = 5", {
  # Wilkinson ill-conditioning: the degree-2m polynomial's coefficients
  # span many orders of magnitude. The warning is mandatory, not advisory.
  k <- c(1, rep(5, 6))
  expect_warning(factorise_spectral(k, 0.98), "ill-conditioned")
})

test_that("the Riccati route stays well conditioned at large m", {
  k <- c(1, 10 * 0.6^(0:19))
  fit <- .lq_alpha(k, 0.98)
  expect_true(fit$converged)
  expect_length(fit$alpha, 20L)
  expect_true(all(is.finite(fit$alpha)))
  # Closed-loop eigenvalues are the reciprocals of the roots of A(z) and
  # must lie strictly inside the unit circle.
  expect_lt(fit$maxeig, 1)
})

test_that("inadmissible costs return NULL rather than erroring", {
  # Internal builders must let the optimiser reject a trial theta cheaply;
  # never stop() inside a function the optimiser calls.
  expect_null(.lq_alpha(c(1, -3, 2), 0.98))
  expect_null(.lq_alpha(c(1, 0, 0), 0.98))
  expect_null(factorise_spectral(c(1, -3, 2), 0.98))
})

test_that("m = 20 Riccati converges quickly", {
  skip_on_cran()
  k <- c(1, 10 * 0.6^(0:19))
  fit <- .lq_alpha(k, 0.98)
  # docs/04-testing.md records about 80 iterations at m = 20.
  expect_lt(fit$iterations, 500L)
})

test_that("a Riccati solve that cannot certify is rejected outright", {
  # Value iteration converges linearly at a rate approaching 1 for high m
  # with slowly decaying costs, so the per-step gain change understates the
  # remaining error by roughly 1/(1-rate). At m = 20, psi = 0.79 a per-step
  # change of 1.1e-07 sat 4.5e-05 from the converged gain. Such a result is
  # outside INVARIANT 6's 1e-6 and must NOT be returned as converged.
  k <- c(1, 4.58 * 0.79^(0:19))
  fit <- .lq_alpha(k, 0.995, maxit = 6000)
  expect_false(fit$converged)
})
