# Covers INVARIANTS 1, 2, 7 and the .G_mat / .scalars / .halflife helpers.
# See docs/04-testing.md.

test_that("alpha <-> a round-trips exactly", {          # INVARIANT 1
  set.seed(101)
  for (m in 1:6) {
    a <- c(-runif(1, .05, .5), runif(m - 1, -.4, .4))
    expect_equal(a, .alpha_to_a(.a_to_alpha(a)), tolerance = 1e-12)
  }
})

test_that("the round trip also holds starting from alpha", {  # INVARIANT 1
  set.seed(102)
  for (m in 1:6) {
    alpha <- runif(m, -.6, .3)
    expect_equal(alpha, .a_to_alpha(.alpha_to_a(alpha)), tolerance = 1e-12)
  }
})

test_that("a_0 = -A(1), with the minus sign", {         # INVARIANT 2
  alpha <- c(-0.6, 0.15)
  a <- .alpha_to_a(alpha)
  expect_equal(a[1], -(1 + sum(alpha)), tolerance = 1e-12)
  # A stable equation has a_0 < 0, hence A(1) > 0. The FRB/US "PAC Basics"
  # note prints a_0 = A(1); that is a typo in the source.
  expect_lt(a[1], 0)
})

test_that("G is the companion matrix of A(beta F)", {
  alpha <- c(-0.6, 0.15)
  beta <- .ref_beta
  G <- .G_mat(alpha, beta)
  expect_equal(dim(G), c(2L, 2L))
  expect_equal(G[2, ], c(-alpha[2] * beta^2, -alpha[1] * beta),
               tolerance = 1e-12)
  expect_equal(G[1, ], c(0, 1), tolerance = 1e-12)
})

test_that(".scalars agrees with direct evaluation of A(1) and A(beta)", {
  alpha <- c(-0.6, 0.15)
  beta <- .ref_beta
  s <- .scalars(alpha, beta)
  expect_equal(s$A1, 1 + sum(alpha), tolerance = 1e-12)
  expect_equal(s$Ab, 1 + sum(alpha * beta^seq_along(alpha)),
               tolerance = 1e-12)
})

test_that("the closed-form totals match summed weights", {
  alpha <- c(-0.93, 0.08)
  beta <- .ref_beta
  s <- .scalars(alpha, beta)
  w <- .dweights(alpha, beta, horizon = 3000)
  expect_equal(sum(w$h), s$sum_h, tolerance = 1e-10)
  expect_equal(sum(w$d), s$sum_d, tolerance = 1e-10)
})

test_that("normalised lead weights sum to one", {       # INVARIANT 7
  alpha <- c(-0.93, 0.08)
  w <- .dweights(alpha, .ref_beta, horizon = 2000)
  f <- w$d / sum(w$d)
  expect_equal(sum(f), 1, tolerance = 1e-12)
})

test_that(".halflife is continuous, so its numerical derivative is nonzero", {
  # An integer step would give a delta-method SE of exactly 0, which is
  # what happened before the interpolation was added.
  a <- c(-0.25, 0.30)
  h1 <- .halflife(a)$halflife
  h2 <- .halflife(a + c(1e-6, 0))$halflife
  expect_true(is.finite(h1))
  expect_gt(abs(h2 - h1), 0)
})

test_that(".halflife returns the gap path starting at one", {
  hl <- .halflife(c(-0.25, 0.30), horizon = 40)
  expect_equal(hl$path[1], 1)
  expect_length(hl$path, 41L)
  # A stable rule closes the gap.
  expect_lt(abs(hl$path[41]), 0.5)
})
