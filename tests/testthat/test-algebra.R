test_that("a and alpha are exact inverses", {
  for (a in list(0.35, c(0.4, 0.25), c(0.5, 0.3, -0.1),
                 c(0.45, 0.2, 0.1, 0.05))) {
    expect_equal(recm:::alpha_to_a(recm:::a_to_alpha(a)), a)
  }
})

test_that("a_to_alpha reproduces Dynare's a2alpha", {
  # The last alpha equals the last a, interior ones are first differences of
  # a running backwards, and the first subtracts an extra one.
  a <- c(0.5, 0.3, -0.1)
  expect_equal(recm:::a_to_alpha(a), c(0.5 - 0.3 - 1, 0.3 - (-0.1), -0.1))
  # The m = 1 case Dynare cannot express: alpha_1 = a0 - 1.
  expect_equal(recm:::a_to_alpha(0.35), 0.35 - 1)
})

test_that("G matches the hand-written companion of A(beta F)", {
  alpha <- c(-0.6, 0.15)
  beta <- 0.98
  expect_equal(
    recm:::g_matrix(alpha, beta),
    rbind(c(0, 1), -rev(alpha * beta^(1:2)))
  )
  expect_equal(recm:::g_matrix(-0.7, 0.98), matrix(0.7 * 0.98, 1, 1))
})

test_that("m = 1 reproduces the closed forms of the notes", {
  beta <- 0.98
  theta <- 0.353245
  lambda <- 1 - theta
  alg <- recm:::pac_algebra(theta, beta)
  d <- recm:::lead_weights(alg, 60L)

  expect_equal(d[1], theta, tolerance = 1e-10)
  expect_equal(d[2] / d[1], beta * lambda, tolerance = 1e-10)
  expect_equal(alg$d_sum, theta / (1 - beta * lambda), tolerance = 1e-10)
  expect_equal(alg$a1, theta, tolerance = 1e-12)
  expect_equal(alg$cc, alg$a1 * alg$ab, tolerance = 1e-12)
})

test_that("the worked check of the notes reproduces to six digits", {
  # beta = 0.98, b = 5 gives lambda1 = 0.646755, theta = 0.353245,
  # d0 = 0.353245, decay = 0.633820, sum(d_i) = 0.964676.
  beta <- 0.98
  b <- 5
  lambda <- ((1 + b + beta * b) -
               sqrt((1 + b + beta * b)^2 - 4 * beta * b^2)) / (2 * beta * b)
  expect_equal(lambda, 0.646755, tolerance = 1e-6)
  alg <- recm:::pac_algebra(1 - lambda, beta)
  d <- recm:::lead_weights(alg, 10L)
  expect_equal(d[1], 0.353245, tolerance = 1e-6)
  expect_equal(d[2] / d[1], 0.633820, tolerance = 1e-6)
  expect_equal(alg$d_sum, 0.964676, tolerance = 1e-6)
})

test_that("lead weights sum and mean lead match their closed forms", {
  alg <- recm:::pac_algebra(c(0.4, 0.25), 0.97)
  d <- recm:::lead_weights(alg, 1500L)
  expect_equal(sum(d), alg$d_sum, tolerance = 1e-10)
  expect_equal(sum(seq_along(d) * d - d), alg$lead_sum, tolerance = 1e-10)
  expect_equal(sum(d / alg$d_sum), 1, tolerance = 1e-10)
})

test_that("structural costs invert the lag polynomial", {
  for (a in list(0.35, c(0.4, 0.25), c(0.5, 0.28, 0.06))) {
    for (beta in c(0.95, 0.98, 1)) {
      alg <- recm:::pac_algebra(a, beta)
      cost <- recm:::cost_params(alg$alpha, beta)
      expect_lt(cost$residual, 1e-10)
      # k0 = A(1) A(beta) is the numerical check of the notes.
      expect_equal(cost$k[1], alg$cc, tolerance = 1e-10)
    }
  }
})

test_that("cost positivity is reported, not assumed", {
  # A freely estimated reduced form need not come from a convex adjustment
  # cost problem, so `positive` is a diagnostic rather than a guarantee.
  expect_true(recm:::cost_params(recm:::a_to_alpha(0.35), 0.98)$positive)
  expect_true(
    recm:::cost_params(recm:::a_to_alpha(c(0.4, 0.25)), 0.98)$positive
  )
  expect_false(
    recm:::cost_params(recm:::a_to_alpha(c(0.5, 0.28, 0.06)), 0.98)$positive
  )
})

test_that("at m = 1 the recovered cost matches b = l / ((1-l)(1-beta l))", {
  beta <- 0.98
  lambda <- 0.646755
  alg <- recm:::pac_algebra(1 - lambda, beta)
  cost <- recm:::cost_params(alg$alpha, beta)
  expect_equal(
    cost$b[1], lambda / ((1 - lambda) * (1 - beta * lambda)),
    tolerance = 1e-8
  )
  expect_equal(cost$b[1], 5, tolerance = 1e-5)
})

test_that("admissibility rejects explosive and non-summable parameters", {
  # a0 = 0.35 with a well behaved auxiliary model is admissible.
  expect_true(recm:::is_admissible(recm:::pac_algebra(0.35, 0.98), 1))
  # a0 = 2.5 puts the root of A(z) inside the unit circle.
  expect_false(recm:::is_admissible(recm:::pac_algebra(2.5, 0.98), 1))
  # rho(G) rho(Phi) >= 1 breaks the closed-form forward sum.
  expect_false(recm:::is_admissible(recm:::pac_algebra(0.05, 0.99), 1.2))
})

test_that("half-life interpolates rather than returning the integer crossing", {
  # At m = 1 the gap decays at lambda per period.
  lambda <- 0.7
  expect_equal(recm:::half_life(1 - lambda), log(0.5) / log(lambda),
               tolerance = 1e-8)
  expect_false(recm:::half_life(0.35) == round(recm:::half_life(0.35)))
})
