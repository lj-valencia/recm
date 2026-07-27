# Covers the m = 2 boundary formulas and rejection of inadmissible theta.
# See docs/02-math-spec.md section 9 and docs/04-testing.md.

test_that("at m = 2, k_2 equals a_1 exactly under the cc = 1 scaling", {
  # docs/02 section 9: k_2 / cc = a_1. As with INVARIANT 3, the identity
  # holds for the cost vector scaled so the factorisation constant is 1,
  # not for the k_0 = 1 vector the package estimates in.
  beta <- .ref_beta
  k <- c(1, 12, 2.5)
  fs <- factorise_spectral(k, beta)
  a <- .alpha_to_a(fs$alpha)
  expect_equal(fs$k_norm[3], a[2], tolerance = 1e-10)
  expect_equal(k[3] / fs$cc, a[2], tolerance = 1e-10)
})

test_that("the m = 2 upper bound on a_1 matches the closed form", {
  beta <- .ref_beta
  a0 <- -0.25
  B <- 1 + beta - beta * a0
  ub <- (B - sqrt(B^2 - 4 * beta * (1 + a0))) / (2 * beta)

  # k_1 > 0 <=> (1 + a_0 + a_1)(1 + beta a_1) > 2 (1 + beta) a_1
  k1_of <- function(a1) {
    (1 + a0 + a1) * (1 + beta * a1) - 2 * (1 + beta) * a1
  }
  expect_gt(k1_of(ub * 0.9), 0)
  expect_lt(k1_of(ub * 1.1), 0)
})

test_that("a negative a_1 is not rationalisable at m = 2", {
  # k_2 = a_1 (up to the cc scaling) at m = 2, so a_1 > 0 is required: a
  # negative coefficient on d y_{t-1} cannot be produced by any admissible
  # cost vector. Round-tripping through the factorisation must never yield
  # one.
  set.seed(902)
  for (i in 1:20) {
    k <- c(1, runif(1, .5, 30), runif(1, .5, 30))
    a <- .alpha_to_a(.lq_alpha(k, .ref_beta)$alpha)
    expect_gt(a[2], 0)
  }
})

test_that("inadmissible cost vectors return NULL, not an error", {
  # Internal builders must let the optimiser reject a trial theta cheaply;
  # never stop() inside a function the optimiser calls.
  expect_null(.lq_alpha(c(1, -3, 2), .ref_beta))
  expect_null(.lq_alpha(c(1, 0, 0), .ref_beta))
  expect_null(.scalars(c(NaN, NaN), .ref_beta))
})

test_that("the fitted cost vector is admissible", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       var_lags = 2, restarts = 3, quiet = TRUE)
  # TODO(R-2): map these constraints into the unconstrained
  # parameterisation and supply lower/upper to nlminb, so the optimiser
  # never has to be rescued by the .PENALTY_INADMISSIBLE cliff.
  expect_true(all(fit$k > 0))
  expect_lt(fit$scalars$rhoG, 1)
  expect_gt(fit$a[2], 0)
})

test_that("cost admissibility is reported for general m", {
  # There is no closed form beyond m = 2; check all(k > 0) and report it.
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 4,
                       var_lags = 2, restarts = 1, quiet = TRUE)
  expect_true(all(fit$k > 0))
  out <- paste(capture.output(print(summary(fit))), collapse = " ")
  expect_match(out, "Cost admissibility")
})
