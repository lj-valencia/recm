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

test_that(".var_states builds the same states .var_companion does", {
  # The layout loop was split out of .var_companion() so a caller holding an
  # already-estimated H can build states on new observations. Splitting it is
  # only safe while the two agree exactly at every lag order.
  set.seed(406)
  X <- cbind(dystar = rnorm(120), income = rnorm(120))
  for (p in 1:4) {
    expect_equal(.var_states(X, p), .var_companion(X, p)$states)
  }
  # Fewer observations than lags is a specification error, not a boundary
  # for the optimiser to reject: p:Tn would run backwards and fill garbage.
  expect_error(.var_states(X[1:3, , drop = FALSE], 4L), "lag order")
})

test_that("a mechanism can run the fitted H over different states", {
  set.seed(407)
  X <- cbind(dystar = rnorm(200))
  vc <- .var_companion(X[1:120, , drop = FALSE], p = 2)
  Xnew <- X[121:200, , drop = FALSE]
  Snew <- .var_states(Xnew, vc$p)
  alpha <- .lq_alpha(.ref_k, .ref_beta)$alpha
  s <- .scalars(alpha, .ref_beta)

  zm <- .zmech_var(vc, .ref_beta, states = Snew)
  Slag <- rbind(NA, Snew[-nrow(Snew), , drop = FALSE])
  expect_equal(zm$z(alpha, s),
               drop(Slag %*% .hvec(alpha, .ref_beta, vc$H, 2L)))
  # support follows the states supplied, not the ones vc was fitted on.
  expect_equal(zm$support, stats::complete.cases(Slag))
  expect_error(.zmech_var(vc, .ref_beta, states = Snew[, 1:2, drop = FALSE]),
               "columns but H is")
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


# ---- perfect foresight (roadmap R-4, second mechanism) ----
#
# .zpf() sums realised d(ystar) forward instead of VAR forecasts. It shares
# .dweights() with the VAR route and nothing else, so the checks below pin
# it against three routes that do not overlap with it: two closed forms and
# the VAR closure itself.

test_that("perfect foresight on a constant target is sum_d times it", {
  alpha <- .lq_alpha(.ref_k, .ref_beta)$alpha
  s <- .scalars(alpha, .ref_beta)
  z <- .zpf(alpha, .ref_beta, rep(0.01, 500))

  ok <- seq_len(500 - z$horizon)
  expect_equal(z$Z[ok], rep(0.01 * s$sum_d, length(ok)))
})

test_that("perfect foresight on a geometric target matches its closed form", {
  # dys_t = rho^t makes sum_i d_i dys_{t+i} = dys_t * sum_i d_i rho^i, and
  # the latter is available in closed form as
  # A(1)A(beta) iota'(I-G)^{-1}(I - rho G)^{-1} iota -- the same identity
  # the DGP fixture uses to avoid a truncated forward sum.
  alpha <- .lq_alpha(.ref_k, .ref_beta)$alpha
  rho <- 0.5
  S1 <- .dgp_S(alpha, .ref_beta, rho)

  dys <- rho^(0:499)
  z <- .zpf(alpha, .ref_beta, dys)
  ok <- seq_len(500 - z$horizon)
  expect_equal(z$Z[ok], dys[ok] * S1)
})

test_that("perfect foresight equals the VAR route on a predictable path", {
  # The cross-mechanism check, and the one that pins the t-1 dating. On a
  # path the VAR forecasts EXACTLY, E_{t-1}[dystar_{t+i}] is the realised
  # value, so the two mechanisms must agree. mu + A cos(w t) satisfies
  # x_t = c + 2cos(w) x_{t-1} - x_{t-2}, which a VAR(2) with an intercept
  # recovers exactly.
  alpha <- .lq_alpha(.ref_k, .ref_beta)$alpha
  Tn <- 900
  w <- 2 * pi / 11
  dys <- 0.0104 + 0.004 * cos(w * seq_len(Tn))

  vc <- .var_companion(cbind(d.ystar = dys), p = 2)
  expect_lt(max(abs(vc$resid)), 1e-12)          # the VAR really is exact

  Slag <- rbind(NA, vc$states[-Tn, , drop = FALSE])
  Zvar <- drop(Slag %*% .hvec(alpha, .ref_beta, vc$H, sel = 2L))
  z <- .zpf(alpha, .ref_beta, dys)

  ok <- which(is.finite(Zvar) & is.finite(z$Z))
  expect_gt(length(ok), 700)
  # 1e-9, not machine precision: the VAR route collapses the sum in closed
  # form while .zpf() truncates it at .TOL_ZPF = 1e-10 relative. The gap
  # measured here is 1.1e-10, i.e. the truncation and nothing else.
  expect_lt(max(abs(Zvar[ok] - z$Z[ok]) / abs(Zvar[ok])), 1e-9)
})

test_that("the truncated tail is returned as NA, not padded", {
  alpha <- .lq_alpha(.ref_k, .ref_beta)$alpha
  z <- .zpf(alpha, .ref_beta, rep(0.01, 500))

  expect_equal(sum(is.na(z$Z)), z$horizon)
  expect_equal(z$n_lost, z$horizon)
  expect_true(all(is.na(utils::tail(z$Z, z$horizon))))
  expect_false(any(is.na(utils::head(z$Z, 500 - z$horizon))))
  # A series with no future at all yields all NA rather than a short sum.
  expect_true(all(is.na(.zpf(alpha, .ref_beta, rep(0.01, 20))$Z)))
})

test_that("perfect foresight refuses rather than truncating blind", {
  # rho(G) >= 1: the forward sum does not converge.       [INVARIANT 8]
  bad <- c(0, -1.5)
  expect_gt(.scalars(bad, 1)$rhoG, 1)
  expect_null(.zpf_horizon(bad, 1))
  expect_null(.zpf(bad, 1, rep(0.01, 500)))

  # Convergent, but not to tolerance within the terms allowed. Returning a
  # horizon would claim a bounded remainder that has not been established.
  alpha <- .lq_alpha(.ref_k, .ref_beta)$alpha
  expect_null(.zpf_horizon(alpha, .ref_beta, tol = 1e-12, max_h = 5L))

  # A tolerance below machine epsilon is refused, not granted. Once the
  # terms underflow relative to the running total the measured remainder is
  # exactly 0, which would certify any tolerance at all on rounding alone.
  expect_null(.zpf(alpha, .ref_beta, rep(0.01, 500), tol = 1e-300))
  expect_null(.zpf_horizon(alpha, .ref_beta, tol = 0))
  expect_false(is.null(.zpf_horizon(alpha, .ref_beta, tol = .EPS_ZPF)))
})

# ---- the mechanism seam (roadmap R-4) ----
#
# recm_estimate() builds a mechanism once and the residual function consumes
# its Z. These assert the contract in expectations.R for every mechanism
# that exists, so a third one is checked without a new test being written --
# and so the seam is exercised by something other than its only production
# caller.

.mechanisms <- function() {
  dys <- c(NA, diff(dgp_pac$ystar))
  vc <- .var_companion(cbind(d.ystar = dys), p = 2)
  list(
    var = .zmech_var(vc, 0.995),
    # 130 clears the 119 the reference calibration needs at beta = 0.995.
    perfect = .zmech_pf(dys, 0.995, horizon = 130L)
  )
}

test_that("every expectations mechanism satisfies the contract", {
  alpha <- .lq_alpha(.ref_k, 0.995)$alpha
  s <- .scalars(alpha, 0.995)
  Tn <- nrow(dgp_pac)

  for (nm in names(.mechanisms())) {
    mech <- .mechanisms()[[nm]]
    expect_true(is.character(mech$name) && length(mech$name) == 1L, info = nm)
    expect_true(is.logical(mech$support), info = nm)
    expect_length(mech$support, Tn)
    expect_false(anyNA(mech$support), info = nm)
    expect_gt(sum(mech$support), 0)
    expect_true(is.function(mech$z), info = nm)

    Z <- mech$z(alpha, s)
    expect_true(is.numeric(Z), info = nm)
    expect_length(Z, Tn)
    # The contract is that Z is available everywhere `support` says it is.
    expect_false(anyNA(Z[mech$support]), info = nm)
  }
})

test_that("mechanism support does not move with alpha", {
  # The load-bearing half of the contract: the sample is fixed before the
  # optimiser runs, so a support that followed theta would have the
  # criterion comparing SSRs computed on different observations.
  beta <- 0.995
  a1 <- .lq_alpha(c(1, 5), beta)$alpha
  a2 <- .lq_alpha(.ref_k, beta)$alpha
  expect_gt(.scalars(a2, beta)$rhoG, .scalars(a1, beta)$rhoG)

  for (nm in names(.mechanisms())) {
    mech <- .mechanisms()[[nm]]
    for (al in list(a1, a2)) {
      Z <- mech$z(al, .scalars(al, beta))
      expect_false(is.null(Z), info = nm)
      expect_false(anyNA(Z[mech$support]), info = nm)
    }
  }
})

test_that("the VAR mechanism reproduces the construction it replaced", {
  # estimate.R and diagnostics.R both built this by hand. Pin the mechanism
  # against that exact expression, so the consolidation is verified rather
  # than assumed.
  beta <- 0.995
  alpha <- .lq_alpha(.ref_k, beta)$alpha
  dys <- c(NA, diff(dgp_pac$ystar))
  vc <- .var_companion(cbind(d.ystar = dys), p = 2)
  Tn <- nrow(vc$states)

  Slag <- rbind(NA, vc$states[-Tn, , drop = FALSE])
  by_hand <- drop(Slag %*% .hvec(alpha, beta, vc$H, sel = 2L))
  expect_equal(.zmech_var(vc, beta)$z(alpha, .scalars(alpha, beta)), by_hand)
  # and the sample it implies is the one complete.cases(Slag) gave
  expect_equal(.zmech_var(vc, beta)$support, stats::complete.cases(Slag))
})

test_that("a mechanism refuses an alpha it cannot serve", {
  beta <- 0.995
  dys <- c(NA, diff(dgp_pac$ystar))

  # Perfect foresight at a horizon too short for this alpha: refused, not
  # truncated early, because the remainder would then be unbounded.
  alpha <- .lq_alpha(.ref_k, beta)$alpha
  need <- .zpf_horizon(alpha, beta)
  expect_gt(need, 20)
  s <- .scalars(alpha, beta)
  expect_null(.zmech_pf(dys, beta, horizon = 20L)$z(alpha, s))
  expect_false(is.null(.zmech_pf(dys, beta, horizon = need)$z(alpha, s)))

  # The VAR mechanism refuses on its own condition, rho(G)rho(H) >= 1.
  vc <- .var_companion(cbind(d.ystar = dys), p = 2)
  bad <- c(0, -1.5)
  expect_gt(.scalars(bad, beta)$rhoG, 1)
  expect_null(.zmech_var(vc, beta)$z(bad, .scalars(bad, beta)))
})

test_that("the horizon is where the exact remainder falls below tol", {
  alpha <- .lq_alpha(.ref_k, .ref_beta)$alpha
  s <- .scalars(alpha, .ref_beta)
  hh <- .zpf_horizon(alpha, .ref_beta)
  d <- .dweights(alpha, .ref_beta, horizon = hh)$d

  # At the horizon the remainder is inside tolerance, and one term earlier
  # it is not -- so the horizon is the first such H, not merely some H.
  expect_lte(abs(s$sum_d - sum(d)) / abs(s$sum_d), .TOL_ZPF)
  expect_gt(abs(s$sum_d - sum(d[-length(d)])) / abs(s$sum_d), .TOL_ZPF)

  # A looser tolerance must not need a longer horizon.
  expect_lte(.zpf_horizon(alpha, .ref_beta, tol = 1e-6), hh)
})
