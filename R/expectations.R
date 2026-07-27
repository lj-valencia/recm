# Owns the expectations mechanisms, of which the auxiliary VAR is the
# first. This is the ONLY module that knows the state vector layout. If
# that layout ever changes, it changes here and nowhere else -- see
# docs/01-architecture.md, "Module boundaries".
#
# State vector, fixed by .var_companion():
#
#   z_t = ( 1, X_t, X_{t-1}, ..., X_{t-p+1} )'     length n_z = 1 + k*p
#
# `states` stores these ROW-wise: one row per time period, n_z columns, so
# that Z = states %*% h is the natural call in estimate.R.
#
# The leading constant carries the VAR intercept, so no demeaning is
# required anywhere in the package. H[1,1] = 1, which puts a unit
# eigenvalue in the companion matrix -- expected, not a bug. Report
# rho_dyn (excluding the constant block) for display.
#
# d(ystar) is ALWAYS column 1 of X, hence position 2 of z. Several call
# sites pass sel = 2L on that assumption. If the target is ever allowed to
# move, introduce a named lookup rather than changing the constant.
#
# Implements docs/02-math-spec.md section 6.

# rho(G) * rho(H) must be below this for the forward sum to converge.
# INVARIANT 8. Because rho(H) is always 1 (the constant block), this
# reduces to rho(G) < 1 in practice, which is the correct conservative
# check.
.TOL_SPECTRAL <- 1 - 1e-8


# ---- auxiliary VAR ----

# The state matrix: row t is ( 1, X_t, X_{t-1}, ..., X_{t-p+1} ), and the
# first p-1 rows are NA because they have no full lag history.
#
# This is the loop that encodes the layout, and it is split out of
# .var_companion() so that it can be reached WITHOUT re-estimating the VAR.
# predict.recmfit() needs exactly that on `newdata`: new states, the fitted
# H. It wrote the loop out again to get it, which put the layout in a second
# module and quietly disagreed with this one -- it filled rows the lag
# history does not reach and never blanked the first p-1. Nothing outside
# this file may write it a third time.
.var_states <- function(X, p) {
  X <- as.matrix(X)
  Tn <- nrow(X)
  k <- ncol(X)
  if (Tn < p) {
    stop("cannot build VAR states: ", Tn, " observation",
         if (Tn == 1) "" else "s", " for a lag order of ", p)
  }
  nz <- 1 + k * p
  S <- matrix(NA_real_, Tn, nz)
  S[, 1] <- 1
  for (j in 0:(p - 1)) {
    idx <- p:Tn
    S[idx, (2 + j * k):(1 + (j + 1) * k)] <- X[idx - j, , drop = FALSE]
  }
  if (p > 1) S[seq_len(p - 1), ] <- NA
  colnames(S) <- c("const", paste0(rep(colnames(X), p), ".l",
                                   rep(0:(p - 1), each = k)))
  S
}

# Estimate the VAR by least squares and build its constant-augmented
# companion form.
#
# The VAR is a nuisance parameter and is estimated ONCE, outside the
# optimiser loop. Re-estimating it inside would make the criterion jagged
# -- the VAR coefficients would jitter with each trial theta and which.min
# would pick numerical noise rather than a minimum. The cost is that
# reported standard errors condition on H; that is the generated-regressor
# problem and it is handled by recm_boot().
.var_companion <- function(X, p) {
  X <- as.matrix(X)
  Tn <- nrow(X)
  k <- ncol(X)
  lagmat <- do.call(cbind, lapply(seq_len(p), function(j) {
    rbind(matrix(NA, j, k), X[seq_len(Tn - j), , drop = FALSE])
  }))
  ok <- stats::complete.cases(cbind(X, lagmat))
  Rg <- cbind(1, lagmat[ok, , drop = FALSE])
  B  <- qr.solve(Rg, X[ok, , drop = FALSE])

  nz <- 1 + k * p
  Hm <- matrix(0, nz, nz)
  Hm[1, 1] <- 1
  Hm[2:(k + 1), ] <- t(B)
  if (p > 1) Hm[(k + 2):nz, 2:(nz - k)] <- diag(k * (p - 1))

  S <- .var_states(X, p)

  resid <- X[ok, , drop = FALSE] - Rg %*% B
  list(H = Hm, states = S, coef = B, resid = resid,
       Sigma = crossprod(resid) / nrow(resid),
       k = k, p = p, n = nz, n_z = nz,
       rho = max(Mod(eigen(Hm, only.values = TRUE)$values)),
       rho_dyn = max(Mod(eigen(Hm[-1, -1, drop = FALSE],
                               only.values = TRUE)$values)))
}


# ---- closed-form forward sum ----

# Closed-form h such that Z_t = h' z_{t-1}, collapsing the infinite forward
# sum with the Kronecker identity sum_i G^i (x) H^i = [I - G (x) H]^{-1}:
#
#   Z_t = A(1)A(beta) ( iota'(I-G)^{-1} (x) e_v' )
#         [I - G (x) H]^{-1} ( iota (x) H ) z_{t-1}
#
# No horizon truncation is required.                     [INVARIANT 4]
#
# Returns NULL when rho(G) rho(H) >= 1, so the forward sum does not
# converge. Never let an estimation routine proceed past this. [INVARIANT 8]
.hvec <- function(alpha, beta, Hm, sel = 2L) {
  s <- .scalars(alpha, beta)
  if (is.null(s)) return(NULL)
  rhoH <- max(Mod(eigen(Hm, only.values = TRUE)$values))
  if (!is.finite(s$rhoG * rhoH) || s$rhoG * rhoH >= .TOL_SPECTRAL) {
    return(NULL)
  }
  m <- length(alpha)
  nz <- nrow(Hm)
  ev <- numeric(nz)
  ev[sel] <- 1
  K <- tryCatch(solve(diag(m * nz) - kronecker(s$G, Hm)),
                error = function(e) NULL)
  if (is.null(K)) return(NULL)
  left <- drop(s$im %*% s$IGi)
  out <- drop(s$A1 * s$Ab *
                (kronecker(matrix(left, 1), matrix(ev, 1)) %*% K %*%
                   kronecker(matrix(s$im, ncol = 1), Hm)))
  if (any(!is.finite(out))) return(NULL)
  out
}


# ---- perfect-foresight forward sum ----

# Relative tail mass at which the perfect-foresight sum is truncated, and
# the horizon beyond which we give up looking for that tolerance. The tail
# is measured against the CLOSED-FORM total sum_i d_i from .scalars(), so
# this is an exact remainder rather than a guess at one.
.TOL_ZPF <- 1e-10
.MAX_H_ZPF <- 5000L

# Tightest tolerance the remainder can actually be held to. It is measured
# as a difference of accumulated sums, so once the terms underflow relative
# to the total the measured remainder is exactly 0 -- which would certify
# ANY tolerance, however small, on nothing but rounding. Asking for one
# below this is refused rather than granted for free.
.EPS_ZPF <- .Machine$double.eps

# Horizon H at which sum_{i>H} d_i has fallen below `tol` of sum_i d_i.
#
# NULL when the sum does not converge (rho(G) >= 1) or when `max_h` terms
# are not enough to reach `tol`. Both are refusals, not approximations: a
# caller that got a horizon back can rely on the truncation being below
# tolerance, which is the whole point of returning one.
.zpf_horizon <- function(alpha, beta, tol = .TOL_ZPF, max_h = .MAX_H_ZPF) {
  if (!is.finite(tol) || tol < .EPS_ZPF) return(NULL)
  s <- .scalars(alpha, beta)
  if (is.null(s)) return(NULL)
  if (!is.finite(s$rhoG) || s$rhoG >= .TOL_SPECTRAL) return(NULL)
  if (!is.finite(s$sum_d) || abs(s$sum_d) < .Machine$double.eps) return(NULL)
  dw <- .dweights(alpha, beta, horizon = max_h)
  if (is.null(dw)) return(NULL)
  # cumsum(d)[j] is sum_{i=0}^{j-1} d_i, so truncating after H = j-1 leaves
  # sum_d - cumsum(d)[j]. d_i oscillates in sign under complex roots, so
  # the remainder is taken directly rather than inferred from |d_H|.
  tail_rel <- abs(s$sum_d - cumsum(dw$d)) / abs(s$sum_d)
  hit <- which(is.finite(tail_rel) & tail_rel <= tol)
  if (!length(hit)) return(NULL)
  hit[1L] - 1L
}

# Perfect-foresight expectations: the second mechanism (roadmap R-4).
#
#   Z_t = sum_{i>=0} d_i D ystar_{t+i}
#
# which is the VAR closure of docs/02 section 6 with E_{t-1}[.] replaced by
# the REALISED path. Certainty equivalence is what makes this worth having:
# the analytic decision rule evaluated at this Z must reproduce the exact
# quadratic-program solution, so it cross-checks the whole
# alpha -> d_i -> Z chain against a route that shares none of it.
#
# Two ways it differs from .hvec(), both consequences of summing over
# realised data rather than over a companion form:
#
#   * No Kronecker collapse is available, so the sum is TRUNCATED. The
#     horizon comes from .zpf_horizon() and its remainder is bounded.
#   * The last H observations have no future left to sum over and are
#     returned as NA. They are not padded with the final value or with a
#     forecast -- padding would substitute a fabricated continuation for
#     the perfect foresight this mechanism claims, and NA lets the caller's
#     complete.cases() drop them, which is the honest cost.
#
# Convergence requires rho(G) < 1. That is INVARIANT 8 specialised to a
# mechanism carrying no H, and it is the same condition in practice, since
# rho(H) is always exactly 1 through the constant block.
#
# `dys` is the realised D ystar series, one element per period, aligned
# with the rows of `data` in estimate.R. Its own leading NA propagates.
# `horizon`, when supplied, fixes the truncation instead of letting it
# follow alpha. That is what a mechanism needs (see below): the sample must
# not move with theta. An alpha whose own required horizon exceeds the fixed
# one is REFUSED rather than truncated early -- the bound on the remainder
# is the only thing that makes the sum meaningful. A fixed horizon longer
# than required is harmless: the extra weights are the ones already known to
# be negligible.
.zpf <- function(alpha, beta, dys, tol = .TOL_ZPF, horizon = NULL) {
  hh <- .zpf_horizon(alpha, beta, tol)
  if (is.null(hh)) return(NULL)
  if (!is.null(horizon)) {
    horizon <- as.integer(horizon)
    if (is.na(horizon) || horizon < 0L || hh > horizon) return(NULL)
    hh <- horizon
  }
  dw <- .dweights(alpha, beta, horizon = hh)
  if (is.null(dw)) return(NULL)
  dd <- dw$d
  Tn <- length(dys)
  Z <- rep(NA_real_, Tn)
  # The obvious loop, kept obvious. It runs once per candidate alpha, not
  # once per observation per candidate, and a filter-based route would need
  # its own agreement test to earn the speed.
  if (Tn > hh) {
    for (tt in seq_len(Tn - hh)) Z[tt] <- sum(dd * dys[tt + 0:hh])
  }
  list(Z = Z, horizon = hh, n_lost = min(hh, Tn))
}


# ---- expectations mechanisms ----
#
# Extension point (docs/01-architecture.md, roadmap R-4). A mechanism is
# built ONCE in recm_estimate(), before the optimiser, and supplies the
# forward sum `Z` to the residual function. The residual function no longer
# builds `Z` itself, which is the whole point: it now consumes a series it
# does not know the provenance of.
#
# A mechanism is a list of three things:
#
#   name        character(1), reported by summary()
#   support     logical, length Tn. Rows this mechanism can supply Z for at
#               ANY admissible alpha.
#   z(alpha, s) numeric length Tn, or NULL if this alpha is inadmissible
#               FOR THIS MECHANISM. `s` is .scalars(alpha, beta), passed in
#               rather than recomputed.
#
# **`support` must not depend on theta.** This is the load-bearing part of
# the contract and the reason the seam is shaped this way rather than as a
# bare function. The estimation sample is fixed once, before the optimiser
# runs, so that every trial theta is scored on the same observations; a
# criterion compared across moving samples is not a criterion. The VAR
# mechanism satisfies this for free -- its NA pattern comes from the lag
# structure, not from alpha. Perfect foresight does not, because its horizon
# follows rho(G), which is why .zmech_pf() takes a FIXED horizon and refuses
# any alpha needing more instead of quietly shortening the sum.
#
# The admissibility condition is the mechanism's own, not a shared one. The
# VAR route needs rho(G) rho(H) < 1 [INVARIANT 8]; perfect foresight carries
# no H and needs rho(G) < 1 plus a horizon inside budget. Neither belongs in
# the residual function.

# VAR-based expectations: the original mechanism, Z_t = h' z_{t-1}.
#
# Owns the one-period lag of the state matrix, which used to be built at the
# call site in estimate.R and again in diagnostics.R.
#
# `states` defaults to the ones `vc` was estimated on and is overridden only
# to run the FITTED VAR over different observations -- predict.recmfit() on
# `newdata`. H, and therefore the forecast rule, still comes from `vc`:
# re-estimating it there would predict from a different model than the one
# that was fitted. `support` is recomputed from whatever states are supplied,
# which keeps it a property of the lag structure and not of theta.
.zmech_var <- function(vc, beta, sel = 2L, states = vc$states) {
  if (ncol(states) != nrow(vc$H)) {
    stop("state matrix has ", ncol(states), " columns but H is ",
         nrow(vc$H), " by ", nrow(vc$H))
  }
  Tn <- nrow(states)
  Slag <- rbind(NA, states[-Tn, , drop = FALSE])
  list(
    name = "var",
    support = stats::complete.cases(Slag),
    z = function(alpha, s) {
      if (!is.finite(s$rhoG * vc$rho) || s$rhoG * vc$rho >= 1) return(NULL)
      hv <- .hvec(alpha, beta, vc$H, sel)
      if (is.null(hv)) return(NULL)
      drop(Slag %*% hv)
    }
  )
}

# Perfect-foresight expectations at a FIXED horizon (roadmap R-4).
#
# Not reachable from recm_estimate() yet -- there is no
# `expectations_backend` argument, because how the horizon should be chosen
# is still open. See R-4. It is exercised through the mechanism contract
# test, which is what keeps this seam honest rather than notional.
#
# The fixed horizon buys a fixed sample and costs admissible parameter
# space: alphas whose forward sum has not decayed inside `horizon` are
# rejected. Choose it generously, and note that the last `horizon`
# observations are outside `support` and are lost.
.zmech_pf <- function(dys, beta, horizon, tol = .TOL_ZPF) {
  Tn <- length(dys)
  horizon <- as.integer(horizon)
  sup <- vapply(seq_len(Tn), function(tt) {
    tt + horizon <= Tn && !anyNA(dys[tt + 0:horizon])
  }, logical(1))
  list(
    name = "perfect",
    support = sup,
    z = function(alpha, s) {
      out <- .zpf(alpha, beta, dys, tol = tol, horizon = horizon)
      if (is.null(out)) NULL else out$Z
    }
  )
}
