# Owns the auxiliary VAR and the companion form. This is the ONLY module
# that knows the state vector layout. If that layout ever changes, it
# changes here and nowhere else. See docs/01-architecture.md.
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

  S <- matrix(NA_real_, Tn, nz)
  S[, 1] <- 1
  for (j in 0:(p - 1)) {
    idx <- p:Tn
    S[idx, (2 + j * k):(1 + (j + 1) * k)] <- X[idx - j, , drop = FALSE]
  }
  if (p > 1) S[seq_len(p - 1), ] <- NA
  colnames(S) <- c("const", paste0(rep(colnames(X), p), ".l",
                                   rep(0:(p - 1), each = k)))

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
