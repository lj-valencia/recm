## The forward term
##
##   Z[t] = sum_{i >= 0} d_i E[dystar[t+i]]
##
## under the two expectation options of the FRB/US note.
##
##   "var" - expectations formed from the auxiliary autoregression, information
##           set dated t-1. Z[t] = h' z[t-1], h in closed form, so no horizon
##           truncation enters. This is Dynare's hVectors.m, kind "dd".
##
##   "mce" - model consistent expectations, evaluated on the realised path of
##           the target (extended past the sample end with autoregressive
##           forecasts) by the finite-lead backward recursion.

# h' = A(1)A(beta) (iota'(I-G)^{-1} kron ev')
#        [I - G kron Phi]^{-1} (iota kron Phi)
#
# using sum_i G^i kron Phi^i = [I - G kron Phi]^{-1}. Requires
# rho(G) rho(Phi) < 1, which is_admissible() enforces upstream.
h_vector <- function(alg, phi, ev) {
  m <- alg$m
  nz <- nrow(phi)
  lhs <- kronecker(t(alg$iot) %*% alg$igi, t(ev))
  rhs <- kronecker(alg$iot, phi)
  kmat <- diag(m * nz) - kronecker(alg$g, phi)
  as.numeric(alg$cc * (lhs %*% solve(kmat, rhs)))
}

# Z[t] = h' z[t-1] for every row of zlag.
z_var <- function(alg, phi, ev, zlag) {
  as.numeric(zlag %*% h_vector(alg, phi, ev))
}

# The same forward term as an explicit truncated sum, Z[t] = sum_{i<=n} d_i
# ev' Phi^{i+1} z[t-1]. Kept because the closed form above is worth checking
# against it; the two must agree to machine precision.
z_var_truncated <- function(alg, phi, ev, zlag, n = 500L) {
  d <- lead_weights(alg, n)
  acc <- numeric(nrow(zlag))
  gp <- phi
  for (i in 0:n) {
    acc <- acc + d[i + 1L] * as.numeric(zlag %*% (t(gp) %*% ev))
    gp <- gp %*% phi
  }
  acc
}

# Model consistent expectations. Backward recursion on the realised path,
#
#   Z[t] = -sum_i alpha_i beta^i Z[t+i]
#          + A(1) [ dystar[t] + sum_{k=1}^{m-1} c_k dystar[t+k] ],
#   c_k = -sum_{j>k} alpha_j beta^j
#
# started from Z = 0 at the far end of the padded path. A constant path
# dystar = g must return Z = sum_i(d_i) * g, which is the check that ties this
# branch to the closed form above.
z_mce <- function(alg, path, n_out) {
  m <- alg$m
  beta <- alg$beta
  ab <- alg$alpha * beta^seq_len(m)
  ck <- if (m > 1L) {
    vapply(
      seq_len(m - 1L),
      function(k) -sum(alg$alpha[(k + 1L):m] * beta^((k + 1L):m)),
      numeric(1)
    )
  } else {
    numeric(0)
  }
  np <- length(path)
  padded <- c(path, rep(path[np], m))
  z <- numeric(np + m)
  for (t in np:1) {
    fwd <- if (m > 1L) sum(ck * padded[(t + 1L):(t + m - 1L)]) else 0
    z[t] <- -sum(ab * z[(t + 1L):(t + m)]) + alg$a1 * (padded[t] + fwd)
  }
  z[seq_len(n_out)]
}

# How far past the sample end the model consistent path must run for the
# terminal Z = 0 to be harmless.
mce_pad_length <- function(alg, tol = 1e-12, floor_n = 200L, cap = 5000L) {
  rho <- max(alg$rho_g, 1e-8)
  if (rho >= 1) {
    return(cap)
  }
  max(floor_n, min(cap, ceiling(log(tol) / log(rho))))
}
