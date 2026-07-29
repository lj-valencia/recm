## The shared algebra of the rational error correction model: the map from the
## reduced-form coefficients a = (a0, a1, ..., a_{m-1}) to the stable lag
## polynomial A(L), the companion matrix G of A(beta F), the lead weights d_i,
## and back out to the structural adjustment costs k_j.
##
## Conventions follow Dynare's PAC implementation (a2alpha.m, buildGmatrix.m):
##
##   dy[t] = a0 * (ystar[t-1] - y[t-1]) + sum_{i=1}^{m-1} a[i] * dy[t-i] + Z[t]
##
## so a0 = A(1) = theta > 0. alpha is stored in ascending powers of L, that is
## A(z) = 1 + alpha[1] z + ... + alpha[m] z^m.

# Multiply two polynomials given as coefficient vectors in ascending powers.
poly_mul <- function(p, q) {
  r <- numeric(length(p) + length(q) - 1L)
  for (i in seq_along(p)) {
    idx <- i:(i + length(q) - 1L)
    r[idx] <- r[idx] + p[i] * q
  }
  r
}

# Laurent polynomials are a list(offset, coef) with coef in ascending powers
# starting at L^offset.
laurent_mul <- function(p, q) {
  list(offset = p$offset + q$offset, coef = poly_mul(p$coef, q$coef))
}

# Pad a Laurent polynomial to the fixed window L^-m .. L^m.
laurent_pad <- function(p, m) {
  out <- numeric(2L * m + 1L)
  first <- p$offset + m + 1L
  out[first:(first + length(p$coef) - 1L)] <- p$coef
  out
}

# Reduced-form a -> lag polynomial coefficients alpha. Dynare's a2alpha.m,
# with the m = 1 case (where a has no autoregressive element) treated as
# alpha_1 = a0 - 1 rather than as an indexing error.
a_to_alpha <- function(a) {
  m <- length(a)
  alpha <- numeric(m)
  alpha[m] <- a[m]
  if (m > 2L) {
    alpha[2:(m - 1L)] <- a[2:(m - 1L)] - a[3:m]
  }
  alpha[1L] <- a[1L] - (if (m > 1L) a[2L] else 0) - 1
  alpha
}

# The exact inverse: a0 = A(1), a_i = sum_{j>i} alpha_j.
alpha_to_a <- function(alpha) {
  m <- length(alpha)
  a <- numeric(m)
  a[1L] <- 1 + sum(alpha)
  if (m > 1L) {
    for (i in 2:m) {
      a[i] <- sum(alpha[i:m])
    }
  }
  a
}

# Companion matrix of A(beta F). Identical to Dynare's buildGmatrix.m.
g_matrix <- function(alpha, beta) {
  m <- length(alpha)
  g <- matrix(0, m, m)
  if (m > 1L) {
    g[1:(m - 1L), 2:m] <- diag(m - 1L)
  }
  g[m, ] <- -rev(alpha * beta^seq_len(m))
  g
}

# Everything downstream of a and beta, computed once per iterate.
#
# d_sum is the total forward loading sum_i d_i, lead_sum the weighted sum
# sum_i i * d_i used for the mean lead. Both are closed forms in G, so no
# horizon truncation enters anywhere.
pac_algebra <- function(a, beta) {
  m <- length(a)
  alpha <- a_to_alpha(a)
  ap <- c(1, alpha)
  a1 <- sum(ap)
  ab <- sum(ap * beta^(0:m))
  cc <- a1 * ab
  g <- g_matrix(alpha, beta)
  iot <- matrix(0, m, 1L)
  iot[m, 1L] <- 1
  igi <- solve(diag(m) - g)
  d_sum <- as.numeric(cc * t(iot) %*% igi %*% igi %*% iot)
  lead_sum <- as.numeric(cc * t(iot) %*% igi %*% g %*% igi %*% igi %*% iot)
  list(a = a, m = m, beta = beta, alpha = alpha, ap = ap,
       a1 = a1, ab = ab, cc = cc, g = g, igi = igi, iot = iot,
       d_sum = d_sum, lead_sum = lead_sum,
       rho_g = max(Mod(eigen(g, only.values = TRUE)$values)))
}

# Lead weights d_i = A(1)A(beta) iota' (I-G)^{-1} G^i iota, for i = 0..n.
lead_weights <- function(alg, n = 40L) {
  v <- alg$cc * (t(alg$iot) %*% alg$igi)
  out <- numeric(n + 1L)
  gp <- diag(alg$m)
  for (i in 0:n) {
    out[i + 1L] <- as.numeric(v %*% gp %*% alg$iot)
    gp <- gp %*% alg$g
  }
  out
}

# Admissibility: every root of A(z) strictly outside the unit circle, and
# rho(G) rho(Phi) < 1 so that the closed-form forward sum converges.
is_admissible <- function(alg, rho_phi, tol = 1e-8) {
  if (any(!is.finite(alg$alpha))) {
    return(FALSE)
  }
  if (alg$rho_g * rho_phi >= 1 - tol) {
    return(FALSE)
  }
  nz <- max(which(abs(alg$ap) > tol))
  if (nz == 1L) {
    return(TRUE)
  }
  rts <- polyroot(alg$ap[seq_len(nz)])
  all(Mod(rts) > 1 + tol)
}

# Structural adjustment costs from the lag polynomial.
#
# The Euler equation of the order-m problem is
#   [k_0 + k_1 Q + ... + k_m Q^m] y[t] = k_0 ystar[t],  Q = (1 - beta F)(1 - L)
# and it factorises as A(beta F) A(L) y[t] = A(1) A(beta) ystar[t]. Matching
# the two Laurent polynomials in L is linear in k, so a single least squares
# solve recovers the costs for any m. The relative costs b_j = k_j / k_0 are
# the normalisation k_0 = 1 used in the notes; k_0 must equal A(1)A(beta).
cost_params <- function(alpha, beta) {
  m <- length(alpha)
  ap <- c(1, alpha)
  target <- laurent_mul(
    list(offset = 0L, coef = ap),
    list(offset = -m, coef = rev(ap * beta^(0:m)))
  )
  qq <- list(offset = -1L, coef = c(-beta, 1 + beta, -1))
  mat <- matrix(0, 2L * m + 1L, m + 1L)
  cur <- list(offset = 0L, coef = 1)
  mat[, 1L] <- laurent_pad(cur, m)
  for (j in seq_len(m)) {
    cur <- laurent_mul(cur, qq)
    mat[, j + 1L] <- laurent_pad(cur, m)
  }
  rhs <- laurent_pad(target, m)
  k <- qr.solve(mat, rhs)
  # A reduced form estimated freely need not come from a convex adjustment
  # cost problem. `positive` reports whether it does; it is a diagnostic, not
  # an admissibility condition, because the zig-zag works in reduced-form
  # space and constraining it there would change the estimator.
  list(k = k, b = k[-1L] / k[1L], positive = all(k > 0),
       residual = max(abs(as.numeric(mat %*% k) - rhs)))
}

# Deterministic response of the error correction gap to a unit gap, with the
# target held flat: the path ystar - y would follow if the equation were run
# forward from a one-off displacement and nothing else ever happened.
#
# Returned including the initial gap, so the result has h + 1 elements and
# element i + 1 is the gap after i periods.
gap_path <- function(a, h) {
  m <- length(a)
  gap <- 1
  dy <- numeric(m)
  out <- numeric(h + 1L)
  out[1L] <- gap
  for (i in seq_len(h)) {
    step <- a[1L] * gap
    if (m > 1L) {
      step <- step + sum(a[-1L] * dy[seq_len(m - 1L)])
    }
    dy <- c(step, dy[-length(dy)])
    gap <- gap - step
    out[i + 1L] <- gap
  }
  out
}

# Half-life of the error correction gap, in periods, from that same response.
#
# Interpolated across the crossing rather than reported as the first integer
# that clears it: an integer crossing has zero derivative and would return a
# standard error of exactly zero under the delta method. The interpolation is
# log-linear where the gap is still positive, which is exact at m = 1, and
# falls back to linear if the response overshoots through zero.
half_life <- function(a, max_h = 1000L) {
  path <- gap_path(a, max_h)
  below <- which(path <= 0.5)
  if (!length(below)) {
    return(NA_real_)
  }
  # path[1] is the gap at h = 0, so index i is period i - 1.
  h <- below[1L] - 1L
  gap <- path[h + 1L]
  prev <- path[h]
  if (isTRUE(all.equal(prev, gap))) {
    return(h)
  }
  share <- if (gap > 0) {
    (log(prev) - log(0.5)) / (log(prev) - log(gap))
  } else {
    (prev - 0.5) / (prev - gap)
  }
  h - 1 + share
}
