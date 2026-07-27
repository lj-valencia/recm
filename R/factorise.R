# Maps k -> alpha. Pure; does not know what alpha will be used for.
#
# Two routes are kept deliberately (docs/01-architecture.md, "Two solution
# routes"). They solve the same problem, and keeping both is what makes
# INVARIANT 6 testable. Do not delete the spectral route because it is
# "slower" or "redundant" -- it is the independent implementation that
# certifies the fast one.
#
# Implements docs/02-math-spec.md section 3.

# Above this order the spectral route is ill-conditioned and warns. At
# m = 20 the degree-2m coefficients span roughly 5e12 and the reciprocity
# check |zeta_in * zeta_out / beta - 1| degrades to about 1e-3.
.M_SPECTRAL_MAX <- 5L

# Riccati iteration controls. Convergence is measured on the FEEDBACK GAIN,
# not on P: at m = 20 with slowly decaying costs max|P| reaches 1e8, so a
# P-based threshold is really a scale-dependent one and the iteration
# stalls short of it. The gain is the quantity we actually want -- it IS
# alpha -- and it stays O(1). INVARIANT 6 requires agreement with the
# spectral route to 1e-6, so 1e-10 on the gain leaves four orders of
# margin.
#
# Do NOT loosen this to make slow cases certify. Value iteration converges
# linearly at a rate that approaches 1 for high m with slowly decaying
# costs, so the per-step change UNDERSTATES the remaining error by roughly
# 1/(1-rate). Measured at m = 20, psi = 0.79: a per-step change of 1.1e-07
# sat 4.5e-05 away from the converged gain -- a factor of 400. Cases that
# cannot certify are genuinely inaccurate and must be rejected, not
# accepted with a wider threshold.
.TOL_RICCATI   <- 1e-10
.MAXIT_RICCATI <- 20000L


# ---- Riccati / discounted LQ (production route, all m) ----

# The production route. With state x_t = (y_{t-1}, ..., y_{t-m})' and
# control u_t = y_t, each penalty (1-L)^j y_t is linear in (u_t, x_t), so
# the stage cost is v'Mv with M = sum_j k_j g_j g_j'. Partition M into R
# (on u), N (cross) and Q (on x) and iterate
#
#   S = R + beta B'PB
#   K = N + beta A'PB
#   P <- Q + beta A'PA - KK'/S
#
# to convergence. The feedback gain F = (N' + beta B'PA)/S satisfies
# u_t = -F x_t, so alpha = F.                            [INVARIANT 6]
#
# Preferred over factorise_spectral() at every m: eigenvalues of an m by m
# matrix are far better conditioned than the roots of a degree-2m
# polynomial whose coefficients span many orders of magnitude.
#
# Returns NULL on an inadmissible k rather than erroring, because the
# optimiser calls this.
.lq_alpha <- function(k, beta, tol = .TOL_RICCATI, maxit = .MAXIT_RICCATI) {
  if (any(!is.finite(k)) || any(k <= 0)) return(NULL)
  m <- length(k) - 1
  M <- matrix(0, m + 1, m + 1)
  for (j in 0:m) {
    g <- numeric(m + 1)
    for (l in 0:j) g[l + 1] <- (-1)^l * choose(j, l)
    M <- M + k[j + 1] * (g %o% g)
  }
  R  <- M[1, 1]
  Nv <- matrix(M[2:(m + 1), 1], ncol = 1)
  Q  <- M[-1, -1, drop = FALSE]
  A  <- matrix(0, m, m)
  if (m > 1) A[2:m, 1:(m - 1)] <- diag(m - 1)
  B <- matrix(c(1, rep(0, m - 1)), ncol = 1)

  P <- Q
  Fg <- drop((t(Nv) + beta * crossprod(B, P %*% A)) /
               as.numeric(R + beta * crossprod(B, P %*% B)))
  conv <- FALSE
  for (it in seq_len(maxit)) {
    S  <- as.numeric(R + beta * crossprod(B, P %*% B))
    Kv <- Nv + beta * crossprod(A, P %*% B)
    P  <- Q + beta * crossprod(A, P %*% A) - tcrossprod(Kv) / S
    Sn <- as.numeric(R + beta * crossprod(B, P %*% B))
    Fn <- drop((t(Nv) + beta * crossprod(B, P %*% A)) / Sn)
    if (any(!is.finite(Fn))) return(NULL)
    if (max(abs(Fn - Fg)) < tol * max(1, max(abs(Fn)))) {
      Fg <- Fn
      conv <- TRUE
      break
    }
    Fg <- Fn
  }
  if (any(!is.finite(Fg))) return(NULL)

  # Closed-loop eigenvalues of A - BF are the reciprocals of the roots of
  # A(z) and must lie strictly inside the unit circle.
  maxeig <- max(Mod(eigen(A - B %*% matrix(Fg, 1),
                          only.values = TRUE)$values))
  list(alpha = Fg, P = P, converged = conv, iterations = it,
       maxeig = maxeig)
}


# ---- spectral factorisation (cross-check, m <= 5) ----

# The independent cross-check route. Setting L = z, F = 1/z, z^m P(q(z)) is
# an ordinary polynomial of degree 2m whose roots come in pairs
# (zeta, beta/zeta) -- exactly m outside the unit circle and m inside. Keep
# the outside set and build A(z) = prod_i (1 - z/zeta_i).
#
# Since z*q(z) = -beta + (1+beta)z - z^2, the coefficient vector of
# z^m P(q(z)) is assembled by convolving that quadratic j times and
# shifting by m - j.
#
# On the normalisation constant
# -----------------------------
# Matching leading coefficients gives
#
#     z^m P(q(z)) = cc * z^m A(beta/z) A(z),   cc = k_m * prod_i zeta_i
#
# and evaluating at z = 1 (where q(1) = 0) gives k_0 = cc * A(1) A(beta).
# So A(1)A(beta) = k_0 / cc: the identity k_0 = A(1)A(beta) holds only
# under the cc = 1 normalisation, which is NOT the k_0 = 1 normalisation
# this package estimates in. `k_norm` is the input rescaled to cc = 1, and
# it is `k_norm[1]` -- not `k[1]` -- that equals A(1)A(beta).
# See docs/02-math-spec.md section 3 and CLAUDE.md INVARIANT 3.
#
# Returns NULL if the root split is not m/m.
factorise_spectral <- function(k, beta) {
  if (any(!is.finite(k)) || any(k <= 0)) return(NULL)
  m <- length(k) - 1
  if (m > .M_SPECTRAL_MAX) {
    warning("factorise_spectral() at m = ", m, " > ", .M_SPECTRAL_MAX,
            ": the degree-", 2 * m, " polynomial is ill-conditioned and ",
            "the root pairing degrades. Use .lq_alpha() instead.",
            call. = FALSE)
  }

  quad <- c(-beta, 1 + beta, -1)          # z * q(z), ascending powers
  cf <- numeric(2 * m + 1)
  for (j in 0:m) {
    pj <- 1
    if (j > 0) for (i in seq_len(j)) pj <- .polymul(pj, quad)
    idx <- (m - j) + seq_along(pj)
    cf[idx] <- cf[idx] + k[j + 1] * pj
  }

  rts <- polyroot(cf)
  zo <- rts[Mod(rts) > 1]                 # outside the unit circle
  zi <- rts[Mod(rts) <= 1]
  if (length(zo) != m) return(NULL)

  # Reciprocity: each outside root should pair with an inside root at
  # beta/zeta. This is the conditioning diagnostic -- it stays near 1e-12
  # for small m and degrades to about 1e-3 by m = 20.
  recip <- 0
  for (z0 in zo) {
    partner <- zi[which.min(Mod(zi - beta / z0))]
    recip <- max(recip, Mod(partner * z0 / beta - 1))
  }

  Az <- 1
  for (z0 in zo) Az <- .polymul(Az, c(1, -1 / z0))
  alpha <- Re(Az[-1])

  cc <- Re(k[m + 1] * prod(zo))
  list(alpha = alpha, roots = zo, cc = cc, k_norm = k / cc,
       reciprocity = recip)
}
