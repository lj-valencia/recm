# Pure functions of alpha, beta and m. This file knows nothing about data,
# VARs or estimation -- everything here is testable without loading a
# dataset. See docs/01-architecture.md, "Module boundaries".
#
# Implements docs/02-math-spec.md sections 4, 5 and 7.

# Round-trip error for .alpha_to_a() / .a_to_alpha() must stay below this.
# INVARIANT 1, CLAUDE.md section 4.
.TOL_ROUNDTRIP <- 1e-12


# ---- reduced form maps ----

# alpha -> (a_0, a_1 .. a_{m-1}).
#
#   a_0 = -A(1) = -(1 + sum_i alpha_i)                  [INVARIANT 2]
#   a_i = sum_{j>i} alpha_j            i = 1 .. m-1
#
# Note the minus sign on a_0. The FRB/US "PAC Basics" note prints
# a_0 = A(1), which is inconsistent with its own equations (4) and (5).
# a_0 < 0 for a stable equation, hence A(1) > 0.
.alpha_to_a <- function(alpha) {
  m <- length(alpha)
  a <- numeric(m)
  a[1] <- -(1 + sum(alpha))
  if (m > 1) {
    for (i in seq_len(m - 1)) a[i + 1] <- sum(alpha[(i + 1):m])
  }
  a
}

# The exact inverse of .alpha_to_a(), used to go from a target decision rule
# back to costs:
#
#   c_0 = -a_0 - 1
#   c_i = a_i                          i = 1 .. m-1
#   alpha_i = c_{i-1} - c_i            i = 1 .. m-1
#   alpha_m = c_{m-1}
#
# Round-trip error must stay below .TOL_ROUNDTRIP.       [INVARIANT 1]
.a_to_alpha <- function(a) {
  m <- length(a)
  cc <- numeric(m)
  cc[1] <- -a[1] - 1
  if (m > 1) cc[2:m] <- a[2:m]
  alpha <- numeric(m)
  if (m > 1) {
    for (i in seq_len(m - 1)) alpha[i] <- cc[i] - cc[i + 1]
  }
  alpha[m] <- cc[m]
  alpha
}


# ---- companion form of A(beta F) ----

# Companion matrix of A(beta F), bottom row
# (-alpha_m beta^m, ..., -alpha_1 beta).
.G_mat <- function(alpha, beta) {
  m <- length(alpha)
  phi <- alpha * beta^seq_len(m)
  G <- matrix(0, m, m)
  if (m > 1) G[1:(m - 1), 2:m] <- diag(m - 1)
  G[m, ] <- -rev(phi)
  G
}


# ---- scalars ----

# A(1), A(beta), G, and the closed-form totals
#
#   sum_i h_i = A(1)A(beta) iota'(I-G)^{-1} iota
#   sum_i d_i = A(1)A(beta) iota'(I-G)^{-2} iota
#
# Returns NULL rather than erroring when (I-G) is singular, so the
# optimiser can reject a trial theta cheaply.
.scalars <- function(alpha, beta) {
  if (any(!is.finite(alpha))) return(NULL)
  m <- length(alpha)
  A1 <- 1 + sum(alpha)
  Ab <- 1 + sum(alpha * beta^seq_len(m))
  G  <- .G_mat(alpha, beta)
  im <- c(rep(0, m - 1), 1)
  IGi <- tryCatch(solve(diag(m) - G), error = function(e) NULL)
  if (is.null(IGi)) return(NULL)
  list(A1 = A1, Ab = Ab, G = G, im = im, IGi = IGi,
       sum_h = A1 * Ab * drop(im %*% IGi %*% im),
       sum_d = A1 * Ab * drop(im %*% IGi %*% IGi %*% im),
       rhoG  = max(Mod(eigen(G, only.values = TRUE)$values)))
}


# ---- forward weights ----

# Scalar lead weights, indexed from i = 0:
#
#   h_i = A(1)A(beta) iota' G^i iota                weights on ystar
#   d_i = A(1)A(beta) iota' (I-G)^{-1} G^i iota     weights on d.ystar
#
# The h_i equal A(1)A(beta) times the i-th coefficient of the power series
# 1 / A(beta F), obtainable by polynomial long division.  [INVARIANT 5]
.dweights <- function(alpha, beta, horizon = 100L) {
  s <- .scalars(alpha, beta)
  if (is.null(s)) return(NULL)
  m <- length(alpha)
  Gi <- diag(m)
  h <- numeric(horizon + 1)
  d <- numeric(horizon + 1)
  for (i in 0:horizon) {
    h[i + 1] <- s$A1 * s$Ab * drop(s$im %*% Gi %*% s$im)
    d[i + 1] <- s$A1 * s$Ab * drop(s$im %*% s$IGi %*% Gi %*% s$im)
    Gi <- Gi %*% s$G
  }
  list(h = h, d = d)
}


# ---- summaries of the lag distribution ----

# Response of the (y - ystar) gap to a unit shock, and the interpolated
# half-life.
#
# The 0.5 crossing is interpolated LINEARLY and must stay that way. An
# integer step has zero numerical derivative, so the delta-method standard
# error comes back as exactly 0 -- which is what happened before the fix.
.halflife <- function(a, horizon = 400) {
  m <- length(a)
  gap <- numeric(horizon + 1)
  dy  <- numeric(horizon + 1)
  gap[1] <- 1
  for (tt in 2:(horizon + 1)) {
    lg <- numeric(max(m - 1, 0))
    if (m > 1) {
      for (j in seq_len(m - 1)) lg[j] <- if (tt - j >= 1) dy[tt - j] else 0
    }
    dy[tt] <- a[1] * gap[tt - 1] + if (m > 1) sum(a[-1] * lg) else 0
    gap[tt] <- gap[tt - 1] + dy[tt]
  }
  hit <- which(abs(gap) < 0.5)
  hl <- NA_real_
  if (length(hit)) {
    j <- hit[1]
    hl <- if (j > 1) {
      (j - 2) + (abs(gap[j - 1]) - 0.5) / (abs(gap[j - 1]) - abs(gap[j]))
    } else {
      0
    }
  }
  list(halflife = hl, path = gap)
}
