# ---- polynomial arithmetic ----

# Multiply two polynomials given as ascending-power coefficient vectors.
# Written out rather than using convolve(), which goes through the FFT and
# introduces noise of order 1e-16 into a routine (factorise_spectral) that
# is already the ill-conditioned one.
.polymul <- function(p, q) {
  mode <- if (is.complex(p) || is.complex(q)) "complex" else "numeric"
  out <- vector(mode, length(p) + length(q) - 1L)
  for (i in seq_along(p)) {
    idx <- i + seq_along(q) - 1L
    out[idx] <- out[idx] + p[i] * q
  }
  out
}


# ---- numerical differentiation ----

# Forward-difference Jacobian with a relative step. Used for every
# delta-method standard error in the package, so the functions it is handed
# must be continuous -- see .halflife(), which interpolates for exactly this
# reason.
.jacnum <- function(f, x, h = 1e-5) {
  f0 <- f(x)
  out <- matrix(NA_real_, length(f0), length(x))
  for (j in seq_along(x)) {
    hh <- h * max(1, abs(x[j]))
    xp <- x
    xp[j] <- xp[j] + hh
    out[, j] <- (f(xp) - f0) / hh
  }
  out
}


# ---- covariance ----

# Newey-West HAC covariance of a matrix of moment contributions, Bartlett
# kernel. lags = NULL uses the usual floor(4 * (T/100)^(2/9)).
.nw <- function(M, lags = NULL) {
  M <- as.matrix(M)
  Tn <- nrow(M)
  if (is.null(lags)) lags <- floor(4 * (Tn / 100)^(2 / 9))
  S <- crossprod(M) / Tn
  for (l in seq_len(lags)) {
    w <- 1 - l / (lags + 1)
    G <- crossprod(M[(l + 1):Tn, , drop = FALSE],
                   M[1:(Tn - l), , drop = FALSE]) / Tn
    S <- S + w * (G + t(G))
  }
  S
}


# ---- data shaping ----

# Lag a vector, padding the head with NA. j = 0 returns x unchanged.
.lagv <- function(v, j) {
  if (j <= 0) return(v)
  c(rep(NA_real_, j), v[seq_len(length(v) - j)])
}

# Resolve a character vector or one-sided formula to a model matrix with no
# intercept. Returns NULL for a NULL spec so callers can test for absence.
#
# na.pass is essential, not cosmetic. model.matrix()'s default na.omit drops
# incomplete ROWS, returning a matrix shorter than `data` with no warning,
# which then misaligns against every other series in the caller. Any lagged
# column has leading NAs -- and the GMM dating rule tells users to supply
# exactly that. Row alignment is the caller's contract; NA handling is done
# once, by the complete.cases() sample selection in recm_estimate().
.getcols <- function(spec, data) {
  if (is.null(spec)) return(NULL)
  if (is.character(spec)) spec <- stats::reformulate(spec)
  spec <- stats::update(spec, ~ . + 0)
  mf <- stats::model.frame(spec, data, na.action = stats::na.pass)
  out <- stats::model.matrix(spec, mf)
  if (nrow(out) != nrow(data)) {
    stop("instrument/regressor specification returned ", nrow(out),
         " rows for ", nrow(data), " rows of `data`")
  }
  out
}
