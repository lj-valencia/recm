## Shared fixtures and data generating processes for the rational error
## correction model.
##
## The two expectation options need different data generating processes and
## are not interchangeable. Under "var" the agent forms expectations from the
## auxiliary autoregression using information dated t-1, so the equation error
## is a predetermined adjustment cost shock and iterative OLS is consistent.
## Under "mce" the forward term is evaluated on the realised path; the equation
## error then contains forecast errors that are correlated with y[t-1], and
## therefore with the error correction term, so a series generated one way
## must not be used to test the other branch.

# The companion matrix of a stationary AR(p) with a constant, matching the
# layout fit_aux_var() builds.
make_companion <- function(const, ar) {
  p <- length(ar)
  nz <- p + 1L
  phi <- matrix(0, nz, nz)
  phi[1L, 1L] <- 1
  phi[2L, 1L] <- const
  phi[2L, 2:nz] <- ar
  if (p > 1L) {
    phi[3:nz, 2:p] <- diag(p - 1L)
  }
  phi
}

# An AR(p) path for the differenced target, and the target itself.
simulate_target <- function(nn, ar, const, sd_target) {
  p <- length(ar)
  dystar <- numeric(nn)
  for (t in (p + 1L):nn) {
    dystar[t] <- const + sum(ar * dystar[(t - 1L):(t - p)]) +
      stats::rnorm(1, 0, sd_target)
  }
  list(dystar = dystar, ystar = 100 + cumsum(dystar))
}

# Simulate y from the rational error correction decision rule, so the
# coefficients used here are exactly the ones estimation must recover.
simulate_recm <- function(n, a, beta, ar, const,
                          expectations = c("var", "mce"),
                          sd_target = 0.3, sd_eq = 0.05, burn = 200L,
                          extra = NULL, seed = NULL) {
  expectations <- match.arg(expectations)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  nn <- n + burn
  p <- length(ar)
  m <- length(a)
  tgt <- simulate_target(nn, ar, const, sd_target)
  dystar <- tgt$dystar
  ystar <- tgt$ystar

  alg <- recm:::pac_algebra(a, beta)
  if (expectations == "var") {
    phi <- make_companion(const, ar)
    ev <- matrix(0, p + 1L, 1L)
    ev[2L, 1L] <- 1
    h <- recm:::h_vector(alg, phi, ev)
    forward <- function(t) sum(h * c(1, dystar[(t - 1L):(t - p)]))
  } else {
    d <- recm:::lead_weights(alg, 600L)
    forward <- function(t) {
      horizon <- min(600L, nn - t)
      sum(d[seq_len(horizon + 1L)] * dystar[t + (0:horizon)])
    }
  }

  w <- if (is.null(extra)) NULL else stats::rnorm(nn, 0, 1)
  start <- p + m + 2L
  y <- numeric(nn)
  dy <- numeric(nn)
  y[seq_len(start - 1L)] <- ystar[seq_len(start - 1L)]
  for (t in start:nn) {
    step <- a[1L] * (ystar[t - 1L] - y[t - 1L]) + forward(t) +
      stats::rnorm(1, 0, sd_eq)
    if (m > 1L) {
      step <- step + sum(a[-1L] * dy[(t - 1L):(t - m + 1L)])
    }
    if (!is.null(w)) {
      step <- step + extra * w[t]
    }
    dy[t] <- step
    y[t] <- y[t - 1L] + step
  }

  idx <- (burn + 1L):nn
  out <- data.frame(y = y[idx], ystar = ystar[idx])
  if (!is.null(w)) {
    out$shock <- w[idx]
  }
  out
}
