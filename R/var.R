## The auxiliary model: a univariate autoregression in the first difference of
## the target, as in the FRB/US PAC basics note. The target is treated as I(1),
## which is Dynare's "dd" case in hVectors.m.
##
## The companion state is z[t] = (1, dystar[t], ..., dystar[t-p+1])' with
## z[t] = Phi z[t-1]. The leading 1 carries the intercept, so no demeaning is
## required, and the drift of the target is transmitted to the forward term
## automatically.
##
## Phi is estimated once and held fixed. Re-estimating it inside the zig-zag
## would make the fixed point chase noise in the auxiliary model.

# Ordinary least squares on an autoregression with a constant, over a fixed
# effective sample so that information criteria across orders are comparable.
ar_fit_ols <- function(x, p, start) {
  n <- length(x)
  idx <- start:n
  design <- matrix(1, length(idx), p + 1L)
  for (j in seq_len(p)) {
    design[, j + 1L] <- x[idx - j]
  }
  y <- x[idx]
  fit <- stats::lm.fit(design, y)
  rss <- sum(fit$residuals^2)
  nobs <- length(idx)
  list(coef = fit$coefficients, rss = rss, nobs = nobs,
       sigma2 = rss / (nobs - p - 1L),
       aic = nobs * log(rss / nobs) + 2 * (p + 1L))
}

# Fit the auxiliary autoregression, choosing the order by AIC.
#
# The order is not an argument of recm(), so it has to be selected. Candidate
# orders share the sample implied by p_max to keep the AIC comparison honest;
# the retained order is then re-estimated on its own maximal sample.
fit_aux_var <- function(dystar, p_max = NULL) {
  x <- dystar[!is.na(dystar)]
  n <- length(x)
  if (is.null(p_max)) {
    p_max <- max(1L, min(8L, n %/% 10L))
  }
  p_max <- max(1L, min(as.integer(p_max), (n - 2L) %/% 2L))
  if (n < p_max + 3L) {
    stop(
      "not enough observations on the target to fit the auxiliary ",
      "autoregression: ", n, " differences available.",
      call. = FALSE
    )
  }

  aic <- vapply(
    seq_len(p_max),
    function(p) ar_fit_ols(x, p, p_max + 1L)$aic,
    numeric(1)
  )
  p <- which.min(aic)
  fit <- ar_fit_ols(x, p, p + 1L)

  const <- unname(fit$coef[1L])
  ar <- unname(fit$coef[-1L])
  roots <- if (all(abs(ar) < 1e-12)) complex(0) else polyroot(c(1, -ar))
  if (length(roots) && any(Mod(roots) <= 1)) {
    warning(
      "the auxiliary autoregression in the differenced target has a root on ",
      "or inside the unit circle; the forward term may not converge.",
      call. = FALSE
    )
  }

  nz <- p + 1L
  phi <- matrix(0, nz, nz)
  phi[1L, 1L] <- 1
  phi[2L, 1L] <- const
  phi[2L, 2:nz] <- ar
  if (p > 1L) {
    phi[3:nz, 2:p] <- diag(p - 1L)
  }

  ev <- matrix(0, nz, 1L)
  ev[2L, 1L] <- 1

  list(order = p, const = const, ar = ar, sigma2 = fit$sigma2,
       nobs = fit$nobs, aic = stats::setNames(aic, seq_len(p_max)),
       companion = phi, selector = ev, roots = roots,
       rho = max(Mod(eigen(phi, only.values = TRUE)$values)),
       mean = if (sum(ar) < 1) const / (1 - sum(ar)) else NA_real_)
}

# The state z[t-1] for every row of the sample, aligned with the equation's
# time index: row t holds (1, dystar[t-1], ..., dystar[t-p]).
aux_state_lagged <- function(dystar, p) {
  n <- length(dystar)
  z <- matrix(NA_real_, n, p + 1L)
  z[, 1L] <- 1
  for (j in seq_len(p)) {
    z[(j + 1L):n, j + 1L] <- dystar[seq_len(n - j)]
  }
  z
}

# Extend the realised path of dystar with autoregressive forecasts, which the
# model consistent expectations recursion needs beyond the sample end.
aux_forecast_path <- function(dystar, aux, n_ahead) {
  obs <- dystar[!is.na(dystar)]
  p <- aux$order
  z <- c(1, rev(utils::tail(obs, p)))
  out <- numeric(n_ahead)
  for (i in seq_len(n_ahead)) {
    z <- as.numeric(aux$companion %*% z)
    out[i] <- z[2L]
  }
  out
}
