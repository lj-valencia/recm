## Case III of the notes: iterative ordinary least squares, the FRB/US zig-zag,
## with the growth neutrality restriction imposed at each step.
##
## The estimator works in reduced-form space and never forms an objective
## function. Each sweep rebuilds the forward term at the current coefficients,
## moves it to the left-hand side as an offset - which is how the coefficient
## on Z is fixed at 1, and therefore how the cross-equation restriction is
## imposed - and runs one linear least squares step.

# Least squares, optionally subject to one linear equality restriction
# cmat %*% gamma == rhs, in closed form via the Lagrangian.
constrained_ls <- function(x, y, cmat = NULL, rhs = NULL) {
  xtx <- crossprod(x)
  xtx_inv <- tryCatch(solve(xtx), error = function(e) NULL)
  if (is.null(xtx_inv)) {
    stop(
      "the regressor matrix of the error correction equation is singular; ",
      "check for a constant or perfectly collinear columns in `data`.",
      call. = FALSE
    )
  }
  gamma <- as.numeric(xtx_inv %*% crossprod(x, y))
  if (is.null(cmat)) {
    return(list(gamma = gamma, xtx_inv = xtx_inv))
  }
  cxc <- cmat %*% xtx_inv %*% t(cmat)
  adj <- xtx_inv %*% t(cmat) %*% solve(cxc, rhs - cmat %*% gamma)
  list(gamma = gamma + as.numeric(adj), xtx_inv = xtx_inv)
}

# Newey-West long run variance of the moment contributions, with the
# Bartlett bandwidth floor(4 (T/100)^(2/9)).
hac_variance <- function(g) {
  tt <- nrow(g)
  lag_max <- floor(4 * (tt / 100)^(2 / 9))
  s <- crossprod(g) / tt
  if (lag_max >= 1L) {
    for (j in seq_len(lag_max)) {
      gj <- crossprod(
        g[(j + 1L):tt, , drop = FALSE],
        g[seq_len(tt - j), , drop = FALSE]
      ) / tt
      s <- s + (1 - j / (lag_max + 1L)) * (gj + t(gj))
    }
  }
  list(s = s, lag = lag_max)
}

# The zig-zag. `forward` is a closure mapping a pac_algebra object to the
# vector Z, so the "var" and "mce" options are interchangeable here.
iterative_ols <- function(dy, x, m, beta, forward, rho_phi,
                          restrict, tol = 1e-8, max_iter = 200L) {
  k <- ncol(x)
  n_w <- k - m

  admissible_a <- function(a) {
    alg <- tryCatch(pac_algebra(a, beta), error = function(e) NULL)
    if (is.null(alg) || !is_admissible(alg, rho_phi)) NULL else alg
  }

  start <- constrained_ls(x, dy)$gamma
  if (is.null(admissible_a(start[seq_len(m)]))) {
    start[seq_len(m)] <- c(0.2, numeric(m - 1L))
    if (is.null(admissible_a(start[seq_len(m)]))) {
      stop(
        "could not find an admissible starting value for the lag polynomial; ",
        "check `m` and `discount`.",
        call. = FALSE
      )
    }
  }

  gamma <- start
  alg <- admissible_a(gamma[seq_len(m)])
  converged <- FALSE
  crit <- NA_real_
  iter <- 0L

  for (iter in seq_len(max_iter)) {
    z <- forward(alg)
    offset <- dy - z

    cmat <- NULL
    rhs <- NULL
    if (restrict) {
      a_now <- gamma[seq_len(m)]
      grad <- gn_grad(a_now, beta)
      cmat <- matrix(c(grad, numeric(n_w)), nrow = 1L)
      rhs <- matrix(sum(grad * a_now) - gn_gap_alg(alg), 1L, 1L)
    }

    step <- constrained_ls(x, offset, cmat, rhs)
    proposal <- step$gamma

    lambda <- 1
    repeat {
      trial <- lambda * proposal + (1 - lambda) * gamma
      alg_try <- admissible_a(trial[seq_len(m)])
      if (!is.null(alg_try)) {
        break
      }
      lambda <- lambda / 2
      if (lambda < 1e-10) {
        stop(
          "no admissible damped step exists at iteration ", iter,
          "; the lag polynomial leaves the stable region.",
          call. = FALSE
        )
      }
    }

    crit <- max(abs(trial - gamma))
    gamma <- trial
    alg <- alg_try
    if (crit < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "iterative OLS did not converge in ", max_iter,
      " iterations (criterion ", format(crit, digits = 3), ").",
      call. = FALSE
    )
  }

  z <- forward(alg)
  resid <- dy - z - as.numeric(x %*% gamma)

  list(gamma = gamma, alg = alg, z = z, residuals = resid,
       ssr = sum(resid^2), iterations = iter, converged = converged,
       criterion = crit)
}

# Covariance of the iterative OLS estimator.
#
# The fixed point solves X'(dy - X gamma - Z(a)) = 0, which is just-identified
# GMM with the regressors as their own instruments. The moment condition
# depends on a through Z, so the OLS standard errors of the final sweep are
# wrong twice over: they omit dZ/da' and they condition on the estimated
# auxiliary model. The sandwich below repairs the first of those.
iterative_ols_vcov <- function(fit, dy, x, m, beta, forward, cmat = NULL,
                               step = 1e-6) {
  tt <- nrow(x)
  k <- ncol(x)
  gamma <- fit$gamma

  zgrad <- matrix(0, tt, k)
  for (j in seq_len(m)) {
    up <- gamma[seq_len(m)]
    dn <- up
    up[j] <- up[j] + step
    dn[j] <- dn[j] - step
    z_up <- forward(pac_algebra(up, beta))
    z_dn <- forward(pac_algebra(dn, beta))
    zgrad[, j] <- (z_up - z_dn) / (2 * step)
  }

  jac <- -crossprod(x, x + zgrad) / tt
  jac_inv <- tryCatch(solve(jac), error = function(e) NULL)
  if (is.null(jac_inv)) {
    return(list(vcov = matrix(NA_real_, k, k), hac_lag = NA_integer_))
  }

  hac <- hac_variance(x * fit$residuals)
  v <- jac_inv %*% hac$s %*% t(jac_inv) / tt

  if (!is.null(cmat)) {
    vc <- v %*% t(cmat)
    v <- v - vc %*% solve(cmat %*% vc, t(vc))
  }

  dimnames(v) <- list(colnames(x), colnames(x))
  list(vcov = v, hac_lag = hac$lag, z_gradient = zgrad)
}
