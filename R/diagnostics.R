# Weak identification is the normal case for this estimator, not the
# exception. Everything in this file exists because point estimates alone
# are not a defensible output. See docs/02-math-spec.md section 10.

# Below this relative range of the criterion over the profiled grid, the
# parameter is effectively unidentified and a note fires. Given how
# routinely this triggers, recm_profile() should be treated as part of the
# standard output rather than an optional diagnostic.
.TOL_FLAT_PROFILE <- 0.02


# ---- profile ----

#' Profile the criterion in one parameter
#'
#' Re-optimises the remaining parameters at each point of a grid in one
#' parameter, using the `objective` closure stored on the fit. No
#' re-specification of the model is required.
#'
#' @param object an object of class `recmfit`.
#' @param which integer, index of the parameter to profile.
#' @param span numeric, half-width of the grid in estimated standard errors.
#' @param ngrid integer, number of grid points.
#' @param level numeric(1) confidence level for the profile interval.
#' @param x an object of class `recm_profile`.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class `recm_profile`, a list with `grid`, `value`,
#'   `a0`, `cutoff`, `estimate`, `level`, `name`, `interval` and
#'   `rel_range`.
#'
#' @details
#' `rel_range` is the range of the criterion over the grid, relative to its
#' minimum. Below 0.02 the parameter is effectively unidentified and the
#' print method says so. Expect it: at `m = 2` a fourfold change in
#' `k_1/k_0` cost 2.2% of SSR, and the true parameters sat 1.7% above the
#' optimum. This is inverse optimal control -- many cost configurations
#' rationalise nearly identical observed behaviour.
#'
#' @seealso [recm_boot()], [recm_estimate()]; `docs/02-math-spec.md`
#'   section 10.
#'
#' @examples
#' set.seed(1)
#' n <- 160
#' dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5,
#'                                 method = "recursive"))
#' dat <- data.frame(ystar = cumsum(dys))
#' dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))
#' fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
#'                      var_lags = 2, restarts = 1, quiet = TRUE)
#' recm_profile(fit, which = 1, ngrid = 7)
#'
#' @export
recm_profile <- function(object, which = 1, span = 3, ngrid = 25,
                         level = 0.95) {
  # TODO(R-5): surface this from summary.recmfit() via profile = TRUE.
  th <- object$theta
  obj <- object$objective
  se <- sqrt(pmax(diag(object$vcov), 0))[which]
  if (!is.finite(se) || se <= 0) se <- 0.5
  grid <- seq(th[which] - span * se, th[which] + span * se,
              length.out = ngrid)
  free <- setdiff(seq_along(th), which)

  vals <- vapply(grid, function(g) {
    t2 <- th
    t2[which] <- g
    if (!length(free)) return(obj(t2))
    o <- stats::optim(th[free],
                      function(v) {
                        t3 <- t2
                        t3[free] <- v
                        obj(t3)
                      },
                      method = if (length(free) == 1) "Brent" else
                        "Nelder-Mead",
                      lower = if (length(free) == 1) th[free] - 10 else -Inf,
                      upper = if (length(free) == 1) th[free] + 10 else Inf,
                      control = list(maxit = 800))
    o$value
  }, numeric(1))

  a0s <- vapply(grid, function(g) {
    t2 <- th
    t2[which] <- g
    al <- object$alpha_fn(t2)
    if (any(!is.finite(al))) NA_real_ else .alpha_to_a(al)[1]
  }, numeric(1))

  vmin <- min(vals, na.rm = TRUE)
  cut <- vmin * (1 + stats::qchisq(level, 1) / (object$n - object$npar))
  inside <- grid[!is.na(vals) & vals <= cut]
  structure(list(grid = grid, value = vals, a0 = a0s, cutoff = cut,
                 estimate = th[which], level = level,
                 name = names(object$par)[which],
                 interval = if (length(inside)) range(inside) else
                   c(NA_real_, NA_real_),
                 rel_range = diff(range(vals, na.rm = TRUE)) / vmin),
            class = "recm_profile")
}

#' @rdname recm_profile
#' @export
print.recm_profile <- function(x, ...) {
  cat("\nProfile of the criterion in ", x$name, "\n", sep = "")
  cat("  estimate            : ", signif(x$estimate, 5), "\n", sep = "")
  cat("  ", 100 * x$level, "% profile interval: [",
      signif(x$interval[1], 5), ", ", signif(x$interval[2], 5), "]\n",
      sep = "")
  cat("  relative range of criterion: ", signif(x$rel_range, 3),
      if (x$rel_range < .TOL_FLAT_PROFILE) {
        "   <-- nearly flat, weakly identified"
      } else {
        ""
      }, "\n\n", sep = "")
  invisible(x)
}

#' @rdname recm_profile
#' @export
plot.recm_profile <- function(x, ...) {
  graphics::plot(x$grid, x$value, type = "l", lwd = 2, col = "steelblue",
                 xlab = x$name, ylab = "profiled criterion", ...)
  graphics::abline(h = x$cutoff, lty = 3)
  graphics::abline(v = x$estimate, lty = 2)
  invisible(x)
}


# ---- bootstrap ----

#' Bootstrap standard errors that propagate VAR estimation error
#'
#' The conditional standard errors reported by [summary()] treat the
#' auxiliary VAR as known. This function relaxes that by resampling both the
#' VAR and the equation error and re-running the whole pipeline per
#' replication.
#'
#' @param object an object of class `recmfit`.
#' @param R integer, number of bootstrap replications.
#' @param var_error logical; recursively resample the auxiliary VAR, giving
#'   a fresh `H*` each replication.
#' @param eq_error logical; wild-bootstrap the equation error (Rademacher
#'   weights).
#' @param seed integer or `NULL`; passed to [set.seed()].
#'
#' @return A list of class `recm_boot` with `cov` (the bootstrap covariance
#'   of the parameter vector), `replicates` (an `R` by `npar` matrix, failed
#'   replications removed), `se`, `interval` (percentile), `R_ok` and
#'   `R_failed`.
#'
#' @details
#' Each replication simulates the auxiliary VAR forward from resampled
#' residuals, rebuilds `ystar` by cumulation, regenerates `y` from the
#' fitted decision rule with a wild-bootstrapped error, and re-estimates.
#'
#' **Deviation from `docs/03-api-reference.md`:** that document describes
#' the equation step as a *fixed-design* wild bootstrap. A fixed design is
#' not available here -- the error-correction regressor is built from `y`,
#' so perturbing `y` necessarily moves the design. This implementation is
#' therefore a recursive (dynamic) bootstrap, which is the correct
#' construction for a dynamic equation. Roadmap item R-1.
#'
#' @section What this does not cover:
#' Uncertainty in `ystar` is **not** propagated beyond its VAR
#' representation. `ystar` normally comes from a cointegrating regression
#' estimated outside the package. Superconsistency makes this second-order
#' asymptotically, but at `T` around 170 quarters it is not zero. Roadmap
#' item R-1.
#'
#' In Monte Carlo, the conditional standard errors came in at roughly half
#' the true sampling standard deviation; this is the intended fix, but its
#' own calibration is not yet verified (docs/04-testing.md, "What is not
#' tested yet").
#'
#' @seealso [recm_profile()], [recm_estimate()].
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5,
#'                                 method = "recursive"))
#' dat <- data.frame(ystar = cumsum(dys))
#' dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))
#' fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
#'                      var_lags = 2, restarts = 0, quiet = TRUE)
#' bs <- recm_boot(fit, R = 3, seed = 1)
#' bs$R_ok
#'
#' @export
recm_boot <- function(object, R = 299, var_error = TRUE, eq_error = TRUE,
                      seed = NULL) {
  # TODO(R-1): accept an optional ystar_fit and redraw the cointegrating
  # coefficients per replication. Cache the VAR lag matrices -- do not
  # recompute the design.
  if (!is.null(seed)) set.seed(seed)
  if (!inherits(object, "recmfit")) stop("`object` must be a `recmfit`")

  env <- environment(object$objective)
  data <- env$data
  if (is.null(data)) {
    stop("the fit's closures no longer carry their data; a `recmfit` is ",
         "not portable across sessions (see docs/01-architecture.md)")
  }
  if (object$growth) {
    stop("recm_boot() does not support `growth` yet; the correction ",
         "depends on sum(d_i), which moves with every replication")
  }

  vc <- env$vc
  Xe <- env$Xe
  ok <- env$ok
  yv <- env$yv
  ysv <- env$ysv
  Wm <- env$W
  p <- vc$p
  kv <- vc$k
  Tn <- length(yv)
  a <- object$a
  m <- object$m
  vres <- vc$resid
  eres <- object$residuals

  # Column names of Xe beyond d.ystar, in the order .var_companion saw them.
  exp_names <- colnames(Xe)[-1]
  diff_exp <- isTRUE(env$diff_exp)
  # The user's `expectations` columns, in levels, matched by the names
  # .getcols() produced.
  exp_cols <- if (length(exp_names)) {
    intersect(exp_names, names(data))
  } else {
    character(0)
  }
  if (length(exp_names) && length(exp_cols) != length(exp_names)) {
    stop("recm_boot() cannot regenerate the expectations block: ",
         "columns ", paste(setdiff(exp_names, exp_cols), collapse = ", "),
         " are not plain columns of `data`")
  }

  reps <- matrix(NA_real_, R, object$npar,
                 dimnames = list(NULL, names(object$par)))

  for (b in seq_len(R)) {
    Xb <- Xe

    if (var_error) {
      # Recursive residual bootstrap of the auxiliary VAR.
      draw <- vres[sample.int(nrow(vres), nrow(vres), replace = TRUE), ,
                   drop = FALSE]
      start <- p + 1L
      Xb[seq_len(start), ] <- Xe[seq_len(start), , drop = FALSE]
      for (tt in seq.int(start + 1L, Tn)) {
        lags <- unlist(lapply(seq_len(p), function(j) Xb[tt - j, ]))
        Xb[tt, ] <- drop(c(1, lags) %*% vc$coef) +
          draw[((tt - start - 1L) %% nrow(draw)) + 1L, ]
      }
    }

    # Rebuild ystar (and any expectations levels) by cumulation.
    dys_b <- Xb[, 1]
    ystar_b <- ysv
    ystar_b[-1] <- ysv[1] + cumsum(ifelse(is.na(dys_b[-1]), 0, dys_b[-1]))

    data_b <- data
    data_b[[object$ystar]] <- ystar_b
    if (length(exp_cols) && diff_exp) {
      for (j in seq_along(exp_cols)) {
        v <- Xb[, 1 + j]
        lev <- data[[exp_cols[j]]]
        lev[-1] <- lev[1] + cumsum(ifelse(is.na(v[-1]), 0, v[-1]))
        data_b[[exp_cols[j]]] <- lev
      }
    } else if (length(exp_cols)) {
      for (j in seq_along(exp_cols)) data_b[[exp_cols[j]]] <- Xb[, 1 + j]
    }

    # Forward loading under the new H*, at the FITTED alpha.
    vcb <- tryCatch(.var_companion(Xb, p), error = function(e) NULL)
    if (is.null(vcb)) next
    hb <- .hvec(object$alpha, object$beta, vcb$H, 2L)
    if (is.null(hb)) next
    Slag_b <- rbind(NA, vcb$states[-Tn, , drop = FALSE])
    Zb <- drop(Slag_b %*% hb)
    if (object$free_forward) Zb <- object$a_f * Zb / object$scalars$sum_d

    # Wild-bootstrap the equation error and regenerate y recursively.
    eb <- numeric(Tn)
    if (eq_error) {
      eb[ok] <- eres * sample(c(-1, 1), length(eres), replace = TRUE)
    }
    y_b <- yv
    dy_b <- c(NA, diff(yv))
    for (tt in ok) {
      v <- a[1] * (y_b[tt - 1] - ystar_b[tt - 1])
      if (m > 1) {
        for (j in seq_len(m - 1)) v <- v + a[j + 1] * dy_b[tt - j]
      }
      v <- v + Zb[tt] + eb[tt]
      if (length(object$vars_names)) {
        v <- v + drop(Wm[tt, , drop = FALSE] %*% object$delta)
      }
      dy_b[tt] <- v
      y_b[tt] <- y_b[tt - 1] + v
    }
    data_b[[object$y]] <- y_b

    cl <- object$call
    cl$data <- quote(data_b)
    cl$quiet <- TRUE
    cl$subset <- NULL
    fb <- tryCatch(eval(cl, list(data_b = data_b), environment(recm_boot)),
                   error = function(e) NULL,
                   warning = function(w) NULL)
    if (!is.null(fb) && length(fb$par) == object$npar) reps[b, ] <- fb$par
  }

  good <- stats::complete.cases(reps)
  if (sum(good) < 2) {
    stop("recm_boot(): only ", sum(good), " of ", R,
         " replications produced a usable fit")
  }
  rr <- reps[good, , drop = FALSE]
  structure(list(cov = stats::var(rr),
                 replicates = rr,
                 se = apply(rr, 2, stats::sd),
                 interval = t(apply(rr, 2, stats::quantile,
                                    probs = c(0.025, 0.975))),
                 R_ok = sum(good), R_failed = R - sum(good)),
            class = "recm_boot")
}
