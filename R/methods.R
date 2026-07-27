# Presentation only. These methods must NOT recompute anything that changes
# the answer. They may recompute derived quantities -- delta-method
# Jacobians, lead weights at new horizons -- from stored inputs.
# See docs/01-architecture.md, "Module boundaries".

# Printed unconditionally by print.summary.recmfit(). Do NOT add a switch
# to suppress it: it is always true, and it is listed under "explicitly out
# of scope" in docs/05-roadmap.md.
.CAVEAT_GENERATED_REGRESSOR <- paste0(
  "Note: Z_t is a generated regressor -- it depends on the auxiliary ",
  "VAR\nand on %s, both estimated. The standard errors above are ",
  "conditional\nand optimistic; in Monte Carlo they came in at roughly ",
  "half the true\nsampling standard deviation. Use recm_boot(), and check ",
  "recm_profile()\nfor curvature.\n"
)


# ---- basic accessors ----

#' Methods for fitted rational error correction models
#'
#' @param object,x an object of class `recmfit`, or of class
#'   `summary.recmfit` for the corresponding print method.
#' @param digits integer, significant digits used in printing.
#' @param parm,level passed to [confint()]; `parm` selects parameters and
#'   `level` is the confidence level.
#' @param horizon integer, highest lead index shown by `plot()`.
#' @param ... further arguments passed to or from other methods.
#'
#' @return `coef()` returns the named parameter vector, `vcov()` its
#'   covariance, `residuals()` and `fitted()` numeric vectors over the
#'   retained sample `object$index`, `confint()` a two-column matrix, and
#'   `logLik()` an object of class `logLik` computed under Gaussian errors.
#'   `print()` returns its argument invisibly. `summary()` returns an object
#'   of class `summary.recmfit`. `plot()` is called for its side effect and
#'   returns the lead-weight data frame invisibly.
#'
#' @details
#' `summary()` prints, in one call:
#'
#' 1. estimated parameters on the optimisation scale, `printCoefmat` style
#' 2. cost parameters on the natural scale, delta method
#' 3. **Format 1, explicit**: `a_0`, `a_1 ... a_{m-1}`, `delta`, with the
#'    forward coefficient fixed at 1
#' 4. **Format 2, compressed**: `a_0`, `sum a_i`, `a_f`, mean lead,
#'    half-life, and the `f_i` table
#' 5. fit and diagnostics: residual standard error, R-squared, Ljung-Box,
#'    Durbin-Watson, `rho(G)*rho(H)`, cost admissibility, Hansen J,
#'    first-stage F
#'
#' The generated-regressor caveat is printed unconditionally.
#'
#' @seealso [recm_estimate()], [equation()], [lead_weights()],
#'   [recm_profile()].
#' @name recmfit-methods
NULL

#' @rdname recmfit-methods
#' @export
coef.recmfit <- function(object, ...) object$par

#' @rdname recmfit-methods
#' @export
vcov.recmfit <- function(object, ...) object$vcov

#' @rdname recmfit-methods
#' @export
residuals.recmfit <- function(object, ...) object$residuals

#' @rdname recmfit-methods
#' @export
fitted.recmfit <- function(object, ...) object$fitted

#' @rdname recmfit-methods
#' @export
confint.recmfit <- function(object, parm, level = 0.95, ...) {
  cf <- object$par
  se <- sqrt(pmax(diag(object$vcov), 0))
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(cf - z * se, cf + z * se)
  colnames(out) <- paste(format(100 * c((1 - level) / 2,
                                        1 - (1 - level) / 2),
                                trim = TRUE, digits = 3), "%")
  rownames(out) <- names(cf)
  if (!missing(parm)) out <- out[parm, , drop = FALSE]
  out
}

#' @rdname recmfit-methods
#' @export
logLik.recmfit <- function(object, ...) {
  n <- object$n
  ssr <- sum(object$residuals^2)
  structure(-0.5 * n * (log(2 * pi) + log(ssr / n) + 1),
            df = object$npar + 1, nobs = n, class = "logLik")
}


# ---- print and summary ----

#' @rdname recmfit-methods
#' @export
print.recmfit <- function(x, ...) {
  cat("\nREC equation: ", x$y, " ~ ", x$ystar,
      if (length(x$vars_names)) {
        paste(" +", paste(x$vars_names, collapse = " + "))
      }, "\n", sep = "")
  cat("m = ", x$m, ", beta = ", x$beta, ", ", toupper(x$method), "\n",
      sep = "")
  cat("a0 = ", signif(x$a[1], 4), "   sum a_i = ",
      signif(if (x$m > 1) sum(x$a[-1]) else 0, 4),
      "   a_f = ", signif(x$a_f, 4), "\n", sep = "")
  cat("residual sd = ", signif(stats::sd(x$residuals), 4), " on ", x$n,
      " observations\n\n", sep = "")
  invisible(x)
}

#' @rdname recmfit-methods
#' @export
summary.recmfit <- function(object, digits = 4, ...) {
  th <- object$theta
  V <- object$vcov
  e <- object$residuals
  nn <- object$n
  df <- nn - object$npar
  Vth <- V[seq_along(th), seq_along(th), drop = FALSE]

  se <- sqrt(pmax(diag(V), 0))
  tv <- object$par / se
  ctab <- cbind(Estimate = object$par, `Std. Error` = se, `t value` = tv,
                `Pr(>|t|)` = 2 * stats::pt(-abs(tv), df))

  nat_fn <- function(t2) .theta_natural(t2, object$m, object$cost)
  ne <- nat_fn(th)
  Jn <- .jacnum(nat_fn, th)
  ns <- sqrt(pmax(diag(Jn %*% Vth %*% t(Jn)), 0))
  ntab <- cbind(Estimate = ne, `Std. Error` = ns,
                `2.5%` = ne - 1.96 * ns, `97.5%` = ne + 1.96 * ns)

  # ---- explicit format: a0, a_1..a_{m-1} ----
  a_fn <- function(t2) {
    al <- object$alpha_fn(t2)
    if (any(!is.finite(al))) rep(NA_real_, object$m) else .alpha_to_a(al)
  }
  ae <- a_fn(th)
  Ja <- .jacnum(a_fn, th)
  as_ <- sqrt(pmax(diag(Ja %*% Vth %*% t(Ja)), 0))
  atab <- cbind(Estimate = ae, `Std. Error` = as_,
                `2.5%` = ae - 1.96 * as_, `97.5%` = ae + 1.96 * as_)
  lag_names <- if (object$m > 1) {
    paste0("a", seq_len(object$m - 1),
           " (dy.l", seq_len(object$m - 1), ")")
  } else {
    character(0)
  }
  rownames(atab) <- c("a0 (ecm)", lag_names)
  if (length(object$vars_names)) {
    iv <- object$npar - length(object$vars_names) +
      seq_along(object$vars_names)
    atab <- rbind(atab,
                  cbind(Estimate = object$par[iv], `Std. Error` = se[iv],
                        `2.5%` = object$par[iv] - 1.96 * se[iv],
                        `97.5%` = object$par[iv] + 1.96 * se[iv]))
  }

  # ---- compressed format: a0, sum a_i, a_f, f_i ----
  comp_fn <- function(t2) {
    al <- object$alpha_fn(t2)
    if (any(!is.finite(al))) return(rep(NA_real_, 5))
    aa <- .alpha_to_a(al)
    s <- .scalars(al, object$beta)
    w <- .dweights(al, object$beta, .HORIZON_REPORT)
    ff <- w$d / s$sum_d
    c(`a0` = aa[1],
      `sum a_i` = if (length(aa) > 1) sum(aa[-1]) else 0,
      `a_f` = s$sum_d,
      `mean lead` = sum((0:.HORIZON_REPORT) * ff),
      `half-life` = .halflife(aa)$halflife)
  }
  ce <- comp_fn(th)
  Jc <- .jacnum(comp_fn, th)
  cs_ <- sqrt(pmax(diag(Jc %*% Vth %*% t(Jc)), 0))
  if (object$free_forward) {
    ii <- length(th) + 1L
    ce["a_f"] <- object$par[ii]
    cs_[3] <- se[ii]
  }
  ctab2 <- cbind(Estimate = ce, `Std. Error` = cs_,
                 `2.5%` = ce - 1.96 * cs_, `97.5%` = ce + 1.96 * cs_)

  lw <- lead_weights(object, 12, normalised = TRUE)

  ssr <- sum(e^2)
  tss <- sum((object$dep - mean(object$dep))^2)
  lb <- stats::Box.test(e, lag = 8, type = "Ljung-Box")

  ans <- list(call = object$call, method = object$method,
              cost = object$cost, m = object$m, beta = object$beta,
              n = nn, df = df, free_forward = object$free_forward,
              growth = object$growth, y = object$y, ystar = object$ystar,
              resid_q = stats::quantile(e, c(0, .25, .5, .75, 1)),
              coefficients = ctab, natural = ntab,
              explicit = atab, compressed = ctab2, f = lw,
              sigma = sqrt(ssr / df), r2 = 1 - ssr / tss,
              lb = lb, dw = sum(diff(e)^2) / ssr,
              rhoG = object$scalars$rhoG, rhoH = object$var$rho_dyn,
              k = object$k, extras = object$extras,
              converged = object$converged, digits = digits,
              vars_names = object$vars_names)
  class(ans) <- "summary.recmfit"
  ans
}

#' @rdname recmfit-methods
#' @export
print.summary.recmfit <- function(x, digits = x$digits, ...) {
  dg <- digits
  cat("\nCall:\n")
  print(x$call)
  cat("\nRational error correction equation, polynomial adjustment costs\n")
  cat("  dependent: ", x$y, "     target: ", x$ystar, "\n", sep = "")
  cat("  m = ", x$m, "   beta = ", x$beta,
      " (calibrated)   costs: ", x$cost, "\n", sep = "")
  cat("  estimator: ", toupper(x$method),
      if (x$free_forward) {
        "   forward loading FREE"
      } else {
        "   forward loading RESTRICTED to sum(d_i)"
      },
      if (x$growth) "   + growth correction" else "", "\n", sep = "")
  if (!x$converged) {
    cat("  ** optimiser did not report convergence **\n")
  }

  cat("\nResiduals:\n")
  rq <- x$resid_q
  names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(signif(rq, dg))

  cat("\n--- Estimated parameters (optimisation scale) ---\n")
  stats::printCoefmat(x$coefficients, digits = dg, signif.stars = TRUE)
  cat("\n--- Cost parameters (natural scale, delta method) ---\n")
  print(signif(x$natural, dg))

  cat("\n================ FORMAT 1: EXPLICIT ================\n")
  cat("d.", x$y, "_t = a0*(", x$y, "-", x$ystar, ")_{t-1}", sep = "")
  if (x$m > 1) {
    cat(" + sum_{i=1}^", x$m - 1, " a_i*d.", x$y, "_{t-i}", sep = "")
  }
  cat(" + Z_t")
  if (length(x$vars_names)) cat(" + delta'W_t")
  cat(" + e_t\n")
  cat("Z_t = sum_i d_i E_{t-1}[d.", x$ystar,
      "_{t+i}]   (coefficient 1)\n\n", sep = "")
  print(signif(x$explicit, dg))

  cat("\n=============== FORMAT 2: COMPRESSED ===============\n")
  cat("d.", x$y, "_t = a0*(", x$y, "-", x$ystar, ")_{t-1}", sep = "")
  if (x$m > 1) cat(" + sum a_i*d.", x$y, "_{t-i}", sep = "")
  cat(" + a_f * sum_i f_i E_{t-1}[d.", x$ystar, "_{t+i}]", sep = "")
  if (length(x$vars_names)) cat(" + delta'W_t")
  cat(" + e_t,   sum f_i = 1\n\n")
  print(signif(x$compressed, dg))
  cat("\nnormalised lead weights f_i:\n")
  fm <- matrix(round(x$f$weight, 4), nrow = 1,
               dimnames = list("f_i", paste0("i=", x$f$i)))
  print(fm)

  cat("\n--- Fit and diagnostics ---\n")
  cat("Residual standard error: ", signif(x$sigma, dg), " on ", x$df,
      " df\n", sep = "")
  cat("Multiple R-squared: ", signif(x$r2, dg), "\n", sep = "")
  cat("Ljung-Box(8): ", signif(x$lb$statistic, dg), " (p = ",
      signif(x$lb$p.value, 3), ")   Durbin-Watson: ", signif(x$dw, dg),
      "\n", sep = "")
  cat("Stability: rho(G) = ", signif(x$rhoG, dg), ", rho(H) = ",
      signif(x$rhoH, dg), ", product = ", signif(x$rhoG * x$rhoH, dg),
      if (x$rhoG * x$rhoH >= 1) "  ** >= 1 **" else "  (< 1, ok)", "\n",
      sep = "")
  cat("Cost admissibility (all k_j > 0): ", all(x$k > 0), "\n", sep = "")
  ex <- x$extras
  if (length(ex)) {
    if (!is.null(ex$J)) {
      cat("Hansen J: ", signif(ex$J, dg), " on ", ex$Jdf, " df (p = ",
          signif(ex$hansen_p, 3), ")\n", sep = "")
    }
    if (!is.null(ex$first_stage_F)) {
      cat("First-stage F: ", signif(ex$first_stage_F, dg),
          if (ex$first_stage_F < 10) "   ** weak instruments **" else "",
          "\n", sep = "")
    }
    if (!is.null(ex$iv_lag)) {
      cat("Instruments: ", ex$n_instruments, " used of ",
          ex$n_instruments_supplied, ", automatic set dated t-",
          ex$iv_lag, "\n", sep = "")
    }
  }
  cat("\n")
  cat(sprintf(.CAVEAT_GENERATED_REGRESSOR, x$ystar))
  cat("\n")
  invisible(x)
}


# ---- plot ----

#' @rdname recmfit-methods
#' @export
plot.recmfit <- function(x, horizon = 24, level = 0.95, ...) {
  lw <- lead_weights(x, horizon, level)
  op <- graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.8, 1.2))
  on.exit(graphics::par(op))

  graphics::plot(lw$i, lw$weight, type = "n",
                 ylim = range(c(lw$lo, lw$hi, 0)),
                 xlab = "lead i", ylab = expression(d[i]),
                 main = "Lead weights")
  graphics::polygon(c(lw$i, rev(lw$i)), c(lw$lo, rev(lw$hi)),
                    col = grDevices::adjustcolor("steelblue", .2),
                    border = NA)
  graphics::lines(lw$i, lw$weight, lwd = 2.4, col = "steelblue")
  graphics::abline(h = 0, col = "grey60")

  # TODO(R-3): this band sums pointwise standard errors, which is
  # conservative in the same way that adding standard deviations instead of
  # variances is wrong. Accumulate attr(lw, "jacobian") first.
  cs <- cumsum(lw$weight)
  graphics::plot(lw$i, cs, type = "l", lwd = 2.4, col = "darkgreen",
                 xlab = "lead i", ylab = "cumulative",
                 main = "Cumulative lead weight")
  graphics::abline(h = x$a_f, lty = 3)

  if (x$m > 1) {
    graphics::plot(seq_len(x$m - 1), x$a[-1], type = "h", lwd = 3,
                   col = "grey30", xlab = "lag i",
                   ylab = expression(a[i]), main = "Lag distribution")
    graphics::abline(h = 0, col = "grey60")
  } else {
    graphics::plot.new()
    graphics::title("no lags (m = 1)")
  }

  hl <- .halflife(x$a, 60)
  graphics::plot(0:60, hl$path, type = "l", lwd = 2.4, col = "firebrick",
                 xlab = "quarters", ylab = "gap",
                 main = "Gap response to a unit shock")
  graphics::abline(h = c(0, .5), lty = c(1, 3),
                   col = c("grey60", "grey40"))
  invisible(lw)
}
