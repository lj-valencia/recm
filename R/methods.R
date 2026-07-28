## Extractors and display for fitted rational error correction models.

#' Methods for fitted rational error correction models
#'
#' @param object,x An object of class `"recm"` returned by [recm()].
#' @param digits Number of significant digits used when printing.
#' @param ... Ignored, present for consistency with the generics.
#'
#' @name recm-methods
NULL

#' @rdname recm-methods
#' @export
coef.recm <- function(object, ...) object$coefficients

#' @rdname recm-methods
#' @export
vcov.recm <- function(object, ...) object$vcov

#' @rdname recm-methods
#' @export
residuals.recm <- function(object, ...) object$residuals

#' @rdname recm-methods
#' @export
fitted.recm <- function(object, ...) object$fitted.values

#' @rdname recm-methods
#' @export
nobs.recm <- function(object, ...) object$nobs

#' @rdname recm-methods
#' @export
logLik.recm <- function(object, ...) {
  n <- object$nobs
  val <- -0.5 * n * (log(2 * pi) + log(object$ssr / n) + 1)
  structure(val, df = length(object$coefficients) + 1L, nobs = n,
            class = "logLik")
}

# The estimated equation, written out with the fitted coefficients.
recm_equation <- function(object, digits = 4L) {
  cf <- object$coefficients
  v <- object$variables
  parts <- sprintf(
    "%s * (%s[t-1] - %s[t-1])",
    format(cf[["ec"]], digits = digits), v$y_star, v$y
  )
  if (object$m > 1L) {
    for (i in seq_len(object$m - 1L)) {
      nm <- paste0("dy_lag", i)
      parts <- c(parts, sprintf(
        "%s * d%s[t-%d]", format(cf[[nm]], digits = digits), v$y, i
      ))
    }
  }
  parts <- c(parts, "Z[t]")
  for (nm in v$w) {
    parts <- c(parts, sprintf(
      "%s * %s[t]", format(cf[[nm]], digits = digits), nm
    ))
  }
  paste0("d", v$y, "[t] = ", paste(parts, collapse = " + "))
}

#' @rdname recm-methods
#' @export
print.recm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Rational error correction model\n\n")
  cat("Call:\n  ", paste(deparse(x$call), collapse = "\n  "), "\n\n", sep = "")
  cat(recm_equation(x, digits), "\n\n")
  cat(sprintf(
    "  m = %d (%d lag%s of the dependent variable)   discount = %s\n",
    x$m, x$m - 1L, if (x$m == 2L) "" else "s", format(x$discount)
  ))
  cat(sprintf(
    "  expectations: %s   estimator: iterative OLS (%d iterations)\n",
    x$expectations, x$iterations
  ))
  cat(sprintf(
    "  observations: %d   residual std. error: %s\n",
    x$nobs, format(sqrt(x$sigma2), digits = digits)
  ))
  invisible(x)
}

#' @rdname recm-methods
#' @export
summary.recm <- function(object, ...) {
  structure(list(fit = object), class = "summary.recm")
}

#' @rdname recm-methods
#' @export
print.summary.recm <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...) {
  obj <- x$fit
  print.recm(obj, digits)

  cat("\nCoefficients:\n")
  tab <- cbind(
    Estimate = obj$coefficients,
    `Std. Error` = obj$std.error,
    `t value` = obj$statistic,
    `Pr(>|t|)` = obj$p.value
  )
  stats::printCoefmat(tab, digits = digits, signif.stars = FALSE)
  cat(
    "Standard errors: just-identified GMM, HAC lag ", obj$hac_lag,
    ", including dZ/da'.\n",
    "They condition on the auxiliary autoregression. In simulation that is\n",
    "not a small matter: with the target process resampled too, they came in\n",
    "at roughly a third of the true sampling standard deviation. Bootstrap\n",
    "the auxiliary model before believing an interval.\n",
    sep = ""
  )

  cat("\nImplied dynamics:\n")
  cat(sprintf("  A(L) roots (modulus)  %s\n",
              paste(format(Mod(obj$roots), digits = digits),
                    collapse = ", ")))
  cat(sprintf("  total forward loading sum(d_i)  %s\n",
              format(obj$d_sum, digits = digits)))
  cat(sprintf("  mean lead  %s quarters   half-life  %s periods\n",
              format(obj$mean_lead, digits = digits),
              format(obj$half_life, digits = digits)))
  cat(sprintf("  adjustment costs k  %s   (k0 = A(1)A(beta), residual %.1e)\n",
              paste(format(obj$structural$k, digits = digits),
                    collapse = ", "),
              obj$structural$residual))
  if (!obj$structural$positive) {
    cat("  NOTE: some k_j <= 0, so the estimated reduced form does not come\n",
        "        from a convex polynomial adjustment cost problem.\n", sep = "")
  }

  cat("\nGrowth neutrality:\n")
  if (obj$growth$restricted) {
    cat(sprintf(
      "  R(a) = 1 - sum(a_i) - sum(d_i) imposed; residual gap %.2e\n",
      obj$growth$gap
    ))
  } else {
    cat(sprintf(
      "  slack at discount = %s: R(a) = %.2e for any coefficients\n",
      format(obj$discount), obj$growth$gap
    ))
  }
  if (!is.na(obj$growth$implied_level_gap)) {
    cat(sprintf(
      "  target drift %s, implied permanent level gap %s\n",
      format(obj$growth$drift, digits = digits),
      format(obj$growth$implied_level_gap, digits = digits)
    ))
  }

  fl <- obj$forward_loading
  cat("\nForward loading, restricted against free (levels 2 vs 3):\n")
  cat(sprintf(
    "  restricted %s (s.e. %s)   free %s (s.e. %s)\n",
    format(fl$restricted, digits = digits),
    format(fl$restricted_se, digits = digits),
    format(fl$free, digits = digits),
    format(fl$free_se, digits = digits)
  ))
  cat(
    "  The two standard errors omit different terms, so no test statistic\n",
    "  is formed; compare the intervals by eye.\n", sep = ""
  )

  cat("\nAuxiliary model: AR(", obj$aux$order, ") in d", obj$variables$y_star,
      " with a constant (order by AIC)\n", sep = "")
  cat(sprintf("  constant %s   ar %s   rho(Phi) %s\n",
              format(obj$aux$const, digits = digits),
              paste(format(obj$aux$ar, digits = digits), collapse = ", "),
              format(obj$aux$rho, digits = digits)))

  if (!obj$converged) {
    cat("\nWARNING: the zig-zag did not converge.\n")
  }
  invisible(x)
}
