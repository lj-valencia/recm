# Presentation only, like methods.R. Implements the dual-format reporting
# described in docs/02-math-spec.md section 7.
#
# Normalisation note: f_i = d_i / sum_j d_j uses the CLOSED-FORM total
# sum_d, not the sum of a truncated vector. That keeps a_f * f_i = d_i
# exact under the Euler restriction; the price is that summing a truncated
# f over a finite horizon reaches 1 only up to the tail, which decays like
# rho(G)^i.

# Horizon used for the internal d/f vectors quoted by equation() and
# summary(). Long enough that the truncated tail is negligible at any
# rho(G) the estimator will accept.
.HORIZON_REPORT <- 200L


# ---- the fitted equation ----

#' Write out the fitted equation
#'
#' Prints the estimated decision rule with coefficients substituted, in the
#' explicit format, the compressed (LENS / FRB-US) format, or both.
#'
#' @param object an object of class `recmfit`.
#' @param digits integer, significant digits used in printing.
#' @param format one of `"both"`, `"explicit"` or `"compressed"`. In the
#'   explicit format the coefficient on the forward sum is 1 by
#'   construction; in the compressed format it is `a_f` and the lead weights
#'   are normalised to `f_i` summing to 1.
#' @param max_lags integer, truncates the printed lag list. The count of
#'   suppressed lags is shown.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Invisibly, a list with `a`, `a_f`, `d` and `f`.
#' @seealso `docs/02-math-spec.md` section 7; [lead_weights()].
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
#' equation(fit, format = "compressed")
#'
#' @export
equation <- function(object, ...) {
  UseMethod("equation")
}

#' @rdname equation
#' @export
equation.recmfit <- function(object, digits = 4,
                             format = c("both", "explicit", "compressed"),
                             max_lags = 12, ...) {
  format <- match.arg(format)
  a <- object$a
  d <- .dweights(object$alpha, object$beta, .HORIZON_REPORT)$d
  f <- d / object$scalars$sum_d
  nm <- object$y
  ns <- object$ystar
  fmt <- function(v) sprintf("%+.*f", digits, v)

  if (format %in% c("both", "explicit")) {
    cat("\nEXPLICIT\n")
    cat("d.", nm, "_t = ", sprintf("%.*f", digits, a[1]),
        "*(", nm, " - ", ns, ")_{t-1}\n", sep = "")
    if (object$m > 1) {
      kk <- min(max_lags, object$m - 1)
      for (i in seq_len(kk)) {
        cat("          ", fmt(a[i + 1]), "*d.", nm, "_{t-", i, "}\n",
            sep = "")
      }
      if (object$m - 1 > kk) {
        cat("          ... (", object$m - 1 - kk,
            " further lags suppressed)\n", sep = "")
      }
    }
    if (length(object$vars_names)) {
      for (j in seq_along(object$vars_names)) {
        cat("          ", fmt(object$delta[j]), "*", object$vars_names[j],
            "_t\n", sep = "")
      }
    }
    cat("          + Z_t + e_t\n")
    cat("   Z_t = sum_i d_i E_{t-1}[d.", ns, "_{t+i}],  d = (",
        paste(sprintf("%.*f", digits, d[1:6]), collapse = ", "), ", ...)\n",
        sep = "")
  }

  if (format %in% c("both", "compressed")) {
    cat("\nCOMPRESSED\n")
    cat("d.", nm, "_t = ", sprintf("%.*f", digits, a[1]),
        "*(", nm, " - ", ns, ")_{t-1}", sep = "")
    if (object$m > 1) {
      cat(" + sum_{i=1}^", object$m - 1, " a_i*d.", nm, "_{t-i}", sep = "")
    }
    cat("\n          ", sprintf("%+.*f", digits, object$a_f),
        " * sum_{i>=0} f_i E_{t-1}[d.", ns, "_{t+i}]", sep = "")
    if (length(object$vars_names)) cat(" + delta'W_t")
    cat(" + e_t\n")
    cat("   sum a_i = ",
        sprintf("%.*f", digits, if (object$m > 1) sum(a[-1]) else 0),
        ",  a_f = ", sprintf("%.*f", digits, object$a_f), "\n", sep = "")
    cat("   f = (", paste(sprintf("%.*f", digits, f[1:6]), collapse = ", "),
        ", ...),  sum f_i = 1\n", sep = "")
  }

  cat("   sigma_e = ", signif(stats::sd(object$residuals), digits),
      ",  n = ", object$n, "\n\n", sep = "")
  invisible(list(a = a, a_f = object$a_f, d = d, f = f))
}


# ---- lead weights ----

#' Lead weights on expected target growth
#'
#' Returns the weights \eqn{d_i} on \eqn{E_{t-1}[\Delta y^*_{t+i}]}, or the
#' normalised \eqn{f_i = d_i / \sum_j d_j}, with delta-method standard
#' errors. Nothing is re-estimated: \eqn{d_i(\theta)} is differentiated
#' numerically and sandwiched with the stored covariance.
#'
#' @param object an object of class `recmfit`.
#' @param horizon integer, highest lead index returned.
#' @param level numeric(1) confidence level for `lo` and `hi`.
#' @param normalised logical; if `TRUE` return `f_i` instead of `d_i`. These
#'   are divided by the closed-form total \eqn{\sum_j d_j}, so a truncated
#'   vector sums to 1 only up to a tail decaying like \eqn{\rho(G)^i}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A data frame with columns `i`, `weight`, `se`, `lo` and `hi`. The
#'   Jacobian used for the standard errors is attached as attribute
#'   `"jacobian"`.
#' @seealso `docs/02-math-spec.md` section 5; [equation()].
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
#' head(lead_weights(fit, horizon = 6))
#'
#' @export
lead_weights <- function(object, ...) {
  UseMethod("lead_weights")
}

#' @rdname lead_weights
#' @export
lead_weights.recmfit <- function(object, horizon = 24, level = 0.95,
                                 normalised = FALSE, ...) {
  f <- function(th) {
    al <- object$alpha_fn(th)
    if (any(!is.finite(al))) return(rep(NA_real_, horizon + 1))
    w <- .dweights(al, object$beta, horizon)
    if (is.null(w)) return(rep(NA_real_, horizon + 1))
    if (normalised) w$d / .scalars(al, object$beta)$sum_d else w$d
  }
  th <- object$theta
  est <- f(th)
  J <- .jacnum(f, th)
  Vth <- object$vcov[seq_along(th), seq_along(th), drop = FALSE]
  se <- sqrt(pmax(diag(J %*% Vth %*% t(J)), 0))
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(i = 0:horizon, weight = est, se = se,
                    lo = est - z * se, hi = est + z * se)
  # TODO(R-3): plot.recmfit() should accumulate this before sandwiching
  # rather than summing pointwise standard errors.
  attr(out, "jacobian") <- J
  out
}
