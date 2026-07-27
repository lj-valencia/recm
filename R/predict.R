# The predict method for recmfit objects: in-sample and out-of-sample
# prediction of the fitted decision rule, in differences or in levels, with
# optional confidence and prediction intervals.
#
# Prediction is one step ahead and CONDITIONAL: every right-hand-side
# quantity -- the error-correction term, the lagged differences, the VAR
# states behind the forward sum -- is read from the data supplied, never
# from the model's own earlier predictions. There is no dynamic simulation
# here; recm_boot() is what iterates the rule forward.


# ---- design ----

# Everything the decision rule needs from a data frame, none of which
# depends on the parameter vector. Built once and reused across the Jacobian
# steps and the bootstrap replicates, which is why it is separate from the
# evaluation below.
#
# The specification is read from recm_estimate()'s own frame rather than
# re-derived from `object$call`: the call holds unevaluated expressions, so
# `growth = gcol` with `gcol <- "trend"` would have looked for a column
# literally named `gcol`.
.recm_predict_design <- function(object, data) {
  env <- .fit_env(object, c("Xe", "diff_exp"))

  yv <- data[[object$y]]
  if (is.null(yv)) {
    stop("column '", object$y, "' (decision variable) not found in `data`")
  }
  ysv <- data[[object$ystar]]
  if (is.null(ysv)) {
    stop("column '", object$ystar,
         "' (frictionless target) not found in `data`")
  }
  Tn <- length(yv)
  m <- object$m

  dy  <- c(NA_real_, diff(yv))
  dys <- c(NA_real_, diff(ysv))
  ecm <- .lagv(yv - ysv, 1)
  DYL <- if (m > 1) {
    matrix(vapply(seq_len(m - 1), function(j) .lagv(dy, j), numeric(Tn)),
           nrow = Tn)
  } else {
    NULL
  }

  # The expectations block, assembled exactly as recm_estimate() assembled
  # it: d(ystar) first -- position 2 of the state vector is load-bearing --
  # then the user's variables, differenced or not as they were at
  # estimation.
  exp_names <- colnames(env$Xe)[-1]
  Xe <- cbind(d.ystar = dys)
  if (length(exp_names)) {
    Ex <- .getcols(exp_names, data)
    if (isTRUE(env$diff_exp)) {
      Ex <- rbind(NA_real_, apply(Ex, 2, diff))
    }
    Xe <- cbind(Xe, Ex[, setdiff(colnames(Ex), colnames(Xe)), drop = FALSE])
  }
  if (ncol(Xe) != object$var$k) {
    stop("the expectations block built from `data` has ", ncol(Xe),
         " column", if (ncol(Xe) == 1) "" else "s", " but the fitted VAR has ",
         object$var$k, ": `data` does not match the estimated specification")
  }

  # Through the mechanism, and at the FITTED H. The states are new -- they
  # are built from `data` -- but nothing about the VAR is re-estimated:
  # re-estimating it would predict from a different model than the one that
  # was fitted, which is the whole reason the VAR sits outside the optimiser
  # in the first place.
  zmech <- .zmech_var(object$var, object$beta,
                      states = .var_states(Xe, object$var$p))

  gr <- if (object$growth) {
    g <- data[[object$growth_name]]
    if (is.null(g)) {
      stop("column '", object$growth_name, "' (trend growth) not found in ",
           "`data`; the fit carries a growth correction")
    }
    g
  } else {
    rep(0, Tn)
  }

  nw <- length(object$vars_names)
  list(Tn = Tn, m = m, nw = nw, ylag = .lagv(yv, 1),
       ecm = ecm, DYL = DYL, zmech = zmech, gr = gr,
       W = if (nw) .getcols(object$vars_names, data) else NULL)
}


# ---- the rule at an arbitrary parameter vector ----

# `p` is on the c(theta, a_f?, delta) layout that recm_estimate() names its
# `par` with, and is read POSITIONALLY, the same way the estimator reads
# `r$lin` -- a name-based lookup would return NA rather than fail if the
# names were ever dropped.
#
# Returns all-NA for an inadmissible `p` rather than stopping, because its
# callers are .jacnum(), which steps into the boundary by design, and the
# bootstrap loop, which averages over replicates that may not all be
# admissible. predict.recmfit() checks the point evaluation itself, where a
# refusal is meaningful and is raised.
.recm_predict_eval <- function(object, dg, p, type = c("diff", "level")) {
  type <- match.arg(type)
  bad <- rep(NA_real_, dg$Tn)

  nth <- length(object$theta)
  kk <- .k_from_theta(p[seq_len(nth)], object$m, object$cost)
  if (any(!is.finite(kk)) || any(kk <= 0)) return(bad)
  lq <- .lq_alpha(kk, object$beta)
  if (is.null(lq) || !lq$converged || lq$maxeig >= 1) return(bad)
  al <- lq$alpha
  s <- .scalars(al, object$beta)
  if (is.null(s)) return(bad)
  Zm <- dg$zmech$z(al, s)
  if (is.null(Zm)) return(bad)

  aa <- .alpha_to_a(al)
  corr <- if (object$growth) (1 - sum(aa[-1]) - s$sum_d) * dg$gr else 0
  Zraw <- Zm + corr
  lin <- p[-seq_len(nth)]

  out <- aa[1] * dg$ecm
  if (dg$m > 1) out <- out + drop(dg$DYL %*% aa[-1])
  out <- out + if (object$free_forward) lin[[1]] * Zraw / s$sum_d else Zraw
  if (dg$nw) {
    out <- out + drop(dg$W %*% lin[length(lin) - dg$nw + seq_len(dg$nw)])
  }
  if (type == "level") out <- dg$ylag + out
  unname(out)
}


# Row-wise quantile that tolerates an all-NA row. Those rows are
# observations outside the mechanism's support, where NA is the answer, and
# stats::quantile() would otherwise be handed a zero-length vector.
.rowquant <- function(M, prob) {
  apply(M, 1, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) NA_real_ else unname(stats::quantile(v, prob))
  })
}


#' Predict method for rational error correction models
#'
#' One-step-ahead predictions from a fitted `recmfit`, with optional
#' confidence or prediction intervals, following [stats::predict.lm()]
#' conventions.
#'
#' @param object an object of class `recmfit`.
#' @param newdata optional `data.frame` of new observations, ordered in time
#'   and with no gaps. It must carry the decision variable, the frictionless
#'   target, and any `vars`, `expectations` and `growth` columns the fit was
#'   estimated with, under the same names. `NULL` (default) predicts on the
#'   data the fit was estimated on, after any `subset`.
#' @param interval `"none"` (default), `"confidence"` or `"prediction"`.
#' @param level numeric(1) in (0, 1), the interval's coverage.
#' @param type `"diff"` (default) predicts \eqn{\Delta y_t}; `"level"`
#'   predicts \eqn{y_t = y_{t-1} + \Delta y_t}, taking \eqn{y_{t-1}} from the
#'   data rather than from the previous prediction.
#' @param boot_object optional object of class `recm_boot`. Supplying one
#'   propagates the auxiliary VAR's estimation error into the interval;
#'   without it the interval conditions on the VAR and is too narrow. See
#'   the Standard errors section.
#' @param quiet logical; suppress the conditional-standard-error warning
#'   raised when an interval is requested without a `boot_object`.
#' @param ... ignored, for consistency with the generic.
#'
#' @return With `interval = "none"`, a numeric vector of length
#'   `nrow(newdata)` (or of the fitted data). Otherwise a data frame with
#'   columns `fit`, `se.fit`, `lwr` and `upr`, one row per observation.
#'
#'   Observations the model cannot be evaluated at are `NA`: the leading
#'   rows consumed by the differencing, the error-correction lag, the `m - 1`
#'   lags of \eqn{\Delta y} and the auxiliary VAR's lag structure. Rows of
#'   the output line up with rows of the input throughout, so no realignment
#'   is needed; `object$index` gives the rows that were estimated on.
#'
#' @details
#' Prediction is one step ahead and conditional. Every right-hand-side
#' quantity is read from the data supplied, never from the model's own
#' earlier predictions, so `type = "level"` is \eqn{y_{t-1}} plus the
#' predicted difference and not a simulated path. On the estimation sample
#' and at the fitted parameters this reproduces `fitted(object)` exactly.
#'
#' The auxiliary VAR is **not** re-estimated on `newdata`. The forward sum
#' uses the fitted `H` applied to states built from the new observations,
#' which is what makes an out-of-sample prediction a prediction from the
#' fitted model rather than from a second, differently fitted one.
#'
#' @section Intervals:
#' Both interval types cover the fitted decision rule only. Neither covers
#' error in `ystar`, which is normally itself an estimate -- see roadmap
#' item R-1.
#'
#' Without a `boot_object` the standard error is the delta method applied to
#' the numerical Jacobian of the prediction with respect to the full
#' parameter vector, sandwiched with `vcov(object)`.
#'
#' With one, the parameter uncertainty comes from the bootstrap replicates:
#' `se.fit` is their standard deviation, and a `"confidence"` interval is
#' their percentile interval, so it is not required to be symmetric about
#' `fit` and need not contain it if the bootstrap distribution is biased.
#'
#' A `"prediction"` interval adds the residual variance
#' \eqn{\hat\sigma^2 = \sum e_t^2 / (n - k)} in quadrature and uses a normal
#' critical value on both routes. It is deliberately not formed by simulating
#' residual draws: that would make repeated calls on one fit return different
#' numbers.
#'
#' @section Standard errors:
#' `Z_t` is a generated regressor. It depends on the estimated auxiliary VAR
#' *and* on `ystar`, which usually comes from a first-stage cointegrating
#' regression. In Monte Carlo the conditional standard errors came in at
#' roughly half the true sampling standard deviation for the PAC parameters.
#' Intervals computed without a `boot_object` inherit that: they are too
#' narrow, and the method says so unless `quiet = TRUE`.
#'
#' @seealso [recm_estimate()], [recm_boot()], [stats::predict.lm()].
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
#'
#' # In-sample, in differences: equals fitted() on the estimated rows.
#' pr <- predict(fit)
#' all.equal(pr[fit$index], fitted(fit))
#'
#' # In levels, with a conditional confidence interval.
#' head(predict(fit, type = "level", interval = "confidence", quiet = TRUE))
#'
#' @export
predict.recmfit <- function(object, newdata = NULL,
                            interval = c("none", "confidence", "prediction"),
                            level = 0.95,
                            type = c("diff", "level"),
                            boot_object = NULL,
                            quiet = FALSE,
                            ...) {
  if (!inherits(object, "recmfit")) stop("`object` must be a `recmfit`")
  interval <- match.arg(interval)
  type <- match.arg(type)
  ok_level <- is.numeric(level) && length(level) == 1L &&
    is.finite(level) && level > 0 && level < 1
  if (!ok_level) {
    stop("`level` must be a single number strictly between 0 and 1")
  }

  data <- if (is.null(newdata)) .fit_env(object)$data else newdata
  dg <- .recm_predict_design(object, data)
  fit_vals <- .recm_predict_eval(object, dg, object$par, type)
  if (all(is.na(fit_vals))) {
    stop("no prediction could be formed on this data: either it has too few ",
         "usable rows for the lag structure (", object$m - 1L, " lags of dy ",
         "and ", object$var$p, " VAR lags), or the fitted alpha is ",
         "inadmissible against its state vector")
  }
  if (interval == "none") return(fit_vals)

  # The residual standard error with the degrees-of-freedom correction the
  # NLS covariance uses -- not sd(), which demeans residuals that are not
  # required to have mean zero.
  sigma2 <- sum(object$residuals^2) / (object$n - object$npar)
  tail_p <- (1 - level) / 2
  zc <- stats::qnorm(1 - tail_p)

  if (!is.null(boot_object)) {
    if (!inherits(boot_object, "recm_boot")) {
      stop("`boot_object` must be of class `recm_boot`")
    }
    reps <- boot_object$replicates
    if (ncol(reps) != length(object$par)) {
      stop("`boot_object` carries ", ncol(reps), " parameters but `object` ",
           "has ", length(object$par), "; they are not from the same fit")
    }
    P <- vapply(seq_len(nrow(reps)),
                function(b) .recm_predict_eval(object, dg, reps[b, ], type),
                numeric(dg$Tn))
    se_fit <- apply(P, 1, stats::sd, na.rm = TRUE)
    if (interval == "confidence") {
      lwr <- .rowquant(P, tail_p)
      upr <- .rowquant(P, 1 - tail_p)
    } else {
      se_pred <- sqrt(se_fit^2 + sigma2)
      lwr <- fit_vals - zc * se_pred
      upr <- fit_vals + zc * se_pred
    }
  } else {
    # Never remove this. It is always true, and it is the difference between
    # an interval and a correct one -- see the Standard errors section.
    if (!quiet) {
      warning("standard errors and intervals condition on the auxiliary VAR ",
              "and on `ystar`; in Monte Carlo the conditional standard ",
              "errors came in at roughly half the true sampling standard ",
              "deviation. Pass a `recm_boot` object via `boot_object`.",
              call. = FALSE)
    }
    V <- object$vcov
    if (all(is.na(V))) {
      stop("the fit has no usable covariance matrix, so no interval can be ",
           "formed; pass a `recm_boot` object via `boot_object`")
    }
    J <- .jacnum(function(pp) .recm_predict_eval(object, dg, pp, type),
                 object$par)
    # rowSums((J V) * J) is diag(J V J') without forming the T by T matrix.
    v_fit <- pmax(rowSums((J %*% V) * J), 0)
    se_fit <- sqrt(v_fit)
    se_use <- if (interval == "confidence") se_fit else sqrt(v_fit + sigma2)
    lwr <- fit_vals - zc * se_use
    upr <- fit_vals + zc * se_use
  }

  data.frame(fit = fit_vals, se.fit = se_fit, lwr = lwr, upr = upr)
}
