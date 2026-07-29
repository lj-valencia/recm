## Forecasting from a fitted rational error correction model.
##
## The forecast is a dynamic simulation of the decision rule, not a one step
## ahead conditional prediction. y is iterated forward from the end of the
## estimation sample, so the error correction term at t + h is formed from the
## model's own y[t + h - 1] and never from a realised one. `newdata` therefore
## supplies the future path of the target and of the exogenous regressors only.
##
## The equation is linear in its regressors and the coefficient on the forward
## term is fixed at 1, so the decomposition
##
##   dy[t] = a0 * ec[t] + sum_i a_i dy[t-i] + Z[t] + delta' W[t]
##
## is exact term by term, with nothing left over. Cumulating it decomposes the
## level as well, since y[T+h] = y[T] + sum_{s <= h} dy[T+s].
##
## What the decomposition is not is a set of partial derivatives. The error
## correction term is formed from the simulated y, so its bar carries the
## feedback from every term that came before it: holding one column fixed and
## re-running the recursion would not move the total by that column's
## contribution.

# Continue a time index past the end of the sample. A numeric index, which is
# what a ts supplies, and a Date index are extended by their own last step;
# anything else falls back to the row number.
recm_extend_index <- function(idx, h, n_hist) {
  fallback <- n_hist + seq_len(h)
  n <- length(idx)
  if (n < 2L || !(is.numeric(idx) || inherits(idx, "Date"))) {
    return(fallback)
  }
  last <- idx[n]
  step <- as.numeric(unclass(last) - unclass(idx[n - 1L]))
  if (!is.finite(step) || step <= 0) {
    return(fallback)
  }
  last + step * seq_len(h)
}

# The future target projected from the auxiliary autoregression, for the case
# where no `newdata` is given. The autoregression is the one estimated by the
# fit; nothing about it is re-estimated here.
recm_future_projected <- function(object, n_ahead) {
  hist <- object$model$data
  n_hist <- nrow(hist)
  ystar <- hist[, object$variables$y_star]
  dfut <- aux_forecast_path(c(NA_real_, diff(ystar)), object$aux, n_ahead)
  list(
    h = n_ahead,
    ystar = ystar[n_hist] + cumsum(dfut),
    w = NULL,
    index = recm_extend_index(object$model$index, n_ahead, n_hist),
    projected = TRUE
  )
}

# The future target and exogenous regressors read from `newdata`. Any column
# holding the decision variable is ignored: the whole point of the forward
# recursion is that y is produced rather than consumed.
recm_future_supplied <- function(object, newdata, n_ahead) {
  ys_name <- object$variables$y_star
  w_names <- object$variables$w
  nd <- as_model_data(newdata)
  cols <- colnames(nd$x)

  if (!ys_name %in% cols) {
    stop(
      "`newdata` has no column `", ys_name, "`. A forecast iterates the ",
      "decision rule forward, so the future path of the target has to be ",
      "supplied; only `", object$variables$y, "` is produced rather than ",
      "read. Columns found: ", paste(cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  absent <- setdiff(w_names, cols)
  if (length(absent)) {
    stop(
      "`newdata` is missing the exogenous regressor",
      if (length(absent) == 1L) " " else "s ",
      paste(absent, collapse = ", "),
      ", which the fit was estimated with. Supply a future path for ",
      if (length(absent) == 1L) "it" else "each of them", ".",
      call. = FALSE
    )
  }

  h <- nrow(nd$x)
  if (!is.null(n_ahead)) {
    if (n_ahead > h) {
      stop(
        "`n_ahead` is ", n_ahead, " but `newdata` has only ", h, " row",
        if (h == 1L) "" else "s", ".",
        call. = FALSE
      )
    }
    h <- n_ahead
  }
  rows <- seq_len(h)

  ystar <- nd$x[rows, ys_name]
  if (anyNA(ystar)) {
    stop(
      "`", ys_name, "` is missing in ", sum(is.na(ystar)), " of the ", h,
      " forecast rows of `newdata`; the path must be complete.",
      call. = FALSE
    )
  }
  w <- if (length(w_names)) nd$x[rows, w_names, drop = FALSE] else NULL
  if (!is.null(w) && anyNA(w)) {
    stop(
      "the exogenous regressors are missing in ",
      sum(!stats::complete.cases(w)), " of the ", h,
      " forecast rows of `newdata`; the paths must be complete.",
      call. = FALSE
    )
  }

  list(
    h = h, ystar = ystar, w = w,
    # A data frame with no label column and a plain matrix both index as
    # 1..nrow, which would restart the clock; only a genuine index is kept.
    index = if (is.null(nd$index_name)) {
      recm_extend_index(object$model$index, h, nrow(object$model$data))
    } else {
      nd$index[rows]
    },
    projected = FALSE
  )
}

# The exogenous regressors as they enter the equation. Under `tr_exog` the
# first new difference is taken against the last level of the estimation
# sample, which is the reason the fit carries it.
recm_future_w <- function(object, fut) {
  w_names <- object$variables$w
  if (!length(w_names)) {
    return(matrix(0, fut$h, 0L))
  }
  if (!object$tr_exog) {
    return(unname(as.matrix(fut$w)))
  }
  hist <- object$model$data
  last <- hist[nrow(hist), w_names, drop = FALSE]
  if (anyNA(last)) {
    stop(
      "the last observation of the exogenous regressors is missing, so the ",
      "first new difference cannot be formed. Either trim `data` back to the ",
      "last complete row before fitting, or difference the regressors ",
      "yourself and fit with `tr_exog = FALSE`.",
      call. = FALSE
    )
  }
  unname(diff(rbind(last, as.matrix(fut$w))))
}

# The forward term over the forecast rows, at the fitted coefficients and the
# fitted auxiliary autoregression. Under "var" it is h'z[t-1] on states built
# from the new path of the target; under "mce" it is the backward recursion on
# that path, extended past the horizon by autoregressive forecasts so the
# terminal condition is harmless.
recm_forward_forecast <- function(object, alg, dystar_full, n_hist, h) {
  rows <- n_hist + seq_len(h)
  if (object$expectations == "var") {
    zlag <- aux_state_lagged(dystar_full, object$aux$order)
    return(z_var(alg, object$aux$companion, object$aux$selector,
                 zlag[rows, , drop = FALSE]))
  }
  pad <- 2000L
  path <- c(dystar_full, aux_forecast_path(dystar_full, object$aux, pad))
  z_mce(alg, path, n_hist + h)[rows]
}

# Iterate the decision rule forward, recording every term of the equation
# separately. The columns of `contrib` sum to the forecast difference exactly.
recm_simulate <- function(object, fut, zf, w_mat) {
  m <- object$m
  h <- fut$h
  cf <- object$coefficients
  a <- unname(cf[seq_len(m)])
  w_terms <- object$variables$w_terms
  delta <- if (length(w_terms)) unname(cf[w_terms]) else numeric(0)

  hist <- object$model$data
  n_hist <- nrow(hist)
  y_hist <- hist[, object$variables$y]
  dy_hist <- c(NA_real_, diff(y_hist))

  lag_names <- if (m > 1L) paste0("dy_lag", seq_len(m - 1L)) else character(0)
  # `forward` is not a design column, so it can collide with an exogenous
  # regressor of that name; make.unique keeps the frames addressable.
  terms <- make.unique(c("ec", lag_names, "forward", w_terms))
  contrib <- matrix(0, h, length(terms), dimnames = list(NULL, terms))

  # dy_state[i] holds dy[t - i] for the row being formed, so at the first
  # forecast row it runs backwards from the last observed difference.
  dy_state <- if (m > 1L) {
    dy_hist[n_hist + 1L - seq_len(m - 1L)]
  } else {
    numeric(0)
  }
  if (anyNA(dy_state)) {
    stop(
      "the estimation sample ends too close to its start to supply ", m - 1L,
      " lag", if (m == 2L) "" else "s", " of the dependent variable.",
      call. = FALSE
    )
  }

  y_prev <- y_hist[n_hist]
  ystar_prev <- hist[n_hist, object$variables$y_star]
  fit <- numeric(h)
  level <- numeric(h)

  for (i in seq_len(h)) {
    contrib[i, 1L] <- a[1L] * (ystar_prev - y_prev)
    if (m > 1L) {
      contrib[i, 1L + seq_len(m - 1L)] <- a[-1L] * dy_state
    }
    contrib[i, m + 1L] <- zf[i]
    if (length(delta)) {
      contrib[i, m + 1L + seq_along(delta)] <- delta * w_mat[i, ]
    }
    step <- sum(contrib[i, ])
    y_prev <- y_prev + step
    ystar_prev <- fut$ystar[i]
    if (m > 1L) {
      dy_state <- c(step, dy_state[-(m - 1L)])
    }
    fit[i] <- step
    level[i] <- y_prev
  }

  list(contrib = contrib, fit = fit, level = level, terms = terms,
       base = unname(y_hist[n_hist]))
}

# Column-wise cumulative sum, written out rather than through apply(), which
# drops to a vector at a horizon of one.
recm_cumulate <- function(x) {
  for (j in seq_len(ncol(x))) {
    x[, j] <- cumsum(x[, j])
  }
  x
}

#' Forecast from a rational error correction model
#'
#' Iterates the estimated decision rule forward and decomposes the result,
#' term by term, into the contribution of each variable. Follows the
#' [stats::predict()] conventions where they apply: `newdata` is optional and
#' carries the future path of the right-hand side.
#'
#' @param object An object of class `"recm"` returned by [recm()].
#' @param newdata Future values of the frictionless target `y_star` and of any
#'   exogenous regressors, under the names they were estimated with. A `ts`,
#'   `data.frame` or named numeric matrix, ordered in time and with no gaps.
#'   A column holding the decision variable is ignored; see Details.
#' @param n_ahead Forecast horizon. Supplied with `newdata` it truncates it to
#'   the first `n_ahead` rows. Supplied on its own it projects the target from
#'   the fitted auxiliary autoregression, which is possible only when the fit
#'   carries no exogenous regressors, as those cannot be projected.
#' @param ... Ignored, present for consistency with the generic.
#'
#' @details
#' The forecast is a dynamic simulation, not a one step ahead conditional
#' prediction. The decision variable is iterated forward from the end of the
#' estimation sample: the error correction term at \eqn{t+h} is formed from
#' the model's own \eqn{y_{t+h-1}}, and the lagged differences from its own
#' \eqn{\Delta y}. A column of `newdata` holding the decision variable is
#' therefore ignored rather than used, because using it would make the result
#' a conditional prediction of something already known.
#'
#' The auxiliary autoregression is not re-estimated on `newdata`. The forward
#' term uses the fitted \eqn{h} vector applied to states built from the new
#' path of the target, which is what makes this a forecast from the fitted
#' model rather than from a second, differently fitted one.
#'
#' @section Contributions:
#' The equation is linear in its regressors and the coefficient on the forward
#' term is fixed at 1, so
#'
#' \deqn{\Delta y_t = a_0 \, \mathrm{ec}_t + \sum_i a_i \Delta y_{t-i} + Z_t
#'   + \delta' W_t}
#'
#' splits exactly, with nothing left over: the term columns of
#' `contributions` sum to `fit` to machine precision. Cumulating them
#' decomposes the level too, since \eqn{y_{T+h} = y_T + \sum_{s \le h}
#' \Delta y_{T+s}}, which is what `levels` holds, with the forecast origin
#' \eqn{y_T} carried in its own `base` column.
#'
#' These are contributions to the simulated path, not partial derivatives. The
#' error correction term is formed from the simulated \eqn{y}, so its column
#' carries the feedback from every term that came before it, and holding one
#' variable fixed and re-running the recursion would not move the total by
#' that variable's contribution.
#'
#' No interval is reported. The standard errors the fit carries condition on
#' the estimated auxiliary autoregression, and [summary.recm()][recm-methods]
#' says what that is worth: against a resampled target process they came in at
#' roughly a third of the true sampling standard deviation. An interval built
#' on them would inherit that and be too narrow by about as much.
#'
#' @return An object of class `"recm_forecast"`, a list whose elements are
#'   \describe{
#'     \item{`fit`}{The forecast of \eqn{\Delta y}, named by time index.}
#'     \item{`level`}{The forecast of \eqn{y} itself.}
#'     \item{`base`}{The forecast origin \eqn{y_T}, the last observation of
#'       the decision variable in the estimation sample.}
#'     \item{`contributions`}{A data frame with one row per horizon: the
#'       horizon, the time index, one column per term of the equation, and
#'       `fit`, their sum.}
#'     \item{`levels`}{The same decomposition in levels: `base`, the
#'       cumulated term columns, and `fit`, the forecast level.}
#'     \item{`forward_term`}{\eqn{Z_t} over the horizon.}
#'     \item{`terms`}{Names of the term columns of both frames.}
#'   }
#'
#' @seealso [recm()], [plot.recm_forecast()], [stats::predict.lm()].
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' dystar <- as.numeric(stats::filter(rnorm(n, 0.5, 0.4), 0.6, "recursive"))
#' ystar <- 100 + cumsum(dystar)
#' y <- numeric(n)
#' y[1] <- ystar[1]
#' for (t in 2:n) {
#'   y[t] <- y[t - 1] + 0.3 * (ystar[t - 1] - y[t - 1]) +
#'     0.7 * dystar[t] + rnorm(1, 0, 0.1)
#' }
#' fit <- recm(y, ystar, data.frame(y = y, ystar = ystar))
#'
#' # Twelve periods ahead, with the target projected from the auxiliary model.
#' fc <- predict(fit, n_ahead = 12)
#' fc
#'
#' # The contributions add up to the forecast, by construction.
#' all.equal(rowSums(fc$contributions[, fc$terms]), unname(fc$fit))
#'
#' # A target path of your own, and the same decomposition in levels.
#' future <- data.frame(ystar = ystar[n] + cumsum(rep(0.5, 8)))
#' head(predict(fit, newdata = future)$levels)
#'
#' @export
predict.recm <- function(object, newdata = NULL, n_ahead = NULL, ...) {
  if (!inherits(object, "recm")) {
    stop("`object` must be a fitted model of class \"recm\".", call. = FALSE)
  }
  if (is.null(object$model$data)) {
    stop(
      "this fit does not carry the sample the decision rule has to be ",
      "iterated forward from; it was made by an earlier version of recm(). ",
      "Refit the model.",
      call. = FALSE
    )
  }
  if (!is.null(n_ahead)) {
    if (!is.numeric(n_ahead) || length(n_ahead) != 1L || is.na(n_ahead) ||
          n_ahead != round(n_ahead) || n_ahead < 1) {
      stop("`n_ahead` must be a single positive whole number or NULL.",
           call. = FALSE)
    }
    n_ahead <- as.integer(n_ahead)
  }
  if (is.null(newdata) && is.null(n_ahead)) {
    stop(
      "nothing to forecast over. Supply `newdata` with the future path of `",
      object$variables$y_star, "`, or `n_ahead` to project that path from ",
      "the fitted auxiliary autoregression.",
      call. = FALSE
    )
  }
  if (is.null(newdata) && length(object$variables$w)) {
    stop(
      "the fit carries the exogenous regressor",
      if (length(object$variables$w) == 1L) " " else "s ",
      paste(object$variables$w, collapse = ", "),
      ", which the auxiliary autoregression cannot project. Supply `newdata` ",
      "with a future path for ",
      if (length(object$variables$w) == 1L) "it" else "each of them",
      " and for `", object$variables$y_star, "`.",
      call. = FALSE
    )
  }

  fut <- if (is.null(newdata)) {
    recm_future_projected(object, n_ahead)
  } else {
    recm_future_supplied(object, newdata, n_ahead)
  }

  hist <- object$model$data
  n_hist <- nrow(hist)
  ystar_hist <- hist[, object$variables$y_star]
  # One contiguous path of differences across the join, so the first forecast
  # difference is taken against the last observed level of the target.
  dystar_full <- c(
    NA_real_, diff(c(ystar_hist, fut$ystar))
  )

  alg <- pac_algebra(recm_a(object), object$discount)
  zf <- recm_forward_forecast(object, alg, dystar_full, n_hist, fut$h)
  w_mat <- recm_future_w(object, fut)
  sim <- recm_simulate(object, fut, zf, w_mat)

  lab <- as.character(fut$index)
  contributions <- data.frame(
    horizon = seq_len(fut$h), time = fut$index, sim$contrib,
    fit = sim$fit, check.names = FALSE, stringsAsFactors = FALSE
  )
  levels_df <- data.frame(
    horizon = seq_len(fut$h), time = fut$index, base = sim$base,
    recm_cumulate(sim$contrib), fit = sim$level,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  # Guard the one added name against a regressor that already carries it.
  base_name <- make.unique(c(sim$terms, "base"))[length(sim$terms) + 1L]
  names(levels_df)[3L] <- base_name

  structure(
    list(
      call = match.call(),
      horizon = fut$h,
      fit = stats::setNames(sim$fit, lab),
      level = stats::setNames(sim$level, lab),
      base = sim$base,
      base_name = base_name,
      contributions = contributions,
      levels = levels_df,
      forward_term = stats::setNames(as.numeric(zf), lab),
      terms = sim$terms,
      time = fut$index,
      origin = object$model$index[n_hist],
      projected = fut$projected,
      expectations = object$expectations,
      m = object$m,
      variables = object$variables
    ),
    class = "recm_forecast"
  )
}

#' @rdname predict.recm
#' @param x An object of class `"recm_forecast"`.
#' @param digits Number of significant digits used when printing.
#' @export
print.recm_forecast <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("Forecast from a rational error correction model\n\n")
  cat("Call:\n  ", paste(deparse(x$call), collapse = "\n  "), "\n\n", sep = "")
  cat(sprintf(
    "  horizon: %d   origin: %s   expectations: %s\n",
    x$horizon, format(x$origin), x$expectations
  ))
  cat(sprintf(
    "  target %s: %s\n",
    x$variables$y_star,
    if (x$projected) {
      "projected from the auxiliary autoregression"
    } else {
      "supplied in `newdata`"
    }
  ))
  cat(
    "\nContributions to d", x$variables$y,
    " (the term columns sum to fit):\n", sep = ""
  )
  print(x$contributions, digits = digits, row.names = FALSE)
  cat("\nThe path is simulated, so the `ec` column carries the feedback from",
      "\nthe terms before it rather than a partial derivative.\n")
  invisible(x)
}
