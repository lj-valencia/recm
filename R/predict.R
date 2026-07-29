## The forecasting entry point. Everything it does is in R/path.R; what is
## here is what makes it a forecast rather than a simulation.
##
## The distinction against simulate.recm() is one of purpose and not of
## arithmetic. predict() answers "where is y going", so it runs from the end of
## the estimation sample to a horizon h, it needs a real path for every input
## it cannot project, and it reports the path itself. simulate() answers "how
## does the decision rule behave", so it will happily be handed a path that
## never happened, and given a shock it reports a deviation rather than a
## level.

#' h-step ahead forecasts from a rational error correction model
#'
#' Forecasts the decision variable over the next `n_ahead` periods by iterating
#' the estimated decision rule forward from the end of the estimation sample,
#' and decomposes the forecast, term by term, into the contribution of each
#' variable.
#'
#' @param object An object of class `"recm"` returned by [recm()].
#' @param newdata Future values of the frictionless target `y_star` and of any
#'   exogenous regressors, under the names they were estimated with. A `ts`,
#'   `data.frame` or named numeric matrix, ordered in time and with no gaps.
#'   `NULL`, the default, projects the target from the fitted auxiliary
#'   autoregression, which is possible only when the fit carries no exogenous
#'   regressors, as those cannot be projected. A column holding the decision
#'   variable is ignored; see Details.
#' @param n_ahead Forecast horizon \eqn{h}, a single positive whole number.
#'   Defaults to `1`, the one step ahead forecast. When `newdata` is supplied
#'   and `n_ahead` is not, the horizon is the number of rows of `newdata`;
#'   supplying both truncates `newdata` to its first `n_ahead` rows.
#' @param ... Ignored, present for consistency with the generic.
#'
#' @details
#' The forecast is dynamic, not a sequence of one step ahead conditional
#' predictions. The decision variable is iterated forward: the error correction
#' term at \eqn{t+h} is formed from the model's own \eqn{y_{t+h-1}}, and the
#' lagged differences from its own \eqn{\Delta y}. A column of `newdata`
#' holding the decision variable is therefore ignored rather than used, because
#' using it would make the result a prediction of something already known, and
#' `predict(object)` does not reproduce [fitted()][recm-methods] in sample.
#'
#' Leaving `newdata` out gives the model's own forecast, with the target
#' projected from the auxiliary autoregression. Supplying it gives a
#' conditional forecast, on the path of the target you assume. Either way the
#' auxiliary autoregression is not re-estimated: the forward term uses the
#' fitted \eqn{h} vector applied to states built from the new path of the
#' target, which is what makes this a forecast from the fitted model rather
#' than from a second, differently fitted one.
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
#' These are contributions to the forecast, not partial derivatives. The error
#' correction term is formed from the forecast \eqn{y}, so its column carries
#' the feedback from every term that came before it, and holding one variable
#' fixed and re-running the recursion would not move the total by that
#' variable's contribution.
#'
#' No interval is reported. The standard errors the fit carries condition on
#' the estimated auxiliary autoregression, and [summary.recm()][recm-methods]
#' says what that is worth: against a resampled target process they came in at
#' roughly a third of the true sampling standard deviation. An interval built
#' on them would inherit that and be too narrow by about as much.
#'
#' @return An object of class `"recm_forecast"`, a list whose elements are
#'   \describe{
#'     \item{`fit`}{The forecast path of \eqn{\Delta y}, named by time index.}
#'     \item{`level`}{The forecast path of \eqn{y} itself.}
#'     \item{`base`}{The forecast origin \eqn{y_T}, the last observation of the
#'       decision variable in the estimation sample.}
#'     \item{`contributions`}{A data frame with one row per period: the
#'       horizon, the time index, one column per term of the equation, and
#'       `fit`, their sum.}
#'     \item{`levels`}{The same decomposition in levels: `base`, the
#'       cumulated term columns, and `fit`, the forecast level.}
#'     \item{`forward_term`}{\eqn{Z_t} over the horizon.}
#'     \item{`terms`}{Names of the term columns of both frames.}
#'     \item{`projected`}{Whether the target was projected from the auxiliary
#'       autoregression rather than supplied.}
#'   }
#'
#' @seealso [simulate.recm()], which iterates the same decision rule but for a
#'   different purpose: over a path of your choosing, and with a shock, to
#'   trace an impulse response rather than to forecast. [plot.recm_forecast()],
#'   [recm()].
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
#' # One step ahead, the default.
#' predict(fit)$fit
#'
#' # Twelve periods, with the target projected from the auxiliary model.
#' fc <- predict(fit, n_ahead = 12)
#' fc
#'
#' # The contributions add up to the forecast, by construction.
#' all.equal(rowSums(fc$contributions[, fc$terms]), unname(fc$fit))
#'
#' # A conditional forecast, on a target path you assume.
#' future <- data.frame(ystar = ystar[n] + cumsum(rep(0.5, 8)))
#' head(predict(fit, newdata = future)$levels)
#'
#' @export
predict.recm <- function(object, newdata = NULL, n_ahead = 1L, ...) {
  # `newdata` on its own means "as far as newdata goes", which is the useful
  # reading and not what a default of 1 would give; only an n_ahead the caller
  # actually typed truncates it.
  horizon <- if (missing(n_ahead) && !is.null(newdata)) NULL else n_ahead
  fut <- recm_future(object, newdata, horizon)
  run <- recm_run(object, recm_algebra(object), fut)
  recm_assemble(object, fut, run, match.call(), "recm_forecast")
}

#' @rdname predict.recm
#' @param x An object of class `"recm_forecast"`.
#' @param digits Number of significant digits used when printing.
#' @export
print.recm_forecast <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("h-step ahead forecast from a rational error correction model\n\n")
  recm_print_body(
    x, digits, paste0("d", x$variables$y),
    sprintf(
      "target %s: %s",
      x$variables$y_star,
      if (x$projected) {
        "projected from the auxiliary autoregression"
      } else {
        "supplied in `newdata`, so the forecast is conditional on it"
      }
    )
  )
  invisible(x)
}
