## The simulation entry point, and the shocks it can be driven with.
##
## An impulse response is computed as the difference between two runs of the
## recursion in R/path.R: one over a baseline path, one over the same path with
## the shock added to it. Nothing about the fit changes between them.
##
## Two runs rather than a closed form because the decision rule is exactly
## linear in its inputs, which makes the difference the impulse response
## itself, computed by the same code that is already tested for the level. It
## also makes the difference exact term by term, so the decomposition carries
## over to the response for free, and it makes the answer independent of the
## baseline - which is why a baseline can be invented here (target projected
## from the auxiliary autoregression, exogenous regressors held flat) in cases
## where predict() rightly refuses to invent a forecast.

#' Specify a shock for an impulse response
#'
#' Describes a perturbation to the path of one of the inputs of the decision
#' rule, for [simulate.recm()] to trace the response to. The shock is always
#' stated in the units and the level of the variable it names: `type` and
#' `period` say when the level is displaced and for how long, not how the
#' displacement is transformed on its way into the equation.
#'
#' @param variable Name of the variable shocked, either the frictionless
#'   target `y_star` or one of the exogenous regressors, as named in the fit.
#'   The decision variable itself cannot be shocked: it is produced by the
#'   recursion rather than supplied to it.
#' @param size Size of the displacement, in the units of `variable`. A single
#'   number, or, with `type = "path"`, the whole displacement path.
#' @param period Horizon at which the shock lands, counting the first simulated
#'   period as 1.
#' @param type `"permanent"` holds the level `size` higher from `period`
#'   onwards; `"transitory"` displaces it in `period` alone and returns it to
#'   baseline after; `"path"` takes `size` as an arbitrary displacement path
#'   beginning at `period`.
#'
#' @return An object of class `"recm_shock"`.
#'
#' @seealso [simulate.recm()].
#'
#' @examples
#' # A permanent one unit rise in the target.
#' recm_shock("ystar")
#'
#' # A blip, three periods in.
#' recm_shock("ystar", size = 2, period = 3, type = "transitory")
#'
#' # An arbitrary hump.
#' recm_shock("ystar", size = c(0.5, 1, 0.5), type = "path")
#'
#' @export
recm_shock <- function(variable, size = 1, period = 1,
                       type = c("permanent", "transitory", "path")) {
  type <- match.arg(type)
  if (!is.character(variable) || length(variable) != 1L || is.na(variable) ||
        !nzchar(variable)) {
    stop("`variable` must be a single non-empty variable name.", call. = FALSE)
  }
  if (!is.numeric(size) || !length(size) || anyNA(size) ||
        !all(is.finite(size))) {
    stop("`size` must be a finite numeric vector.", call. = FALSE)
  }
  if (type != "path" && length(size) != 1L) {
    stop(
      "`size` must be a single number unless `type = \"path\"`, which is the ",
      "way to give a displacement path of your own.",
      call. = FALSE
    )
  }
  if (!is.numeric(period) || length(period) != 1L || is.na(period) ||
        period != round(period) || period < 1) {
    stop("`period` must be a single positive whole number.", call. = FALSE)
  }
  structure(
    list(variable = variable, size = as.numeric(size),
         period = as.integer(period), type = type),
    class = "recm_shock"
  )
}

# A shock in one line, for the print methods and the plot titles.
recm_shock_label <- function(shock) {
  if (shock$type == "path") {
    return(paste0(shock$variable, ", path of ", length(shock$size),
                  " periods, from period ", shock$period))
  }
  paste0(shock$variable, " ", if (shock$size >= 0) "+" else "",
         format(shock$size, digits = 4), ", ", shock$type,
         ", from period ", shock$period)
}

#' @rdname recm_shock
#' @param x An object of class `"recm_shock"`.
#' @param ... Ignored, present for consistency with the generic.
#' @export
print.recm_shock <- function(x, ...) {
  cat("Shock: ", recm_shock_label(x), "\n", sep = "")
  invisible(x)
}

# Accept either a shock object or the named-number shorthand, and check the
# variable against the fit, which recm_shock() cannot see.
recm_resolve_shock <- function(object, shock) {
  if (!inherits(shock, "recm_shock")) {
    ok <- is.numeric(shock) && length(shock) == 1L &&
      !is.null(names(shock)) && nzchar(names(shock))
    if (!ok) {
      stop(
        "`shock` must be an object from recm_shock(), or a single named ",
        "number such as `c(", object$variables$y_star, " = 1)` for a ",
        "permanent unit shock.",
        call. = FALSE
      )
    }
    shock <- recm_shock(names(shock), unname(shock))
  }
  allowed <- c(object$variables$y_star, object$variables$w)
  if (!shock$variable %in% allowed) {
    stop(
      "`", shock$variable, "` is not an input of this decision rule, so a ",
      "shock to it has nowhere to enter. Shockable: ",
      paste(allowed, collapse = ", "),
      if (identical(shock$variable, object$variables$y)) {
        paste0(". `", object$variables$y, "` is produced by the recursion ",
               "rather than supplied to it")
      } else {
        ""
      },
      ".",
      call. = FALSE
    )
  }
  shock
}

# The displacement of the shocked variable's level, period by period.
recm_shock_path <- function(shock, h) {
  if (shock$period > h) {
    stop(
      "the shock lands at period ", shock$period, " but the horizon is ", h,
      ", so nothing would happen. Lengthen the horizon or move the shock.",
      call. = FALSE
    )
  }
  dev <- numeric(h)
  if (shock$type == "permanent") {
    dev[shock$period:h] <- shock$size
  } else if (shock$type == "transitory") {
    dev[shock$period] <- shock$size
  } else {
    last <- shock$period + length(shock$size) - 1L
    if (last > h) {
      stop(
        "the shock path runs to period ", last, " but the horizon is ", h,
        ". Lengthen the horizon, or shorten `size`.",
        call. = FALSE
      )
    }
    dev[shock$period + seq_along(shock$size) - 1L] <- shock$size
  }
  dev
}

# The baseline future with the displacement added to the shocked input.
recm_shift_future <- function(object, fut, shock) {
  dev <- recm_shock_path(shock, fut$h)
  if (identical(shock$variable, object$variables$y_star)) {
    fut$ystar <- fut$ystar + dev
  } else {
    fut$w[, shock$variable] <- fut$w[, shock$variable] + dev
  }
  fut
}

# Two runs, differenced. Exact term by term because the recursion is linear in
# its inputs, so the term columns of the response add up as the levels do.
recm_difference <- function(shocked, baseline) {
  list(
    contrib = shocked$contrib - baseline$contrib,
    fit = shocked$fit - baseline$fit,
    level = shocked$level - baseline$level,
    forward = shocked$forward - baseline$forward,
    terms = shocked$terms,
    # The response starts from nothing, so there is no level to carry.
    base = 0
  )
}

#' Dynamic simulation and impulse responses of the decision rule
#'
#' Iterates the estimated decision rule forward over a path of your choosing
#' and decomposes the result, term by term, into the contribution of each
#' variable. Given a `shock`, reports the impulse response: how the path
#' differs from the one the same rule would have produced without it.
#'
#' @param object An object of class `"recm"` returned by [recm()].
#' @param nsim Number of paths. Only `1` is accepted: the path is the
#'   deterministic one, the decision rule iterated forward with the equation
#'   error set to zero, so a second replicate would be identical to the first.
#'   The argument is present because [stats::simulate()] defines it.
#' @param seed Not used, and refused when supplied rather than silently
#'   ignored, because nothing here is random. Present for the same reason as
#'   `nsim`.
#' @param newdata Values of the frictionless target `y_star` and of any
#'   exogenous regressors over the simulated periods, under the names they were
#'   estimated with. A `ts`, `data.frame` or named numeric matrix, ordered in
#'   time and with no gaps. A column holding the decision variable is ignored;
#'   see Details.
#' @param n_ahead Length of the simulation. Supplied with `newdata` it
#'   truncates it to the first `n_ahead` rows. Supplied on its own it projects
#'   the target from the fitted auxiliary autoregression.
#' @param shock A shock to trace the response to, from [recm_shock()], or the
#'   shorthand `c(ystar = 1)` for a permanent unit shock. `NULL`, the default,
#'   simulates the path itself rather than a response to anything.
#' @param ... Ignored, present for consistency with the generic.
#'
#' @details
#' The simulation is dynamic. The decision variable is iterated forward from
#' the end of the estimation sample: the error correction term at \eqn{t+h} is
#' formed from the model's own \eqn{y_{t+h-1}}, and the lagged differences from
#' its own \eqn{\Delta y}. A column of `newdata` holding the decision variable
#' is therefore ignored rather than used.
#'
#' @section Impulse responses:
#' With a `shock`, the recursion is run twice - once over the baseline path,
#' once over the same path with the shock added to the level of the variable it
#' names - and `fit` and `level` report the difference. The two underlying
#' paths are kept in `baseline` and `shocked`, each an ordinary simulation
#' object in its own right.
#'
#' The response does not depend on the baseline it was measured against. The
#' decision rule is linear in its inputs and the forward term is linear in the
#' states built from the target, so the difference between the two runs is a
#' function of the shock alone. That is why a shocked simulation will project
#' a baseline for itself, holding any exogenous regressors at their last
#' observed level, in cases where [predict.recm()] refuses to project them: as
#' a forecast that flat path would be a fiction, but as the baseline of a
#' response it cancels exactly.
#'
#' Timing follows the model's information set, and the two expectation branches
#' differ because of it. Under `expectations = "var"` the forward term at
#' \eqn{t} is built from states dated \eqn{t-1}, so a shock landing in period
#' \eqn{k} moves nothing until period \eqn{k+1}. Under `"mce"` the forward term
#' is evaluated on the realised path and the response begins in period \eqn{k}
#' itself.
#'
#' @return An object of class `"recm_simulation"`, a list with the same
#'   elements as a [predict.recm()] forecast - `fit`, `level`, `base`,
#'   `contributions`, `levels`, `forward_term`, `terms` - together with
#'   \describe{
#'     \item{`shock`}{The shock traced, or `NULL`.}
#'     \item{`baseline`, `shocked`}{With a shock, the two simulations whose
#'       difference is being reported. Absent otherwise.}
#'   }
#'   With a shock, `fit` and `level` are responses rather than paths: `fit` is
#'   the deviation of \eqn{\Delta y} from baseline period by period, `level`
#'   the deviation of \eqn{y}, and `base` is `0`.
#'
#' @seealso [recm_shock()] for the shock, [plot.recm_simulation()] for the
#'   picture, and [predict.recm()] for the same recursion put to the other
#'   purpose: forecasting \eqn{y} over the next \eqn{h} periods rather than
#'   exploring how the rule behaves.
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
#' # The path itself, over a target of your choosing.
#' sim <- simulate(fit, newdata = data.frame(ystar = ystar[n] + 0.5 * (1:20)))
#' sim
#'
#' # The response to a permanent unit rise in the target. Growth neutrality
#' # means y takes all of it eventually, so the response converges on 1.
#' irf <- simulate(fit, n_ahead = 60, shock = recm_shock("ystar"))
#' round(irf$level[c(1, 5, 10, 20, 60)], 3)
#'
#' # A transitory shock is undone, so the level response returns to zero.
#' blip <- simulate(fit, n_ahead = 60,
#'                  shock = recm_shock("ystar", type = "transitory"))
#' round(blip$level[c(1, 5, 60)], 3)
#'
#' @export
simulate.recm <- function(object, nsim = 1, seed = NULL, newdata = NULL,
                          n_ahead = NULL, shock = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim != 1) {
    stop(
      "`nsim` must be 1. The simulated path is the deterministic one, the ",
      "decision rule iterated forward with the equation error set to zero, ",
      "so a second replicate would be identical to the first. Drawing the ",
      "equation error from the residuals is not implemented.",
      call. = FALSE
    )
  }
  if (!is.null(seed)) {
    stop(
      "`seed` has no effect, so it is refused rather than ignored: the ",
      "simulated path is deterministic. The argument is in the signature ",
      "only because stats::simulate() defines it.",
      call. = FALSE
    )
  }
  call <- match.call()
  alg <- recm_algebra(object)

  if (is.null(shock)) {
    fut <- recm_future(object, newdata, n_ahead)
    run <- recm_run(object, alg, fut)
    return(recm_assemble(object, fut, run, call, "recm_simulation",
                         list(shock = NULL)))
  }

  shock <- recm_resolve_shock(object, shock)
  # hold_w: an invented baseline is admissible here and only here, because the
  # response is measured against it and is invariant to it.
  fut <- recm_future(object, newdata, n_ahead, hold_w = TRUE)
  fut_shocked <- recm_shift_future(object, fut, shock)
  base_run <- recm_run(object, alg, fut)
  shock_run <- recm_run(object, alg, fut_shocked)

  recm_assemble(
    object, fut, recm_difference(shock_run, base_run), call, "recm_simulation",
    list(
      shock = shock,
      baseline = recm_assemble(object, fut, base_run, call, "recm_simulation",
                               list(shock = NULL)),
      shocked = recm_assemble(object, fut_shocked, shock_run, call,
                              "recm_simulation", list(shock = NULL))
    )
  )
}

#' @rdname simulate.recm
#' @param x An object of class `"recm_simulation"`.
#' @param digits Number of significant digits used when printing.
#' @export
print.recm_simulation <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  shocked <- !is.null(x$shock)
  if (shocked) {
    cat("Impulse response of a rational error correction decision rule\n\n")
  } else {
    cat("Dynamic simulation of a rational error correction decision rule\n\n")
  }
  target <- sprintf(
    "target %s: %s", x$variables$y_star,
    if (x$projected) {
      "projected from the auxiliary autoregression"
    } else {
      "supplied in `newdata`"
    }
  )
  lines <- if (shocked) {
    c(paste0("shock: ", recm_shock_label(x$shock)), target)
  } else {
    target
  }
  what <- if (shocked) {
    paste0("the response of d", x$variables$y)
  } else {
    paste0("d", x$variables$y)
  }
  recm_print_body(x, digits, what, lines)
  if (shocked) {
    cat("Reported as the deviation from the baseline in `$baseline`, which",
        "\nthe response does not depend on: the decision rule is linear.\n")
  }
  invisible(x)
}
