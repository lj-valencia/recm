## The machinery both entry points are built on: iterating the estimated
## decision rule forward from the end of the sample, and splitting the path it
## produces into the contribution of each term of the equation.
##
## Nothing here is exported. predict.recm() wraps it as an h step ahead
## forecast and simulate.recm() wraps it as a dynamic simulation, optionally
## differencing two runs against each other to form an impulse response. The
## split is deliberate: the arithmetic is one thing and is tested once, while
## what the two entry points mean is a separate question answered in their own
## files.
##
## The path is a dynamic one, not a sequence of one step ahead conditional
## predictions: y is iterated forward, so the error correction term at t + h is
## formed from the model's own y[t + h - 1] and never from a realised one. The
## future data therefore supply the target and the exogenous regressors only.
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
## correction term is formed from the iterated y, so its bar carries the
## feedback from every term that came before it: holding one column fixed and
## re-running the recursion would not move the total by that column's
## contribution.
##
## That the whole map from inputs to path is *linear* is worth stating on its
## own, because an impulse response leans on it entirely. Every term above is
## linear in (ystar path, W path, initial y, initial lags of dy), and the
## forward term is linear in the states built from the target under both
## expectation branches. The difference between a shocked run and a baseline
## one therefore depends on the shock alone, and not on the baseline it was
## measured against.

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
#
# `hold_w` keeps the exogenous regressors at their last observed level. That is
# not a forecast of them and is never offered as one - it is used only as the
# baseline of an impulse response, where the linearity noted at the top of this
# file makes the baseline irrelevant to the answer.
recm_future_projected <- function(object, n_ahead, hold_w = FALSE) {
  hist <- object$model$data
  n_hist <- nrow(hist)
  ystar <- hist[, object$variables$y_star]
  dfut <- aux_forecast_path(c(NA_real_, diff(ystar)), object$aux, n_ahead)
  w_names <- object$variables$w
  w <- if (hold_w && length(w_names)) {
    matrix(
      rep(hist[n_hist, w_names], each = n_ahead), n_ahead, length(w_names),
      dimnames = list(NULL, w_names)
    )
  } else {
    NULL
  }
  list(
    h = n_ahead,
    ystar = ystar[n_hist] + cumsum(dfut),
    w = w,
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
      "`newdata` has no column `", ys_name, "`. The decision rule is iterated ",
      "forward, so the future path of the target has to be supplied; only `",
      object$variables$y, "` is produced rather than read. Columns found: ",
      paste(cols, collapse = ", "), ".",
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
      " rows of `newdata` that are used; the path must be complete.",
      call. = FALSE
    )
  }
  w <- if (length(w_names)) nd$x[rows, w_names, drop = FALSE] else NULL
  if (!is.null(w) && anyNA(w)) {
    stop(
      "the exogenous regressors are missing in ",
      sum(!stats::complete.cases(w)), " of the ", h,
      " rows of `newdata` that are used; the paths must be complete.",
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

# `n_ahead` as a whole number, or NULL when the horizon comes from `newdata`.
recm_check_horizon <- function(n_ahead) {
  if (is.null(n_ahead)) {
    return(NULL)
  }
  if (!is.numeric(n_ahead) || length(n_ahead) != 1L || is.na(n_ahead) ||
        n_ahead != round(n_ahead) || n_ahead < 1) {
    stop("`n_ahead` must be a single positive whole number or NULL.",
         call. = FALSE)
  }
  as.integer(n_ahead)
}

# Everything the recursion needs to know about the future, from whichever
# source, with the fitted object checked over first.
recm_future <- function(object, newdata, n_ahead, hold_w = FALSE) {
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
  n_ahead <- recm_check_horizon(n_ahead)
  if (is.null(newdata) && is.null(n_ahead)) {
    stop(
      "nothing to iterate over. Supply `newdata` with the future path of `",
      object$variables$y_star, "`, or `n_ahead` to project that path from ",
      "the fitted auxiliary autoregression.",
      call. = FALSE
    )
  }
  if (is.null(newdata) && !hold_w && length(object$variables$w)) {
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

  if (is.null(newdata)) {
    recm_future_projected(object, n_ahead, hold_w)
  } else {
    recm_future_supplied(object, newdata, n_ahead)
  }
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

# The forward term over the future rows, at the fitted coefficients and the
# fitted auxiliary autoregression. Under "var" it is h'z[t-1] on states built
# from the new path of the target; under "mce" it is the backward recursion on
# that path, extended past the horizon by autoregressive forecasts so the
# terminal condition is harmless.
recm_forward_path <- function(object, alg, dystar_full, n_hist, h) {
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

# The recursion itself: step the decision rule forward one period at a time,
# recording every term of the equation separately. The columns of `contrib`
# sum to the path of differences exactly.
recm_iterate <- function(object, fut, zf, w_mat) {
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
  # future row it runs backwards from the last observed difference.
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
       forward = as.numeric(zf), base = unname(y_hist[n_hist]))
}

# The polynomial algebra implied by the fit. Computed once per call and passed
# down, because an impulse response makes two passes over the recursion and
# this part of it depends on the fit alone, not on the path.
recm_algebra <- function(object) {
  pac_algebra(recm_a(object), object$discount)
}

# One complete pass of the decision rule over `fut`.
recm_run <- function(object, alg, fut) {
  hist <- object$model$data
  n_hist <- nrow(hist)
  # One contiguous path of differences across the join, so the first future
  # difference is taken against the last observed level of the target.
  dystar_full <- c(
    NA_real_, diff(c(hist[, object$variables$y_star], fut$ystar))
  )
  zf <- recm_forward_path(object, alg, dystar_full, n_hist, fut$h)
  recm_iterate(object, fut, zf, recm_future_w(object, fut))
}

# Column-wise cumulative sum, written out rather than through apply(), which
# drops to a vector at a horizon of one.
recm_cumulate <- function(x) {
  for (j in seq_len(ncol(x))) {
    x[, j] <- cumsum(x[, j])
  }
  x
}

# Wrap a run in the object the user is handed. `cls` and `extra` are what the
# two entry points differ by; everything else about the result is shared, and
# so is the guarantee that the term columns of both frames add up.
recm_assemble <- function(object, fut, run, call, cls, extra = list()) {
  lab <- as.character(fut$index)
  contributions <- data.frame(
    horizon = seq_len(fut$h), time = fut$index, run$contrib,
    fit = run$fit, check.names = FALSE, stringsAsFactors = FALSE
  )
  levels_df <- data.frame(
    horizon = seq_len(fut$h), time = fut$index, base = run$base,
    recm_cumulate(run$contrib), fit = run$level,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  # Guard the one added name against a regressor that already carries it.
  base_name <- make.unique(c(run$terms, "base"))[length(run$terms) + 1L]
  names(levels_df)[3L] <- base_name

  structure(
    c(
      list(
        call = call,
        horizon = fut$h,
        fit = stats::setNames(run$fit, lab),
        level = stats::setNames(run$level, lab),
        base = run$base,
        base_name = base_name,
        contributions = contributions,
        levels = levels_df,
        forward_term = stats::setNames(run$forward, lab),
        terms = run$terms,
        time = fut$index,
        origin = object$model$index[nrow(object$model$data)],
        projected = fut$projected,
        expectations = object$expectations,
        m = object$m,
        variables = object$variables
      ),
      extra
    ),
    class = cls
  )
}

# The body both print methods share: the call, the one line summary, whatever
# each of them wants to add to it, and the decomposition itself.
recm_print_body <- function(x, digits, what, extra_lines = character(0)) {
  cat("Call:\n  ", paste(deparse(x$call), collapse = "\n  "), "\n\n", sep = "")
  cat(sprintf(
    "  horizon: %d   origin: %s   expectations: %s\n",
    x$horizon, format(x$origin), x$expectations
  ))
  for (line in extra_lines) {
    cat("  ", line, "\n", sep = "")
  }
  cat("\nContributions to ", what, " (the term columns sum to fit):\n",
      sep = "")
  print(x$contributions, digits = digits, row.names = FALSE)
  cat("\nThe path is iterated forward, so the `ec` column carries the feedback",
      "\nfrom the terms before it rather than a partial derivative.\n")
  invisible(NULL)
}
