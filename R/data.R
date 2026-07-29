## Data handling: resolving the y / y_star arguments against `data`, coercing
## the supported container types to a numeric matrix, and building the design
## of the rational error correction equation.

# Resolve an unevaluated `y` / `y_star` argument to a column name of `data`.
#
# Accepts a bare symbol (`gdp`), a character scalar (`"gdp"`), or an expression
# in the caller's environment that evaluates to a character scalar naming a
# column. A column name always wins over a same-named object in the caller.
resolve_column <- function(expr, env, cols, arg) {
  nm <- NULL
  if (is.symbol(expr)) {
    nm <- as.character(expr)
  } else if (is.character(expr) && length(expr) == 1L) {
    nm <- expr
  }
  if (!is.null(nm) && nm %in% cols) {
    return(nm)
  }
  val <- tryCatch(eval(expr, env), error = function(e) NULL)
  if (is.character(val) && length(val) == 1L && val %in% cols) {
    return(val)
  }
  stop(
    "`", arg, "` (", deparse(expr), ") does not name a column of `data`. ",
    "Columns available: ", paste(cols, collapse = ", "), ".",
    call. = FALSE
  )
}

# Coerce a ts / mts, data.frame (including tibble and data.table), or matrix to
# a list holding a numeric matrix and a time index.
#
# At most one non-numeric column is tolerated in a data frame; it is taken as
# the time index rather than as a regressor. Anything more is an error, because
# silently dropping columns would silently drop regressors.
as_model_data <- function(data) {
  if (stats::is.ts(data)) {
    x <- as.matrix(data)
    if (is.null(colnames(x))) {
      stop(
        "`data` is a time series with no column names; ",
        "`y` and `y_star` cannot be resolved.",
        call. = FALSE
      )
    }
    return(list(x = x, index = as.numeric(stats::time(data)),
                index_name = "time"))
  }

  if (is.data.frame(data)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE)
    if (!ncol(df)) {
      stop("`data` has no columns.", call. = FALSE)
    }
    numeric_col <- vapply(df, is.numeric, logical(1))
    other <- names(df)[!numeric_col]
    if (length(other) > 1L) {
      stop(
        "`data` has more than one non-numeric column (",
        paste(other, collapse = ", "),
        "). Keep at most one (a date or period label) and convert or drop ",
        "the rest; every remaining numeric column is treated as an ",
        "exogenous regressor.",
        call. = FALSE
      )
    }
    if (!any(numeric_col)) {
      stop("`data` has no numeric columns.", call. = FALSE)
    }
    index <- if (length(other)) df[[other]] else seq_len(nrow(df))
    return(list(x = as.matrix(df[numeric_col]), index = index,
                index_name = if (length(other)) other else NULL))
  }

  if (is.matrix(data)) {
    if (!is.numeric(data)) {
      stop("`data` is a matrix but is not numeric.", call. = FALSE)
    }
    if (is.null(colnames(data))) {
      stop(
        "`data` is a matrix with no column names; ",
        "`y` and `y_star` cannot be resolved.",
        call. = FALSE
      )
    }
    return(list(x = data, index = seq_len(nrow(data)), index_name = NULL))
  }

  stop(
    "`data` must be a ts, data.frame (tibble and data.table qualify), ",
    "or a named numeric matrix; got ", class(data)[1L], ".",
    call. = FALSE
  )
}

# Build the regressors of the rational error correction equation.
#
# The estimated equation is
#
#   dy[t] = a0 * (ystar[t-1] - y[t-1]) + sum_i a[i] * dy[t-i] + Z[t] + d'W[t]
#
# so the design is the error correction term, m-1 lags of dy, and the columns
# of W dated t. No intercept: a free constant is inconsistent with growth
# neutrality, and the auxiliary autoregression carries the drift instead.
#
# `tr_exog` differences W as y and y_star are differenced, and renames the
# column so the coefficient label states the transform. The dependent variable
# is a difference, so a level regressor is a different object from the rest of
# the equation: were W trending it would put a trend into dy and break the
# growth neutrality the restriction is there to enforce. Differencing costs no
# observations, since the first row is already lost to dy.
build_design <- function(x, y_name, ystar_name, m, tr_exog) {
  n <- nrow(x)
  y <- x[, y_name]
  ystar <- x[, ystar_name]
  w_names <- setdiff(colnames(x), c(y_name, ystar_name))
  # Guarded on length: paste0("d_", character(0)) recycles to "d_" rather than
  # returning nothing, which would invent a regressor when there are none.
  w_terms <- if (tr_exog && length(w_names)) paste0("d_", w_names) else w_names

  lag_vec <- function(v, k) c(rep(NA_real_, k), v[seq_len(n - k)])

  dy <- c(NA_real_, diff(y))
  dystar <- c(NA_real_, diff(ystar))
  ec <- lag_vec(ystar, 1L) - lag_vec(y, 1L)

  cols <- list(ec)
  nms <- "ec"
  if (m > 1L) {
    for (i in seq_len(m - 1L)) {
      cols[[length(cols) + 1L]] <- lag_vec(dy, i)
      nms <- c(nms, paste0("dy_lag", i))
    }
  }
  for (k in seq_along(w_names)) {
    wk <- x[, w_names[k]]
    cols[[length(cols) + 1L]] <- if (tr_exog) c(NA_real_, diff(wk)) else wk
    nms <- c(nms, w_terms[k])
  }

  # A duplicate here would make coef()[["name"]] silently return the wrong
  # column, so name the offender rather than let it through.
  if (anyDuplicated(nms)) {
    stop(
      "the design has duplicate column names (",
      paste(unique(nms[duplicated(nms)]), collapse = ", "),
      "). Rename the offending column of `data`: the error correction term is ",
      "called `ec`, the lags of the dependent variable `dy_lag1` and so on, ",
      "and with `tr_exog = TRUE` each exogenous regressor is prefixed `d_`.",
      call. = FALSE
    )
  }

  design <- matrix(unlist(cols), nrow = n, ncol = length(cols))
  colnames(design) <- nms

  list(dy = dy, dystar = dystar, x = design, w_names = w_names,
       w_terms = w_terms, y = y, ystar = ystar)
}
