test_that("recm recovers the coefficients of its own data generating process", {
  df <- simulate_recm(
    n = 800L, a = 0.3, beta = 1, ar = 0.6, const = 0.2, seed = 4
  )
  fit <- recm(y, ystar, df)
  expect_s3_class(fit, "recm")
  expect_true(fit$converged)
  # Monte Carlo over 40 replications puts the standard deviation of this
  # coefficient near 0.013 at n = 800, so 0.05 is roughly four of them.
  expect_equal(unname(coef(fit)[["ec"]]), 0.3, tolerance = 0.05)
  # At beta = 1 the total forward loading is 1 whatever a0 is.
  expect_equal(fit$d_sum, 1, tolerance = 1e-10)
})

test_that("recm recovers a model with an autoregressive term", {
  df <- simulate_recm(
    n = 1500L, a = c(0.35, 0.25), beta = 1, ar = 0.5, const = 0.2, seed = 7
  )
  fit <- recm(y, ystar, df, m = 2)
  expect_true(fit$converged)
  expect_equal(unname(coef(fit)[["ec"]]), 0.35, tolerance = 0.08)
  expect_equal(unname(coef(fit)[["dy_lag1"]]), 0.25, tolerance = 0.2)
})

test_that("the mce branch recovers a model consistent process", {
  df <- simulate_recm(
    n = 800L, a = 0.3, beta = 1, ar = 0.6, const = 0.2,
    expectations = "mce", sd_eq = 0.02, seed = 8
  )
  fit <- recm(y, ystar, df, expectations = "mce")
  expect_true(fit$converged)
  expect_equal(unname(coef(fit)[["ec"]]), 0.3, tolerance = 0.05)
})

test_that("y and y_star are accepted as objects and as strings", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 2)
  by_symbol <- recm(y, ystar, df)
  by_string <- recm("y", "ystar", df)
  expect_equal(coef(by_symbol), coef(by_string))

  # A character variable in the caller naming a column also resolves.
  target <- "ystar"
  expect_equal(coef(recm(y, target, df)), coef(by_symbol))
})

test_that("ts, matrix and data.frame inputs give identical fits", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 3)
  as_ts <- stats::ts(as.matrix(df), start = c(1980, 1), frequency = 4)
  base <- coef(recm(y, ystar, df))
  expect_equal(coef(recm(y, ystar, as_ts)), base)
  expect_equal(coef(recm(y, ystar, as.matrix(df))), base)

  # A single non-numeric column is taken as the time index, not a regressor.
  dated <- cbind(df, period = as.Date("1980-01-01") + seq_len(nrow(df)) * 90)
  fit_dated <- recm(y, ystar, dated)
  expect_equal(coef(fit_dated), base)
  expect_length(fit_dated$variables$w, 0L)

  # More than one non-numeric column is ambiguous and is refused.
  worse <- cbind(dated, label = rep("a", nrow(df)))
  expect_error(recm(y, ystar, worse), "more than one non-numeric")
})

test_that("a logical column is a dummy regressor, not the time index", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 24)
  n <- nrow(df)
  on <- seq_len(n) %in% 80:140

  # Logical is how R writes a dummy. Left to the non-numeric rule it would be
  # taken as the time index and never enter the equation at all: no error, no
  # warning, and a regressor quietly missing.
  as_logical <- df
  as_logical$recession <- on
  fit <- recm(y, ystar, as_logical)
  expect_identical(fit$variables$w, "recession")
  expect_identical(fit$variables$w_terms, "recession")
  expect_identical(fit$variables$w_diff, FALSE)

  # And it is the same fit as writing the same dummy out as 0 and 1.
  as_numeric <- df
  as_numeric$recession <- as.numeric(on)
  expect_equal(coef(fit), coef(recm(y, ystar, as_numeric)))

  # The conversion also frees the one non-numeric slot for a real index, which
  # a logical column used to occupy.
  with_date <- as_logical
  with_date$period <- as.Date("1980-01-01") + seq_len(n) * 90
  fit_dated <- recm(y, ystar, with_date)
  expect_identical(fit_dated$variables$w, "recession")
  expect_s3_class(fit_dated$index, "Date")
})

test_that("remaining numeric columns enter as exogenous regressors", {
  df <- simulate_recm(n = 800L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      extra = 0.5, seed = 5)
  expect_named(df, c("y", "ystar", "shock"))
  fit <- recm(y, ystar, df)
  expect_identical(fit$variables$w, "shock")
  # tr_exog defaults to TRUE, so the column enters differenced and renamed.
  expect_true(fit$tr_exog)
  expect_identical(fit$variables$w_terms, "d_shock")
  expect_equal(unname(coef(fit)[["d_shock"]]), 0.5, tolerance = 0.05)
  # The design column really is the difference of the supplied one. The kept
  # rows are the trailing ones, every dropped row being a leading NA.
  expect_equal(unname(fit$model$x[, "d_shock"]),
               utils::tail(diff(df$shock), fit$nobs), tolerance = 1e-12)
})

test_that("tr_exog = FALSE leaves the exogenous regressors in levels", {
  df <- simulate_recm(n = 800L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      extra = 0.5, tr_exog = FALSE, seed = 5)
  fit <- recm(y, ystar, df, tr_exog = FALSE)
  expect_false(fit$tr_exog)
  expect_identical(fit$variables$w_terms, "shock")
  expect_equal(unname(coef(fit)[["shock"]]), 0.5, tolerance = 0.05)

  # Differencing a series the equation wants in levels does not recover it.
  wrong <- recm(y, ystar, df)
  expect_gt(abs(unname(coef(wrong)[["d_shock"]]) - 0.5), 0.1)

  expect_error(recm(y, ystar, df, tr_exog = NA), "must be TRUE or FALSE")
  expect_error(recm(y, ystar, df, tr_exog = c(TRUE, TRUE)),
               "must be TRUE or FALSE")
})

test_that("the dummy test wants both values, not merely no others", {
  expect_true(recm:::is_dummy_column(c(0, 1, 1, 0, 1)))
  expect_true(recm:::is_dummy_column(c(0, 1, NA, 1, 0)))
  # A constant is 0/1 valued and is still not a dummy. Left in levels it would
  # be an intercept, which this equation deliberately does not have.
  expect_false(recm:::is_dummy_column(rep(1, 10)))
  expect_false(recm:::is_dummy_column(rep(0, 10)))
  expect_false(recm:::is_dummy_column(c(0, 1, 2)))
  expect_false(recm:::is_dummy_column(c(0, 0.5, 1)))
  expect_false(recm:::is_dummy_column(c(-1, 0, 1)))
  expect_false(recm:::is_dummy_column(c(NA_real_, NA_real_)))
})

test_that("a dummy is entered in levels beside a differenced regressor", {
  df <- simulate_recm(n = 400L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 20)
  n <- nrow(df)
  set.seed(21)
  df$recession <- as.numeric(seq_len(n) %in% 120:160)
  df$oil <- cumsum(stats::rnorm(n))

  fit <- recm(y, ystar, df)
  expect_true(fit$tr_exog)
  expect_identical(fit$variables$w, c("recession", "oil"))
  # The transform is chosen per column, so the two coexist untransformed on
  # one another's terms, and the names say which is which.
  expect_identical(fit$variables$w_diff, c(FALSE, TRUE))
  expect_identical(fit$variables$w_terms, c("recession", "d_oil"))

  # The design column really is the dummy itself and not its difference. The
  # kept rows are the trailing ones, every dropped row being a leading NA.
  expect_equal(unname(fit$model$x[, "recession"]),
               utils::tail(df$recession, fit$nobs))
  expect_equal(unname(fit$model$x[, "d_oil"]),
               utils::tail(diff(df$oil), fit$nobs), tolerance = 1e-12)

  expect_output(print(fit), "recession \\(level\\), oil \\(first difference\\)")
})

test_that("tr_exog = FALSE puts everything in levels, dummy or not", {
  df <- simulate_recm(n = 400L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 22)
  n <- nrow(df)
  set.seed(23)
  df$recession <- as.numeric(seq_len(n) %in% 100:140)
  df$spread <- stats::rnorm(n)

  fit <- recm(y, ystar, df, tr_exog = FALSE)
  expect_identical(fit$variables$w_diff, c(FALSE, FALSE))
  expect_identical(fit$variables$w_terms, c("recession", "spread"))
})

test_that("a design name collision is named rather than silently resolved", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 16)
  # `ec` is the error correction term's own name, so in levels it collides.
  clash <- cbind(df, ec = stats::rnorm(nrow(df)))
  expect_error(recm(y, ystar, clash, tr_exog = FALSE), "duplicate column names")
  expect_error(recm(y, ystar, cbind(df, dy_lag1 = stats::rnorm(nrow(df))),
                    m = 2, tr_exog = FALSE), "duplicate column names")

  # Differencing renames a regressor out of the way, so a differenced one
  # cannot collide: the `d_` prefix is injective and no reserved name carries
  # it. A dummy keeps its own name, so it can collide where that one could not.
  fit <- recm(y, ystar, clash)
  expect_identical(fit$variables$w_terms, "d_ec")
  dummy_clash <- cbind(df, ec = as.numeric(seq_len(nrow(df)) > 100))
  expect_error(recm(y, ystar, dummy_clash), "duplicate column names")
})

test_that("a fit with no exogenous regressors carries no phantom term", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 17)
  fit <- recm(y, ystar, df)
  # paste0("d_", character(0)) is "d_", not character(0).
  expect_length(fit$variables$w_terms, 0L)
  expect_identical(colnames(fit$model$x), "ec")
})

test_that("the fitted equation reproduces its own residuals", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 6)
  fit <- recm(y, ystar, df)
  implied <- fit$model$dy - as.numeric(fit$model$x %*% coef(fit)) -
    as.numeric(fit$forward_term)
  expect_equal(unname(implied), unname(residuals(fit)), tolerance = 1e-12)
  expect_equal(unname(fitted(fit) + residuals(fit)), unname(fit$model$dy),
               tolerance = 1e-12)
  # The coefficient on the forward term is 1 by construction, never estimated.
  expect_false("forward_loading" %in% names(coef(fit)))
})

test_that("the growth neutrality restriction is imposed when it binds", {
  df <- simulate_recm(n = 800L, a = c(0.35, 0.2), beta = 1, ar = 0.5,
                      const = 0.2, seed = 10)
  fit <- recm(y, ystar, df, m = 2, discount = 0.95)
  expect_true(fit$growth$restricted)
  expect_lt(abs(fit$growth$gap), 1e-6)
  expect_equal(fit$growth$implied_level_gap, 0, tolerance = 1e-6)

  slack <- recm(y, ystar, df, m = 2, discount = 1)
  expect_false(slack$growth$restricted)
  expect_lt(abs(slack$growth$gap), 1e-10)
})

test_that("standard errors differ from the final sweep's OLS ones", {
  # The just-identified GMM sandwich includes dZ/da', which the last
  # iteration's lm() standard errors omit.
  df <- simulate_recm(n = 600L, a = 0.3, beta = 1, ar = 0.6, const = 0.2,
                      seed = 12)
  fit <- recm(y, ystar, df)
  naive <- sqrt(fit$sigma2 * diag(solve(crossprod(fit$model$x))))
  expect_true(all(is.finite(fit$std.error)))
  expect_gt(abs(fit$std.error[["ec"]] - naive[1L]) / naive[1L], 0.05)
})

test_that("unimplemented estimators and bad arguments are refused", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 13)
  expect_error(recm(y, ystar, df, method = "nls"), "not implemented")
  expect_error(recm(y, ystar, df, method = "gmm"), "not implemented")
  expect_error(recm(y, ystar, df, m = 1, discount = 0.95),
               "instantaneous adjustment")
  expect_error(recm(y, ystar, df, m = 0), "positive whole number")
  expect_error(recm(y, ystar, df, discount = 1.2), "in \\(0, 1\\]")
  expect_error(recm(nope, ystar, df), "does not name a column")
  expect_error(recm(y, y, df), "same column")
  expect_error(recm(y, ystar, list(y = 1, ystar = 2)), "must be a ts")
  df$y[10] <- NA
  expect_error(recm(y, ystar, df), "must be complete")
})

test_that("print and summary run without error", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 14)
  fit <- recm(y, ystar, df)
  expect_output(print(fit), "Rational error correction model")
  expect_output(print(summary(fit)), "Growth neutrality")
  expect_output(print(summary(fit)), "Auxiliary model")
  expect_output(print(summary(fit)), "Forward loading")
  expect_output(print(summary(fit)), "R-squared \\(uncentered\\)")
})

test_that("the reported R-squared is the uncentered one", {
  df <- simulate_recm(n = 400L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 21)
  fit <- recm(y, ystar, df)
  dy <- fit$model$dy

  # Fitted plus residual reconstructs the differenced sample, and the ratio is
  # taken about zero rather than about mean(dy) because no intercept is fitted.
  expect_equal(unname(fitted(fit) + residuals(fit)), dy, tolerance = 1e-12)
  expect_equal(fit$r.squared, 1 - fit$ssr / sum(dy^2), tolerance = 1e-12)
  expect_equal(fit$adj.r.squared,
               1 - (1 - fit$r.squared) * fit$nobs / fit$df.residual,
               tolerance = 1e-12)

  # The drift is left in, so the uncentered measure sits above the centered
  # one it is deliberately not.
  centered <- 1 - fit$ssr / sum((dy - mean(dy))^2)
  expect_gt(fit$r.squared, centered)
  expect_gt(fit$r.squared, 0.9)
})

test_that("the forward loading is reported both ways, without a test", {
  df <- simulate_recm(n = 600L, a = c(0.35, 0.2), beta = 1, ar = 0.6,
                      const = 0.2, seed = 15)
  fl <- recm(y, ystar, df, m = 2)$forward_loading
  expect_named(fl, c("restricted", "restricted_se", "free", "free_se"))
  expect_true(all(vapply(fl, is.finite, logical(1))))
  expect_gt(fl$restricted_se, 0)
  expect_gt(fl$free_se, 0)
  # No statistic: the two standard errors are not comparable.
  expect_false("statistic" %in% names(fl))
  # Both land near the truth, sum(d_i) = 1 - a1 = 0.8 at beta = 1.
  expect_equal(fl$restricted, 0.8, tolerance = 0.25)
  expect_equal(fl$free, 0.8, tolerance = 0.25)
})
