## predict.recm() iterates the decision rule forward rather than predicting one
## step ahead from realised data, so there are two separate things to check:
## that the arithmetic of the decomposition is exact, and that the simulated
## path obeys the properties the estimated model is supposed to have. The
## balanced growth test below is the second kind, and is the one that would
## catch a forward term wired to the wrong state.

forecast_example <- function(n = 400L, seed = 41L, ...) {
  df <- simulate_recm(n = n, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = seed, ...)
  df
}

test_that("the contributions add up, in differences and in levels", {
  df <- forecast_example()
  fit <- recm("y", "ystar", df)
  fc <- predict(fit, n_ahead = 20L)

  expect_s3_class(fc, "recm_forecast")
  expect_identical(fc$horizon, 20L)
  expect_equal(
    unname(rowSums(fc$contributions[, fc$terms, drop = FALSE])),
    unname(fc$fit)
  )
  expect_equal(
    unname(rowSums(fc$levels[, c(fc$base_name, fc$terms), drop = FALSE])),
    unname(fc$level)
  )
  # The level is the forecast origin plus the accumulated differences.
  expect_equal(unname(fc$level), fc$base + cumsum(unname(fc$fit)))
  expect_equal(fc$base, df$y[nrow(df)])
})

test_that("the first forecast step is the estimated equation by hand", {
  df <- forecast_example(seed = 42L)
  fit <- recm("y", "ystar", df)
  fc <- predict(fit, n_ahead = 3L)
  n <- nrow(df)

  # m = 1 and no exogenous regressors, so dy = a0 * (ystar[T] - y[T]) + Z.
  gap <- df$ystar[n] - df$y[n]
  expect_equal(unname(fc$contributions$ec[1L]),
               unname(coef(fit)[["ec"]]) * gap)
  expect_equal(unname(fc$fit[1L]),
               unname(coef(fit)[["ec"]]) * gap + unname(fc$forward_term[1L]))
  expect_equal(unname(fc$level[1L]), df$y[n] + unname(fc$fit[1L]))
})

test_that("a balanced growth target closes the gap", {
  df <- forecast_example(n = 600L, seed = 43L)
  fit <- recm("y", "ystar", df)
  # At the auxiliary model's own mean the state is a fixed point, so every
  # forecast of dystar equals g and the forward term settles at sum(d_i) * g,
  # which is g at beta = 1. The gap then decays at rate 1 - a0 and the
  # equation is growth neutral: y catches ystar rather than trailing it.
  g <- fit$aux$mean
  n <- nrow(df)
  h <- 300L
  ystar_future <- df$ystar[n] + g * seq_len(h)

  fc <- predict(fit, newdata = data.frame(ystar = ystar_future))
  expect_equal(unname(fc$fit[h]), g, tolerance = 1e-6)
  expect_equal(unname(fc$level[h]), ystar_future[h], tolerance = 1e-6)
  # Getting there monotonically from a sample that ends near the path.
  expect_true(abs(fc$level[h] - ystar_future[h]) <
                abs(fc$level[1L] - ystar_future[1L]))
})

test_that("the autoregressive term is seeded with the last observed change", {
  df <- simulate_recm(n = 800L, a = c(0.35, 0.25), beta = 1, ar = 0.5,
                      const = 0.2, seed = 46L)
  fit <- recm("y", "ystar", df, m = 2)
  fc <- predict(fit, n_ahead = 4L)
  n <- nrow(df)
  b <- unname(coef(fit)[["dy_lag1"]])

  expect_equal(unname(fc$contributions$dy_lag1[1L]),
               b * (df$y[n] - df$y[n - 1L]))
  # The second row uses the model's own first forecast difference, which is
  # what makes this a simulation rather than a one step ahead prediction.
  expect_equal(unname(fc$contributions$dy_lag1[2L]), b * unname(fc$fit[1L]))
})

test_that("the forward term picks up exactly where estimation left it", {
  # The sharpest check available on the state the forward term is built from:
  # hold the fit entirely fixed and move only the point at which its sample is
  # cut. The coefficients, the auxiliary autoregression and the realised path
  # of the target are then identical on both sides, so the Z that predict()
  # forms for the removed rows must equal the Z the estimator formed for those
  # same rows. An off-by-one in the state would shift it by a period and show
  # up here, where the balanced growth test cannot see it: at a fixed point
  # every state is the same and a misalignment is invisible.
  for (exp_kind in c("var", "mce")) {
    df <- simulate_recm(n = 400L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                        expectations = exp_kind, sd_eq = 0.02, seed = 60L)
    fit <- recm("y", "ystar", df, expectations = exp_kind)
    n <- nrow(df)
    h <- 5L
    keep <- n - h
    rows <- (keep + 1L):n

    cut <- fit
    cut$model$data <- fit$model$data[seq_len(keep), , drop = FALSE]
    cut$model$index <- fit$model$index[seq_len(keep)]

    fc <- predict(cut, newdata = df[rows, "ystar", drop = FALSE])
    expect_equal(
      unname(fc$forward_term),
      unname(fit$forward_term[as.character(rows)]),
      info = exp_kind
    )
    # And with the forward term right, the first forecast difference is the
    # estimator's own fitted value for that row.
    expect_equal(unname(fc$fit[1L]),
                 unname(fitted(fit)[as.character(rows[1L])]),
                 info = exp_kind)
  }
})

test_that("n_ahead projects the target from the auxiliary autoregression", {
  df <- forecast_example(seed = 45L)
  fit <- recm("y", "ystar", df)
  n <- nrow(df)
  projected <- recm:::aux_forecast_path(
    c(NA_real_, diff(df$ystar)), fit$aux, 10L
  )
  by_hand <- predict(
    fit, newdata = data.frame(ystar = df$ystar[n] + cumsum(projected))
  )
  auto <- predict(fit, n_ahead = 10L)

  expect_true(auto$projected)
  expect_false(by_hand$projected)
  expect_equal(unname(auto$fit), unname(by_hand$fit))
  expect_equal(unname(auto$level), unname(by_hand$level))
})

test_that("ts, matrix and data.frame newdata agree, and a y column is unused", {
  set.seed(52)
  df <- forecast_example(seed = 44L)
  fit <- recm("y", "ystar", df)
  n <- nrow(df)
  path <- df$ystar[n] + cumsum(rep(0.4, 6))
  base_fc <- predict(fit, newdata = data.frame(ystar = path))

  as_matrix <- as.matrix(data.frame(ystar = path))
  expect_equal(predict(fit, newdata = as_matrix)$fit, base_fc$fit)
  as_ts <- stats::ts(as_matrix, start = c(2001, 1), frequency = 4)
  expect_equal(unname(predict(fit, newdata = as_ts)$fit), unname(base_fc$fit))
  # A ts carries its own index, so the forecast is labelled by it.
  expect_equal(unname(predict(fit, newdata = as_ts)$time[1L]), 2001)

  # The decision variable is produced, not read: a y column of nonsense in
  # newdata must leave the forecast untouched.
  with_y <- data.frame(ystar = path, y = rnorm(6L, 1e4))
  expect_equal(unname(predict(fit, newdata = with_y)$fit), unname(base_fc$fit))
})

test_that("n_ahead truncates newdata", {
  df <- forecast_example(seed = 44L)
  fit <- recm("y", "ystar", df)
  n <- nrow(df)
  nd <- data.frame(ystar = df$ystar[n] + cumsum(rep(0.4, 6)))

  full <- predict(fit, newdata = nd)
  short <- predict(fit, newdata = nd, n_ahead = 3L)
  expect_identical(short$horizon, 3L)
  expect_equal(short$fit, full$fit[1:3])
})

test_that("the first new difference of an exogenous regressor spans the join", {
  df <- forecast_example(seed = 44L, extra = 0.5)
  fit <- recm("y", "ystar", df)
  expect_identical(fit$variables$w, "shock")
  expect_identical(fit$variables$w_terms, "d_shock")

  n <- nrow(df)
  nd <- data.frame(ystar = df$ystar[n] + cumsum(rep(0.3, 5)),
                   shock = df$shock[n] + seq_len(5L))
  fc <- predict(fit, newdata = nd)
  # shock rises by one per period including across the join, so every
  # differenced contribution is the coefficient itself.
  expect_equal(unname(fc$contributions$d_shock),
               rep(unname(coef(fit)[["d_shock"]]), 5L))
})

test_that("a level regressor enters the forecast in levels", {
  df <- forecast_example(seed = 50L, extra = 0.5, tr_exog = FALSE)
  fit <- recm("y", "ystar", df, tr_exog = FALSE)
  expect_identical(fit$variables$w_terms, "shock")

  n <- nrow(df)
  nd <- data.frame(ystar = df$ystar[n] + cumsum(rep(0.3, 4)),
                   shock = c(0.5, -0.5, 1, 0))
  fc <- predict(fit, newdata = nd)
  expect_equal(unname(fc$contributions$shock),
               unname(coef(fit)[["shock"]]) * nd$shock)
})

test_that("the mce branch forecasts and its contributions add up", {
  df <- simulate_recm(n = 500L, a = 0.3, beta = 1, ar = 0.6, const = 0.2,
                      expectations = "mce", sd_eq = 0.02, seed = 47L)
  fit <- recm("y", "ystar", df, expectations = "mce")
  fc <- predict(fit, n_ahead = 15L)

  expect_identical(fc$expectations, "mce")
  expect_true(all(is.finite(fc$fit)))
  expect_equal(
    unname(rowSums(fc$contributions[, fc$terms, drop = FALSE])),
    unname(fc$fit)
  )
})

test_that("predict refuses what it cannot forecast", {
  df <- forecast_example(n = 300L, seed = 48L)
  fit <- recm("y", "ystar", df)

  expect_error(predict(fit), "nothing to forecast over")
  expect_error(predict(fit, n_ahead = 0), "positive whole number")
  expect_error(predict(fit, n_ahead = c(2, 3)), "positive whole number")
  expect_error(predict(fit, n_ahead = 2.5), "positive whole number")
  expect_error(predict(fit, newdata = data.frame(z = 1:5)),
               "has no column `ystar`")
  expect_error(predict(fit, newdata = data.frame(ystar = c(1, NA, 3))),
               "must be complete")
  expect_error(predict(fit, newdata = data.frame(ystar = 1:3), n_ahead = 9L),
               "only 3 rows")

  dfw <- forecast_example(n = 300L, seed = 49L, extra = 0.5)
  fitw <- recm("y", "ystar", dfw)
  expect_error(predict(fitw, n_ahead = 5L), "cannot project")
  expect_error(predict(fitw, newdata = data.frame(ystar = 1:5)),
               "missing the exogenous regressor")
})

test_that("the forecast plot draws and returns its argument invisibly", {
  set.seed(51)
  df <- forecast_example(seed = 51L, extra = 0.5)
  fit <- recm("y", "ystar", df)
  n <- nrow(df)
  fc <- predict(fit, newdata = data.frame(
    ystar = df$ystar[n] + cumsum(rep(0.3, 20L)),
    shock = df$shock[n] + cumsum(rnorm(20L))
  ))

  on_null_device({
    expect_silent(plot(fc))
    expect_silent(plot(fc, type = "level"))
    expect_silent(plot(fc, legend = FALSE, col = c("grey70", "grey40")))
    expect_silent(plot(fc, main = "custom"))
    expect_invisible(plot(fc))
    expect_identical(plot(fc), fc)
  })
  expect_error(plot(fc, legend = NA), "must be TRUE or FALSE")
  expect_error(plot(fc, type = "levels"), "should be one of")
})

test_that("print reports the horizon and the contributions", {
  df <- forecast_example(n = 300L, seed = 53L)
  fit <- recm("y", "ystar", df)
  fc <- predict(fit, n_ahead = 5L)

  expect_output(print(fc), "Forecast from a rational error correction model")
  expect_output(print(fc), "horizon: 5")
  expect_output(print(fc), "projected from the auxiliary autoregression")
  # Nested so the printed frame does not leak into the test log.
  expect_output(expect_invisible(print(fc)))
})
