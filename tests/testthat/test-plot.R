## Plotting is exercised on a null pdf device, so nothing reaches the disk and
## nothing depends on a display being available.

# Draw on a throwaway device and return whatever the expression returned.
on_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# The columns are named as strings rather than as bare symbols: inside a
# function body the non-standard evaluation form reads as an unbound variable,
# and recm() accepts either.
fitted_example <- function(...) {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 31)
  recm("y", "ystar", df, ...)
}

test_that("every panel draws, singly and together", {
  fit <- fitted_example()
  on_null_device({
    for (k in 1:9) {
      expect_silent(plot(fit, which = k, ask = FALSE))
    }
    expect_silent(plot(fit, which = 1:9, ask = FALSE))
  })
})

test_that("plot returns its argument invisibly", {
  fit <- fitted_example()
  on_null_device({
    expect_invisible(plot(fit, which = 1, ask = FALSE))
    expect_identical(plot(fit, which = 1, ask = FALSE), fit)
  })
})

test_that("the panels survive an autoregressive term and the mce branch", {
  df2 <- simulate_recm(n = 400L, a = c(0.35, 0.25), beta = 1, ar = 0.5,
                       const = 0.2, seed = 32)
  fit2 <- recm(y, ystar, df2, m = 2)
  mce <- recm(y, ystar,
              simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.6,
                            const = 0.2, expectations = "mce", sd_eq = 0.02,
                            seed = 33),
              expectations = "mce")
  on_null_device({
    expect_silent(plot(fit2, which = 1:9, ask = FALSE))
    expect_silent(plot(mce, which = 1:9, ask = FALSE))
  })
})

test_that("a Date index and an exogenous regressor both plot", {
  df <- simulate_recm(n = 300L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      extra = 0.5, seed = 34)
  df$period <- as.Date("1980-01-01") + seq_len(nrow(df)) * 90
  fit <- recm(y, ystar, df)
  expect_s3_class(fit$index, "Date")
  on_null_device(expect_silent(plot(fit, which = c(5L, 8L), ask = FALSE)))
})

test_that("a character index is plotted against position with labels", {
  df <- simulate_recm(n = 120L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                      seed = 35)
  df$period <- paste0("q", seq_len(nrow(df)))
  fit <- recm(y, ystar, df)
  expect_type(fit$index, "character")

  # recm_time() cannot plot labels on a numeric axis, so it falls back to
  # position and hands the labels to the axis instead.
  tm <- recm:::recm_time(fit$index, fit$nobs)
  expect_equal(tm$x, seq_len(fit$nobs))
  expect_identical(tm$labels, as.character(fit$index))

  on_null_device(expect_silent(plot(fit, which = c(5L, 8L), ask = FALSE)))
})

test_that("recm_time passes a numeric or Date index straight through", {
  quarters <- seq(1980, by = 0.25, length.out = 10L)
  expect_null(recm:::recm_time(quarters, 10L)$labels)
  expect_null(recm:::recm_time(as.Date("1980-01-01") + 1:10, 10L)$labels)
  # A length mismatch cannot be plotted against, so position is used.
  tm <- recm:::recm_time(letters[1:3], 10L)
  expect_equal(tm$x, seq_len(10L))
  expect_null(tm$labels)
})

test_that("horizon and lag_max are honoured", {
  fit <- fitted_example()
  w <- recm:::recm_weights(fit, horizon = 12L)
  expect_length(w$d, 13L)
  expect_equal(w$horizon, 12L)
  # The default trims the negligible tail rather than returning the cap.
  auto <- recm:::recm_weights(fit)
  expect_lt(auto$horizon, 400L)
  expect_gte(auto$horizon, 8L)
  # Whatever it keeps, essentially all of the loading is inside it.
  expect_equal(sum(auto$d), fit$d_sum, tolerance = 1e-3)

  on_null_device({
    expect_silent(plot(fit, which = c(3L, 4L, 9L), horizon = 15L, ask = FALSE))
    expect_silent(plot(fit, which = 7L, lag_max = 5L, ask = FALSE))
  })
})

test_that("recm_trim floors a fast decay and keeps a slow one", {
  # Everything past the first element is negligible, so the floor binds.
  expect_equal(recm:::recm_trim(c(1, 1e-9, 1e-12)), 8L)
  # A sequence that is still large at the end is kept whole; element i + 1 is
  # horizon i, so 30 elements are horizons 0 to 29.
  expect_equal(recm:::recm_trim(rep(1, 30L)), 29L)
  # A degenerate all-zero sequence has no scale to compare against.
  expect_equal(recm:::recm_trim(numeric(5L)), 8L)
})

test_that("bad arguments are refused with a message naming the range", {
  fit <- fitted_example()
  expect_error(plot(fit, which = 0), "whole numbers in 1:9")
  expect_error(plot(fit, which = 10), "whole numbers in 1:9")
  expect_error(plot(fit, which = 2.5), "whole numbers in 1:9")
  expect_error(plot(fit, which = integer(0)), "whole numbers in 1:9")
  expect_error(plot(fit, which = NA), "whole numbers in 1:9")
  expect_error(plot(fit, which = 1, horizon = 0), "positive whole number")
  expect_error(plot(fit, which = 1, horizon = c(4, 5)), "positive whole number")
  expect_error(plot(fit, which = 1, lag_max = -2), "positive whole number")
  expect_error(plot(fit, which = 1, id_n = -1), "non-negative whole number")
})

test_that("id_n = 0 labels nothing and main overrides the caption", {
  fit <- fitted_example()
  on_null_device({
    expect_silent(plot(fit, which = c(1L, 5L, 6L), id_n = 0L, ask = FALSE))
    expect_silent(plot(fit, which = 1:2, main = "one title", ask = FALSE))
  })
})

test_that("gap_path reproduces the half-life it was factored out of", {
  # At m = 1 the gap decays at lambda per period in closed form.
  lambda <- 0.7
  path <- recm:::gap_path(1 - lambda, 20L)
  expect_equal(path[1L], 1)
  expect_equal(path, lambda^(0:20), tolerance = 1e-12)
  expect_equal(recm:::half_life(1 - lambda), log(0.5) / log(lambda),
               tolerance = 1e-8)

  # And with an autoregressive term the recursion still starts at a unit gap
  # and closes it.
  path2 <- recm:::gap_path(c(0.35, 0.25), 200L)
  expect_equal(path2[1L], 1)
  expect_lt(abs(path2[201L]), 1e-6)

  # A model that never halves the gap has no half-life rather than a wrong one.
  expect_true(is.na(recm:::half_life(1e-9, max_h = 10L)))
})
