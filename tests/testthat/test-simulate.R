## simulate.recm() runs the same recursion as predict.recm(), which is tested
## in test-predict.R, so what is checked here is what simulate() adds: the
## shocks, and the impulse response formed by differencing two runs.
##
## The response is computed as shocked minus baseline. That is only legitimate
## because the decision rule is linear in its inputs, so the two properties
## worth testing hardest are the ones that follow from linearity and would fail
## loudly without it: the response does not depend on the baseline, and it
## still decomposes exactly term by term.

sim_example <- function(n = 400L, seed = 71L, ...) {
  simulate_recm(n = n, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                seed = seed, ...)
}

test_that("an unshocked simulation is the same recursion as a forecast", {
  df <- sim_example(n = 300L, seed = 61L)
  fit <- recm("y", "ystar", df)
  n <- nrow(df)
  nd <- data.frame(ystar = df$ystar[n] + cumsum(rep(0.4, 7)))

  sim <- simulate(fit, newdata = nd)
  fc <- predict(fit, newdata = nd)
  expect_s3_class(sim, "recm_simulation")
  expect_null(sim$shock)
  expect_equal(unname(sim$fit), unname(fc$fit))
  expect_equal(unname(sim$level), unname(fc$level))
  expect_equal(sim$contributions, fc$contributions)

  # match.call() inside an S3 method names the method, not the generic.
  expect_identical(as.character(sim$call[[1L]]), "simulate.recm")
  expect_identical(as.character(fc$call[[1L]]), "predict.recm")
  expect_equal(simulate(fit, newdata = nd, n_ahead = 3L)$fit, sim$fit[1:3])
})

test_that("simulate refuses the arguments it cannot honour", {
  df <- sim_example(n = 300L, seed = 62L)
  fit <- recm("y", "ystar", df)

  # The path is deterministic, so a second replicate would repeat the first
  # and a seed would change nothing. Both are refused rather than ignored.
  expect_error(simulate(fit, nsim = 2, n_ahead = 4L), "`nsim` must be 1")
  expect_error(simulate(fit, nsim = 0, n_ahead = 4L), "`nsim` must be 1")
  expect_error(simulate(fit, seed = 42, n_ahead = 4L), "no effect")
  expect_silent(simulate(fit, nsim = 1, seed = NULL, n_ahead = 4L))
  expect_error(simulate(fit), "nothing to iterate over")
})

test_that("recm_shock validates and prints its specification", {
  s <- recm_shock("ystar", size = 2, period = 3, type = "transitory")
  expect_s3_class(s, "recm_shock")
  expect_identical(s$period, 3L)
  expect_output(print(s), "ystar \\+2, transitory, from period 3")
  expect_output(print(recm_shock("ystar", size = c(1, 2), type = "path")),
                "path of 2 periods")

  expect_error(recm_shock(c("a", "b")), "single non-empty variable name")
  expect_error(recm_shock("ystar", size = NA_real_), "finite numeric")
  expect_error(recm_shock("ystar", size = c(1, 2)), "single number unless")
  expect_error(recm_shock("ystar", period = 0), "positive whole number")
  expect_error(recm_shock("ystar", type = "step"), "should be one of")
})

test_that("a shock has to name an input of the decision rule", {
  df <- sim_example(n = 300L, seed = 64L)
  fit <- recm("y", "ystar", df)

  # y is produced by the recursion, so a shock to it has nowhere to enter.
  expect_error(simulate(fit, n_ahead = 5L, shock = recm_shock("y")),
               "produced by the recursion")
  expect_error(simulate(fit, n_ahead = 5L, shock = recm_shock("zzz")),
               "not an input of this decision rule")
  expect_error(simulate(fit, n_ahead = 5L, shock = 1), "single named number")
  expect_error(simulate(fit, n_ahead = 5L, shock = "ystar"),
               "single named number")
  # The shorthand is a permanent unit shock.
  expect_equal(
    simulate(fit, n_ahead = 5L, shock = c(ystar = 1))$fit,
    simulate(fit, n_ahead = 5L, shock = recm_shock("ystar"))$fit
  )
  expect_error(simulate(fit, n_ahead = 3L, shock = recm_shock("ystar",
                                                              period = 5L)),
               "nothing would happen")
  expect_error(
    simulate(fit, n_ahead = 3L,
             shock = recm_shock("ystar", size = 1:5, type = "path")),
    "shock path runs to period 5"
  )
})

test_that("the response is shocked minus baseline, and still adds up", {
  df <- sim_example(seed = 72L)
  fit <- recm("y", "ystar", df)
  n <- nrow(df)
  nd <- data.frame(ystar = df$ystar[n] + cumsum(rep(0.4, 25)))

  irf <- simulate(fit, newdata = nd, shock = recm_shock("ystar", size = 1.5))
  by_hand_base <- simulate(fit, newdata = nd)
  by_hand_shock <- simulate(fit, newdata = transform(nd, ystar = ystar + 1.5))

  expect_equal(unname(irf$fit),
               unname(by_hand_shock$fit) - unname(by_hand_base$fit))
  expect_equal(unname(irf$level),
               unname(by_hand_shock$level) - unname(by_hand_base$level))
  # The two runs are carried along, and are ordinary simulations themselves.
  expect_equal(unname(irf$baseline$fit), unname(by_hand_base$fit))
  expect_equal(unname(irf$shocked$fit), unname(by_hand_shock$fit))

  # Linearity again: differencing the runs differences the decomposition, so
  # the response splits exactly as the level does.
  expect_equal(unname(rowSums(irf$contributions[, irf$terms, drop = FALSE])),
               unname(irf$fit))
  expect_equal(unname(irf$level), cumsum(unname(irf$fit)))
  expect_identical(irf$base, 0)
  expect_equal(
    unname(rowSums(irf$levels[, c(irf$base_name, irf$terms), drop = FALSE])),
    unname(irf$level)
  )
})

test_that("the response does not depend on the baseline it is measured from", {
  # The claim the whole two-run construction rests on. Three baselines that
  # could hardly be more different - the model's own projection, a flat target
  # and a steeply rising one - must give the same response to the same shock.
  df <- sim_example(seed = 73L)
  fit <- recm("y", "ystar", df)
  n <- nrow(df)
  h <- 30L
  shock <- recm_shock("ystar", size = 1, period = 2)

  projected <- simulate(fit, n_ahead = h, shock = shock)
  flat <- simulate(fit, shock = shock,
                   newdata = data.frame(ystar = rep(df$ystar[n], h)))
  steep <- simulate(fit, shock = shock,
                    newdata = data.frame(ystar = df$ystar[n] + 5 * seq_len(h)))

  expect_equal(unname(flat$fit), unname(projected$fit))
  expect_equal(unname(steep$fit), unname(projected$fit))
  # The baselines themselves are of course nothing like one another.
  expect_false(isTRUE(all.equal(unname(flat$baseline$fit),
                                unname(steep$baseline$fit))))
})

test_that("a permanent shock to the target is fully passed through", {
  # Growth neutrality with beta = 1: in the long run y takes all of a permanent
  # displacement of the target, so the level response converges on the size of
  # the shock, and the response of the difference dies away.
  df <- sim_example(n = 600L, seed = 74L)
  fit <- recm("y", "ystar", df)
  irf <- simulate(fit, n_ahead = 200L, shock = recm_shock("ystar", size = 1))

  expect_equal(unname(irf$level[200L]), 1, tolerance = 1e-6)
  expect_lt(abs(unname(irf$fit[200L])), 1e-6)
  expect_lt(unname(irf$level[2L]), 1)

  # The approach is not asserted to be monotone, and is not: the displacement
  # is one large change in the target, and while it sits in the auxiliary
  # state the autoregression reads it as news about growth and revises that
  # reading period by period, so the forward term can pull the other way. Once
  # the displacement has left the state window there is nothing left for it to
  # respond to and what remains is pure error correction, decaying at exactly
  # 1 - a0 per period.
  a0 <- unname(coef(fit)[["ec"]])
  rows <- (fit$aux$order + 2L):30L
  expect_true(all(abs(unname(irf$contributions$forward[rows])) < 1e-12))
  expect_equal(unname(irf$fit[rows + 1L]), (1 - a0) * unname(irf$fit[rows]),
               tolerance = 1e-6)
})

test_that("a transitory shock to the target is undone", {
  df <- sim_example(n = 600L, seed = 75L)
  fit <- recm("y", "ystar", df)
  irf <- simulate(fit, n_ahead = 200L,
                  shock = recm_shock("ystar", size = 1, type = "transitory"))

  # The target returns to baseline, so the error correction term pulls y back
  # with it and the level response dies out rather than settling anywhere.
  expect_lt(abs(unname(irf$level[200L])), 1e-6)
  expect_gt(max(abs(unname(irf$level))), 1e-3)
})

test_that("the branches respond on the timetables their information sets set", {
  # Timing is where the two branches genuinely differ, and it is not an
  # implementation detail. It is also the thing an impulse response is best
  # placed to expose, because a shock at a known date makes the information
  # set visible in a way a fitted sample never does.
  k <- 4L
  shock <- recm_shock("ystar", size = 1, period = k)
  horizon <- 12L

  # Under "var" the forward term at t is built from states dated t-1, so the
  # shock is unanticipated: nothing whatever moves until the period after it
  # lands. And once the displaced difference has fallen out of the state
  # window, `order` periods later, the forward term is done responding.
  df_var <- simulate_recm(n = 400L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                          expectations = "var", sd_eq = 0.02, seed = 76L)
  fit_var <- recm("y", "ystar", df_var)
  irf_var <- simulate(fit_var, n_ahead = horizon, shock = shock)
  expect_true(all(abs(unname(irf_var$fit[seq_len(k)])) < 1e-12))
  expect_gt(abs(unname(irf_var$fit[k + 1L])), 1e-6)
  spent <- (k + fit_var$aux$order + 1L):horizon
  expect_true(all(abs(unname(irf_var$contributions$forward[spent])) < 1e-12))

  # Under "mce" the forward term is evaluated on the realised path, so the
  # shock is anticipated from the very first period. What Z responds by is
  # exactly the lead weight the displaced period carries: at t the target is
  # displaced k - t periods ahead, so the response is d_{k-t}, and nothing at
  # all once the displacement is in the past.
  df_mce <- simulate_recm(n = 400L, a = 0.3, beta = 1, ar = 0.5, const = 0.2,
                          expectations = "mce", sd_eq = 0.02, seed = 76L)
  fit_mce <- recm("y", "ystar", df_mce, expectations = "mce")
  irf_mce <- simulate(fit_mce, n_ahead = horizon, shock = shock)
  d <- recm:::lead_weights(
    recm:::pac_algebra(recm:::recm_a(fit_mce), fit_mce$discount), k
  )
  expect_equal(unname(irf_mce$contributions$forward[seq_len(k)]),
               rev(d[seq_len(k)]))
  expect_true(all(abs(unname(
    irf_mce$contributions$forward[(k + 1L):horizon]
  )) < 1e-12))
  # In the first period nothing has happened yet to correct, so the whole of
  # the response is that anticipation.
  expect_equal(unname(irf_mce$fit[1L]),
               unname(irf_mce$contributions$forward[1L]))
  expect_gt(abs(unname(irf_mce$fit[1L])), 1e-6)
})

test_that("a shock path of your own lands where it is put", {
  df <- sim_example(seed = 77L)
  fit <- recm("y", "ystar", df)
  bump <- c(0.5, 1, 0.5)
  irf <- simulate(fit, n_ahead = 20L,
                  shock = recm_shock("ystar", size = bump, period = 4L,
                                     type = "path"))
  by_hand <- simulate(fit, n_ahead = 20L, shock = recm_shock(
    "ystar", size = c(numeric(3L), bump, numeric(14L)), type = "path"
  ))
  expect_equal(unname(irf$fit), unname(by_hand$fit))
})

test_that("an exogenous regressor can be shocked, over an invented baseline", {
  df <- sim_example(n = 400L, seed = 78L, extra = 0.5)
  fit <- recm("y", "ystar", df)
  expect_identical(fit$variables$w, "shock")

  # predict() refuses to project this regressor, and rightly so. A response is
  # a different matter: the baseline cancels, so holding it flat is admissible
  # here even though it would be a fiction as a forecast.
  expect_error(predict(fit, n_ahead = 10L), "cannot project")
  irf <- simulate(fit, n_ahead = 10L,
                  shock = recm_shock("shock", size = 1, period = 3L))

  # tr_exog differences the regressor, so a permanent unit step in its level
  # is a single unit change in the period it lands and nothing afterwards.
  expected <- numeric(10L)
  expected[3L] <- unname(coef(fit)[["d_shock"]])
  expect_equal(unname(irf$contributions$d_shock), expected)
  expect_equal(unname(rowSums(irf$contributions[, irf$terms, drop = FALSE])),
               unname(irf$fit))
  # The target is untouched by a shock to the regressor.
  expect_equal(unname(irf$shocked$forward_term),
               unname(irf$baseline$forward_term))
})

test_that("the simulation plot draws paths and responses", {
  df <- sim_example(n = 300L, seed = 79L)
  fit <- recm("y", "ystar", df)
  sim <- simulate(fit, n_ahead = 15L)
  irf <- simulate(fit, n_ahead = 15L, shock = recm_shock("ystar"))

  on_null_device({
    expect_silent(plot(sim))
    expect_silent(plot(sim, type = "level"))
    expect_silent(plot(irf))
    expect_silent(plot(irf, type = "level"))
    expect_silent(plot(irf, legend = FALSE, main = "custom"))
    expect_identical(plot(irf), irf)
  })
  expect_error(plot(irf, legend = NA), "must be TRUE or FALSE")
})

test_that("print says which of the two things it is showing", {
  df <- sim_example(n = 300L, seed = 80L)
  fit <- recm("y", "ystar", df)

  expect_output(print(simulate(fit, n_ahead = 5L)), "Dynamic simulation")
  irf <- simulate(fit, n_ahead = 5L, shock = recm_shock("ystar", size = 2))
  expect_output(print(irf), "Impulse response")
  expect_output(print(irf), "shock: ystar \\+2, permanent, from period 1")
  expect_output(print(irf), "deviation from the baseline")
  expect_output(expect_invisible(print(irf)))
})
