# The cost-parameterisation extension point (docs/01-architecture.md).
#
# The first test is a CONTRACT test over `.COST`: it loops the registry
# rather than naming "geometric" and "free", so a parameterisation added
# later is checked by it without anyone writing a new test. That is the
# whole point of the seam -- if you add an entry and this file starts
# failing, the entry is inconsistent, not the test.

test_that("every registered cost parameterisation satisfies the contract", {
  for (nm in names(.COST)) {
    spec <- .COST[[nm]]
    expect_named(spec, c("npar", "k", "names", "natural", "start", "grid",
                         "warn"),
                 ignore.order = TRUE, info = nm)
    for (m in 1:5) {
      lbl <- paste0(nm, ", m = ", m)
      np <- spec$npar(m)
      expect_true(is.numeric(np) && length(np) == 1 && np >= 1, info = lbl)

      # Everything on the optimisation scale agrees on the parameter count.
      expect_length(spec$names(m), np)
      expect_length(spec$start(m), np)
      expect_equal(ncol(spec$grid(m)), np, info = lbl)
      expect_length(spec$natural(spec$start(m), m), np)
      expect_true(!is.null(names(spec$natural(spec$start(m), m))), info = lbl)

      # k_0 is normalised to 1 and the weights are admissible at the
      # default start -- an unusable default start would send every fit
      # into the fallback grid.
      kk <- spec$k(spec$start(m), m)
      expect_length(kk, m + 1)
      expect_equal(kk[1], 1, info = lbl)
      expect_true(all(is.finite(kk)) && all(kk > 0), info = lbl)

      # Every fallback start must itself be admissible, or it is not a
      # fallback.
      g <- spec$grid(m)
      for (i in seq_len(nrow(g))) {
        kg <- spec$k(as.numeric(g[i, ]), m)
        expect_true(all(is.finite(kg)) && all(kg > 0),
                    info = paste0(lbl, ", grid row ", i))
      }

      w <- spec$warn(m)
      expect_true(is.null(w) || (is.character(w) && length(w) == 1),
                  info = lbl)
    }
  }
})

test_that("the dispatchers agree with the registry", {
  # .k_from_theta/.theta_names/.theta_natural are the names docs/01 gives;
  # they must stay thin wrappers, not a second implementation.
  for (nm in names(.COST)) {
    th <- .COST[[nm]]$start(3)
    expect_equal(.k_from_theta(th, 3, nm), .COST[[nm]]$k(th, 3))
    expect_equal(.theta_names(3, nm), .COST[[nm]]$names(3))
    expect_equal(.theta_natural(th, 3, nm), .COST[[nm]]$natural(th, 3))
    expect_equal(.theta_npar(3, nm), .COST[[nm]]$npar(3))
  }
})

test_that("an unregistered cost parameterisation errors", {
  # It used to fall through to the "free" branch and return c(1, exp(th)) --
  # a plausible cost vector, silently the wrong model.
  expect_error(.k_from_theta(c(1, 0), 3, "hump"), "unknown cost")
  expect_error(.theta_names(3, "hump"), "unknown cost")
  expect_error(.theta_natural(c(1, 0), 3, "hump"), "unknown cost")
  # The message names the registered set so the fix is obvious.
  expect_error(.cost_spec("hump"), "geometric")
})

test_that("the two known parameterisations produce the documented k", {
  # Pins the algebra itself, independently of the registry plumbing.
  # geometric: k_j = kappa * psi^(j-1), j = 1..m, on top of k_0 = 1.
  kappa <- 20
  psi <- 0.5
  th <- c(log(kappa), log(psi / (1 - psi)))
  expect_equal(.k_from_theta(th, 4, "geometric"),
               c(1, kappa * psi^(0:3)))
  expect_equal(unname(.theta_natural(th, 4, "geometric")), c(kappa, psi))
  # free: each k_j is its own parameter.
  expect_equal(.k_from_theta(log(c(3, 7)), 2, "free"), c(1, 3, 7))
  expect_equal(unname(.theta_natural(log(c(3, 7)), 2, "free")), c(3, 7))
})

test_that("cost = 'free' estimates, and gives m free parameters", {
  # The suite previously tested only that "free" WARNS above m = 4, never
  # that it fits. A parameterisation nothing exercises end to end is not
  # covered by the contract test above.
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       cost = "free", var_lags = 2, restarts = 1,
                       quiet = TRUE)
  expect_equal(fit$cost, "free")
  expect_length(fit$theta, 2L)
  expect_equal(names(fit$par)[1:2], c("log k1", "log k2"))
  expect_equal(fit$k[1], 1)
  expect_true(all(fit$k > 0))
  expect_lt(fit$a[1], 0)                       # INVARIANT 2
})

test_that("a wrong-length `start` errors instead of being read positionally", {
  # geometric takes 2; handing it m = 3 values used to run and return
  # numbers built from start[1:2].
  expect_error(
    recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 3,
                  cost = "geometric", start = c(1, 0, 1), var_lags = 2,
                  restarts = 0, quiet = TRUE),
    "takes 2 parameters"
  )
  expect_error(
    recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 3,
                  cost = "free", start = c(1, 0), var_lags = 2,
                  restarts = 0, quiet = TRUE),
    "takes 3 parameters"
  )
})

test_that("a correct-length `start` is honoured", {
  fit <- recm_estimate("y", "ystar", beta = 0.995, data = dgp_pac, m = 2,
                       start = c(log(15), 0.2), var_lags = 2, restarts = 0,
                       quiet = TRUE)
  expect_length(fit$theta, 2L)
  expect_lt(fit$a[1], 0)
})
