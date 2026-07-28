#' Estimate a rational error correction model
#'
#' Fits the rational error correction (REC) representation of polynomial
#' adjustment cost behaviour,
#'
#' \deqn{\Delta y_t = a_0 (y^*_{t-1} - y_{t-1})
#'   + \sum_{i=1}^{m-1} a_i \Delta y_{t-i} + Z_t + \delta' W_t + \epsilon_t}
#'
#' where \eqn{Z_t = \sum_{i \ge 0} d_i E[\Delta y^*_{t+i}]} is the forward
#' term. The lead weights \eqn{d_i}, the speed of adjustment \eqn{a_0} and the
#' autoregressive coefficients \eqn{a_i} are not free of one another: all are
#' functions of the same lag polynomial \eqn{A(L)}, and imposing that
#' dependence is what distinguishes a rational error correction model from an
#' error correction model written down by hand. The coefficient on \eqn{Z_t}
#' is fixed at 1 by construction.
#'
#' @param y The decision variable. A bare column name of `data` (not a
#'   string is required, though a character scalar also works).
#' @param y_star The frictionless target, in levels, given rather than
#'   estimated. Also a bare column name of `data`.
#' @param data A `ts`/`mts`, `data.frame` (tibbles and data tables qualify) or
#'   named numeric matrix holding `y`, `y_star`, and any exogenous regressors.
#'   **Every remaining numeric column is used as an exogenous regressor**
#'   \eqn{W_t}, dated \eqn{t}. At most one non-numeric column is allowed and is
#'   taken as the time index. No intercept is fitted: a free constant is
#'   inconsistent with growth neutrality, and the auxiliary autoregression
#'   carries the drift instead.
#' @param expectations How expectations of the target are formed. `"var"`
#'   (default) uses the auxiliary univariate autoregression with information
#'   dated \eqn{t-1}, giving \eqn{Z_t = h' z_{t-1}} in closed form. `"mce"`
#'   uses model consistent expectations, evaluated on the realised path of the
#'   target by the finite-lead backward recursion.
#' @param method Estimator. Only `"ols"`, the FRB/US iterative ordinary least
#'   squares zig-zag, is implemented in this release; `"nls"` and `"gmm"` are
#'   accepted by the signature and error informatively.
#' @param m Order of the adjustment cost polynomial: the highest difference
#'   penalised in the loss. It fixes the number of lags of \eqn{\Delta y} at
#'   `m - 1`, so `m = 1` is the canonical rational error correction model with
#'   no autoregressive terms. Unrelated to the length of the forward sum, which
#'   is never truncated.
#' @param discount The discount factor \eqn{\beta}, calibrated and never
#'   estimated. Must lie in \eqn{(0, 1]}. At the default `discount = 1` the
#'   equation is growth neutral for any coefficients and any `m`, so the
#'   restriction described below is slack.
#'
#' @section Growth neutrality:
#' On a balanced growth path the equation delivers \eqn{y = y^*} only if
#' \eqn{R(a) = 1 - \sum_i a_i - \sum_i d_i = 0}. FRB/US and Dynare satisfy this
#' by adding \eqn{R(a) g} to the equation. `recm()` instead imposes
#' \eqn{R(a) = 0} as a nonlinear restriction on the estimated coefficients, so
#' nothing is added and neutrality holds by construction. The restriction is
#' linearised at each zig-zag sweep, which is exact at the fixed point.
#'
#' At `discount = 1` the restriction holds identically and is not imposed. At
#' `discount < 1` with `m = 1` it degenerates: \eqn{R(a_0) = 0} forces
#' \eqn{a_0 = 1}, instantaneous adjustment, so that combination is refused.
#'
#' @section Standard errors:
#' The fixed point solves \eqn{X'(\Delta y - X\gamma - Z(a)) = 0}, which is
#' just-identified GMM with the regressors as their own instruments. The
#' reported covariance is the corresponding sandwich, including the numerical
#' \eqn{\partial Z / \partial a'} that the final sweep's OLS standard errors
#' would omit.
#'
#' It still conditions on the estimated auxiliary model, and that is not a
#' small matter. Holding the target path fixed across replications, the
#' reported standard errors matched the sampling standard deviation to within
#' a tenth. Resampling the target process as well, they came in at roughly a
#' third of it. Bootstrap the auxiliary model before believing an interval.
#'
#' @return An object of class `"recm"`. Use [summary()][summary.recm] for the
#'   coefficient table and the diagnostics, and [coef()], [vcov()],
#'   [fitted()] and [residuals()] as usual. Notable list elements are
#'   `lead_weights`, `structural` (the adjustment costs \eqn{k_j}),
#'   `growth` (the neutrality gap and whether it was imposed), `aux` (the
#'   auxiliary autoregression) and `forward_loading` (the free versus
#'   restricted total forward loading).
#'
#' @references
#' Tinsley, P. A. (2002). Rational error correction. *Computational
#' Economics*, 19(2), 197--225.
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
#' summary(fit)
#'
#' @export
recm <- function(y, y_star, data,
                 expectations = c("var", "mce"),
                 method = c("ols", "nls", "gmm"),
                 m = 1,
                 discount = 1) {
  cl <- match.call()
  expectations <- match.arg(expectations)
  method <- match.arg(method)

  if (method != "ols") {
    stop(
      "method = \"", method, "\" is not implemented in this release. ",
      "Only the iterative OLS zig-zag (method = \"ols\") is available; ",
      "\"nls\" and \"gmm\" are planned.",
      call. = FALSE
    )
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
        m < 1 || m != round(m)) {
    stop("`m` must be a single positive whole number.", call. = FALSE)
  }
  m <- as.integer(m)

  if (!is.numeric(discount) || length(discount) != 1L || is.na(discount) ||
        discount <= 0 || discount > 1) {
    stop("`discount` must be a single number in (0, 1].", call. = FALSE)
  }
  beta <- as.numeric(discount)

  md <- as_model_data(data)
  cols <- colnames(md$x)
  env <- parent.frame()
  y_name <- resolve_column(substitute(y), env, cols, "y")
  ystar_name <- resolve_column(substitute(y_star), env, cols, "y_star")
  if (identical(y_name, ystar_name)) {
    stop("`y` and `y_star` resolve to the same column.", call. = FALSE)
  }

  if (anyNA(md$x[, y_name]) || anyNA(md$x[, ystar_name])) {
    stop(
      "`", y_name, "` and `", ystar_name, "` must be complete: the forward ",
      "term and the auxiliary autoregression are built from contiguous ",
      "differences.",
      call. = FALSE
    )
  }

  des <- build_design(md$x, y_name, ystar_name, m)
  aux <- fit_aux_var(des$dystar)
  n <- nrow(md$x)

  # Growth neutrality: identically satisfied at beta = 1, degenerate at m = 1
  # when beta < 1, a genuine restriction otherwise.
  if (beta >= 1 - 1e-10) {
    restrict <- FALSE
  } else if (m == 1L) {
    stop(
      "with m = 1 and discount < 1 the growth neutrality restriction ",
      "R(a0) = 1 - sum(d_i) = 0 is satisfied only by a0 = 1, that is by ",
      "instantaneous adjustment. Use discount = 1, under which the ",
      "restriction holds for any a0, or m >= 2, under which it defines a ",
      "surface rather than a point.",
      call. = FALSE
    )
  } else {
    restrict <- TRUE
  }

  zlag <- aux_state_lagged(des$dystar, aux$order)
  keep <- which(
    stats::complete.cases(des$dy, des$x) &
      stats::complete.cases(zlag) &
      seq_len(n) > 1L
  )
  if (length(keep) <= ncol(des$x) + 1L) {
    stop(
      "only ", length(keep), " usable observations for ",
      ncol(des$x), " coefficients.",
      call. = FALSE
    )
  }

  dy_est <- des$dy[keep]
  x_est <- des$x[keep, , drop = FALSE]

  if (expectations == "var") {
    zlag_est <- zlag[keep, , drop = FALSE]
    forward <- function(alg) {
      z_var(alg, aux$companion, aux$selector, zlag_est)
    }
  } else {
    pad <- 2000L
    path <- c(des$dystar, aux_forecast_path(des$dystar, aux, pad))
    forward <- function(alg) z_mce(alg, path, n)[keep]
  }

  fit <- iterative_ols(
    dy = dy_est, x = x_est, m = m, beta = beta, forward = forward,
    rho_phi = aux$rho, restrict = restrict
  )

  gamma <- stats::setNames(fit$gamma, colnames(x_est))
  alg <- fit$alg
  a_hat <- fit$gamma[seq_len(m)]

  cmat <- NULL
  if (restrict) {
    grad <- gn_grad(a_hat, beta)
    cmat <- matrix(c(grad, numeric(ncol(x_est) - m)), nrow = 1L)
  }
  vc <- iterative_ols_vcov(fit, dy_est, x_est, m, beta, forward, cmat)

  tt <- length(keep)
  df_resid <- tt - ncol(x_est) + as.integer(restrict)
  sigma2 <- fit$ssr / df_resid
  se <- sqrt(pmax(diag(vc$vcov), 0))

  d <- lead_weights(alg, 40L)
  costs <- cost_params(alg$alpha, beta)
  g_drift <- aux$mean

  # Levels 2 versus 3 of the nested tests: the total forward loading implied
  # by the cross-equation restrictions, beside the same loading estimated
  # freely with the weight shape held at a_hat.
  #
  # No test statistic is formed. The two standard errors below omit different
  # things - the restricted one propagates a_hat but not Phi_hat, the free one
  # treats the normalised forward term as fixed data - and their covariance is
  # not available in closed form, so a t-ratio on the difference would be
  # arbitrary. Compare the numbers and their intervals by eye.
  d_sum_grad <- vapply(
    seq_len(ncol(x_est)),
    function(j) {
      if (j > m) {
        return(0)
      }
      up <- a_hat
      dn <- a_hat
      up[j] <- up[j] + 1e-6
      dn[j] <- dn[j] - 1e-6
      (pac_algebra(up, beta)$d_sum - pac_algebra(dn, beta)$d_sum) / 2e-6
    },
    numeric(1)
  )
  d_sum_se <- sqrt(max(0, drop(d_sum_grad %*% vc$vcov %*% d_sum_grad)))

  zn <- fit$z / alg$d_sum
  x_free <- cbind(x_est, forward_loading = zn)
  free_fit <- stats::lm.fit(x_free, dy_est)
  af_free <- unname(free_fit$coefficients[ncol(x_free)])
  free_resid <- free_fit$residuals
  af_se <- sqrt(
    sum(free_resid^2) / (tt - ncol(x_free)) *
      solve(crossprod(x_free))[ncol(x_free), ncol(x_free)]
  )

  structure(
    list(
      call = cl,
      method = method,
      expectations = expectations,
      m = m,
      discount = beta,
      coefficients = gamma,
      vcov = vc$vcov,
      std.error = stats::setNames(se, names(gamma)),
      statistic = stats::setNames(gamma / se, names(gamma)),
      p.value = stats::setNames(
        2 * stats::pt(abs(gamma / se), df_resid, lower.tail = FALSE),
        names(gamma)
      ),
      df.residual = df_resid,
      sigma2 = sigma2,
      ssr = fit$ssr,
      residuals = stats::setNames(fit$residuals, md$index[keep]),
      fitted.values = stats::setNames(dy_est - fit$residuals, md$index[keep]),
      forward_term = stats::setNames(fit$z, md$index[keep]),
      alpha = alg$alpha,
      a_poly = alg$ap,
      roots = if (m >= 1L) polyroot(alg$ap) else complex(0),
      lead_weights = d,
      lead_weights_normalised = d / alg$d_sum,
      d_sum = alg$d_sum,
      mean_lead = alg$lead_sum / alg$d_sum,
      half_life = half_life(a_hat),
      structural = costs,
      growth = list(
        gap = gn_gap_alg(alg),
        restricted = restrict,
        gradient = if (restrict) gn_grad(a_hat, beta) else NULL,
        drift = g_drift,
        implied_level_gap = if (is.na(g_drift)) {
          NA_real_
        } else {
          gn_implied_gap(a_hat, beta, g_drift)
        }
      ),
      forward_loading = list(
        restricted = alg$d_sum, restricted_se = d_sum_se,
        free = af_free, free_se = af_se
      ),
      aux = aux,
      hac_lag = vc$hac_lag,
      iterations = fit$iterations,
      converged = fit$converged,
      criterion = fit$criterion,
      nobs = tt,
      index = md$index[keep],
      variables = list(y = y_name, y_star = ystar_name, w = des$w_names),
      model = list(dy = dy_est, x = x_est)
    ),
    class = "recm"
  )
}
