# The only module that touches `data`. Assembles the design, calls the
# pipeline, runs the optimiser. See docs/01-architecture.md.
#
# Pipeline:
#
#   theta --.k_from_theta()--> k --.lq_alpha()--> alpha
#         --.alpha_to_a()/.scalars()--> a, G, sum_d
#         --zmech$z()--> Z, the forward sum
#         --> residual --> SSR or GMM criterion
#
# The mechanism supplying Z is built once, before the optimiser, and the
# residual function consumes its output without knowing how it was made.
# Under `.zmech_var()` that is Z_t = h' z_{t-1}, the VAR closure.

# Penalty returned to the optimiser for an inadmissible theta. R-2 will
# replace this cliff with a continuous boundary penalty -- Nelder-Mead
# handles the cliff badly.
.PENALTY_INADMISSIBLE <- 1e12

# Above this order, m free cost parameters are meaningless: Var((1-L)^j y)
# grows like choose(2(j-1), j-1), about 3.5e10 at j = 20, so the cost
# weights would have to span ten orders of magnitude to contribute
# comparably. cost = "geometric" is the sane choice above this.
.M_FREE_MAX <- 4L


# ---- cost parameterisations ----
#
# Extension point (docs/01-architecture.md). A parameterisation is one entry
# in .COST plus one word in the `cost` vector of recm_estimate()'s signature.
# Nothing else in the package branches on `cost`.
#
#   "geometric"  k_j = kappa * psi^(j-1), j = 1..m   -> 2 free parameters
#   "free"       each k_j free                       -> m free parameters
#
# k_0 is normalised to 1 throughout; only ratios are identified.
#
# Every entry must supply all seven fields, and they must agree on the
# parameter count: names(), natural(), start() all have length npar(m) and
# grid() has npar(m) columns. test-cost.R asserts exactly that for every
# registered entry, so a new one is checked without writing a new test.
#
#   npar(m)      integer, number of free parameters at order m
#   k(th, m)     theta -> cost vector, length m+1, k[1] == 1
#   names(m)     labels on the OPTIMISATION scale (these reach summary())
#   natural(th, m)  named vector on the interpretable scale
#   start(m)     default starting value
#   grid(m)      npar(m)-column matrix of fallback starts, tried in order
#                when start(m) is inadmissible
#   warn(m)      character(1) warning, or NULL if the order is sensible
#
# This used to be three functions each shaped `if (cost == "geometric") ...
# else ...`, with starting values and the restart grid written the same way
# at their call sites. A name that reached those `else` branches was silently
# treated as "free": .k_from_theta() returned c(1, exp(th)), a perfectly
# plausible cost vector, and estimation ran to completion on it. Dispatch is
# now a lookup that stops on an unknown name.
.cost_spec <- function(cost) {
  spec <- .COST[[cost]]
  if (is.null(spec)) {
    stop("unknown cost parameterisation '", cost, "'. Registered: ",
         paste(names(.COST), collapse = ", "),
         ". A new one is an entry in .COST plus a word in the `cost` ",
         "argument of recm_estimate() -- see docs/01-architecture.md.")
  }
  spec
}

.COST <- list(
  geometric = list(
    npar = function(m) 2L,
    k = function(th, m) {
      kappa <- exp(th[1])
      psi <- 1 / (1 + exp(-th[2]))
      c(1, kappa * psi^(0:(m - 1)))
    },
    names = function(m) c("log kappa", "logit psi"),
    natural = function(th, m) {
      c(kappa = unname(exp(th[1])), psi = unname(1 / (1 + exp(-th[2]))))
    },
    start = function(m) c(log(20), 0),
    grid = function(m) {
      unname(as.matrix(expand.grid(log(c(1, 5, 20, 100)), c(-2, -1, 0, 1))))
    },
    warn = function(m) NULL
  ),
  free = list(
    npar = function(m) as.integer(m),
    k = function(th, m) c(1, exp(th)),
    names = function(m) paste0("log k", seq_len(m)),
    natural = function(th, m) {
      stats::setNames(exp(th), paste0("k", seq_len(m)))
    },
    start = function(m) rep(log(5), m),
    grid = function(m) matrix(rep(log(c(0.5, 2, 10)), m), ncol = m),
    warn = function(m) {
      if (m <= .M_FREE_MAX) return(NULL)
      paste0("cost = \"free\" with m = ", m, " > ", .M_FREE_MAX,
             ": Var((1-L)^j y) grows like choose(2(j-1), j-1), so the k_j ",
             "would have to span ten orders of magnitude to contribute ",
             "comparably. Use cost = \"geometric\".")
    }
  )
)

# ---- estimators ----
#
# Extension point (docs/01-architecture.md). Unlike the cost registry above
# this is a documented contract rather than a lookup table, deliberately:
# the two estimators share the residual function and little else, and GMM
# has to build its instruments BEFORE sample selection -- upstream of
# anywhere a plug-in function could run. See docs/01 for the reasoning.
#
# Adding an estimator is three edits:
#
#   1. the name in recm_estimate()'s `method` argument
#   2. if it needs the automatic instruments, add it to .METHODS_IV
#   3. a branch in the estimator block that sets ALL SEVEN of
#
#        th        theta at the optimum
#        r         resid_full(th); carries $b, $lin, $e
#        par_all   c(th, r$lin), the full parameter vector
#        V         covariance, length(par_all) square
#        fit       optimiser result; only $convergence is read
#        extras    estimator-specific diagnostics, possibly empty list
#        objective closure giving the criterion at any theta
#
# docs/01 used to name only th, V and extras. Setting the other four is not
# optional -- the object assembly at the end of recm_estimate() reads every
# one of them.
#
# Estimators whose branch consumes `Ziv`. This gate used to be written as a
# bare `method == "gmm"` at the instrument block, which made it an invisible
# second edit site: a new estimator reached that block, got Ziv = NULL, and
# died inside apply() with "dim(X) must have a positive length".
.METHODS_IV <- c("gmm")

.k_from_theta <- function(th, m, cost) .cost_spec(cost)$k(th, m)

.theta_names <- function(m, cost) .cost_spec(cost)$names(m)

.theta_natural <- function(th, m, cost) .cost_spec(cost)$natural(th, m)

.theta_npar <- function(m, cost) .cost_spec(cost)$npar(m)


# ---- main entry point ----

#' Estimate a rational error correction model with polynomial adjustment
#' costs
#'
#' Fits the forward-looking error-correction specification used in the
#' Federal Reserve Board's FRB/US model and the Bank of Canada's LENS model.
#' Given a decision variable, a frictionless target and a **calibrated**
#' discount factor, `recm_estimate()` recovers the adjustment cost weights
#' and reports the implied decision rule.
#'
#' @param y character(1), column name of the decision variable, **in levels**
#'   (usually logs). The equation is estimated in first differences.
#' @param ystar character(1), column name of the frictionless target, in the
#'   same units as `y`. Supplied by the user, typically from a cointegrating
#'   regression. **Its estimation error is not propagated** -- see
#'   [recm_boot()] and roadmap item R-1.
#' @param vars formula or character vector of extra regressors, entering
#'   linearly with free coefficients. Concentrated out of the optimiser.
#' @param beta numeric(1) discount factor. **Calibrate this; do not estimate
#'   it.** It is not identified jointly with the cost parameters, entering
#'   only through \eqn{(1-\beta F)(1-L)}. Use `1/(1+r)` with a *real*
#'   quarterly rate: the loss is quadratic in log deviations, which are
#'   unit-free, so the real rate is the relevant one even when `y` and
#'   `ystar` are nominal. The default of `1` is admissible provided
#'   \eqn{\rho(G) < 1}, but it means indifference between the present and the
#'   arbitrarily distant future and it lengthens the effective lead horizon;
#'   it warns unless `quiet = TRUE`.
#' @param data data.frame ordered in time, **with no gaps**. Not checked --
#'   a gap will silently corrupt every lag. Roadmap item R-7.
#' @param m integer, order of the adjustment cost polynomial. The equation
#'   carries `m - 1` lags of `dy`.
#' @param cost `"geometric"` (`k_j = kappa * psi^(j-1)`, two parameters) or
#'   `"free"` (`m` parameters). `"free"` warns above `m = 4`; see the
#'   variance-amplification note in Details.
#' @param expectations formula or character vector of extra VAR variables.
#'   `d(ystar)` is always included and always occupies position 1.
#' @param var_lags integer lag order of the auxiliary VAR. **This also sets
#'   GMM instrument dating** through the default `iv_lag = var_lags + 2`. A
#'   longer VAR buys forecast accuracy and pays in instrument strength: on
#'   178 quarterly observations first-stage `F` was 13.1 at `p = 1`, 6.0 at
#'   `p = 2` and 3.3 at `p = 4`, and dating the instruments correctly costs
#'   more strength again.
#' @param diff_exp logical; difference the `expectations` variables before
#'   they enter the VAR.
#' @param method `"nls"` (default) or `"gmm"` (two-step, Newey-West).
#' @param instruments formula or character vector of *excluded* instruments
#'   for GMM, added to the automatic set. These are taken as supplied -- the
#'   `iv_lag` rule is **not** applied to them, so date them yourself.
#'   Collinear columns are dropped by QR pivot.
#' @param iv_lag integer, the lag at which those *automatic* GMM instruments
#'   that are built from `y` or `ystar` are dated -- the error-correction
#'   term, the lags of `dy`, and the `d(ystar)` states. `NULL` (default)
#'   uses `var_lags + 2`, the smallest lag the dating rule permits. Other
#'   state variables are exogenous with respect to measurement error in
#'   `ystar` and stay at `t-1` whatever `iv_lag` is set to.
#'
#'   The composite error carries measurement error `nu` at `t-1` through
#'   `t - var_lags - 1`: the error-correction term contributes `nu_{t-1}`,
#'   and the forward term contributes `d(ystar)` at `t-1 .. t-var_lags`.
#'   Anything built from `y` or `ystar` is therefore only predetermined at
#'   `t - (var_lags + 2)` or earlier.
#'
#'   `iv_lag = 1` puts the equation's own regressors back in the instrument
#'   set. That is more efficient, and correct **only if `ystar` is measured
#'   without error** -- which it usually is not, since it normally comes
#'   from a first-stage cointegrating regression. It warns. See the
#'   Instrument dating section for the size of the trade.
#' @param free_forward logical. If `TRUE`, `a_f` is estimated freely instead
#'   of being restricted to `sum(d_i)`. Comparing the free estimate to
#'   `sum(d_i)` evaluated at the estimated `alpha` is the direct test of the
#'   Euler restriction, and is the single most informative diagnostic the
#'   package offers.
#' @param growth character(1), column name of trend growth. Adds the
#'   growth-neutrality correction `(1 - sum a_i - sum d_i) * g_t` to the
#'   forward term. See Details -- this is a deliberate departure from the
#'   exact optimum.
#' @param hac_lags integer Newey-West bandwidth; `NULL` uses
#'   `floor(4 * (T/100)^(2/9))`.
#' @param start numeric vector of starting values on the optimisation scale.
#'   Must have one element per free cost parameter — two under
#'   `cost = "geometric"`, `m` under `cost = "free"` — and is checked.
#' @param subset logical or integer vector selecting rows of `data`.
#' @param maxit integer, optimiser iteration limit.
#' @param restarts integer, number of times Nelder-Mead is re-run from its
#'   own solution. This matters: a single pass routinely stops short.
#' @param quiet logical; suppress the `beta = 1` warning. The
#'   variance-amplification warning under `cost = "free"` is not suppressed
#'   by it — that one says the specification is meaningless, not that a
#'   default was taken.
#'
#' @return An object of class `recmfit`, a list with `theta`, `par`, `vcov`,
#'   `alpha`, `a`, `k`, `scalars`, `a_f`, `delta`, `residuals`, `fitted`,
#'   `dep`, `index`, `n`, `npar`, `var`, `extras`, `converged`, the column
#'   names `y`, `ystar`, `vars_names` and `growth_name`, plus two closures:
#'
#'   \describe{
#'     \item{`alpha_fn(theta)`}{maps any `theta` to `alpha`; powers
#'       delta-method standard errors.}
#'     \item{`objective(theta)`}{the criterion at any `theta`; powers
#'       [recm_profile()].}
#'   }
#'
#'   **The closures capture `data` by reference**, so a `recmfit` is not
#'   portable across sessions without its data. That is a deliberate trade.
#'
#' @details
#' The fitted decision rule is reported in two equivalent formats.
#'
#' *Explicit* -- the coefficient on the forward sum is 1 by construction:
#'
#' \deqn{\Delta y_t = a_0 (y - y^*)_{t-1}
#'   + \sum_{i=1}^{m-1} a_i \Delta y_{t-i}
#'   + \sum_{i \geq 0} d_i E_{t-1}[\Delta y^*_{t+i}] + \delta'W_t + e_t}
#'
#' *Compressed* (LENS / FRB-US presentation), with
#' \eqn{f_i = d_i / \sum_j d_j} summing to 1 and \eqn{a_f = \sum_j d_j}
#' under the Euler restriction:
#'
#' \deqn{\Delta y_t = a_0 (y - y^*)_{t-1}
#'   + \sum_{i=1}^{m-1} a_i \Delta y_{t-i}
#'   + a_f \sum_{i \geq 0} f_i E_{t-1}[\Delta y^*_{t+i}] + \delta'W_t + e_t}
#'
#' @section Identification:
#' This is inverse optimal control, and weak identification is the normal
#' case rather than the exception. At `m = 2` a fourfold change in `k_1/k_0`
#' cost 2.2% of SSR. **Point estimates alone are not a defensible output** --
#' run [recm_profile()] and [recm_boot()].
#'
#' High-order differences also explode: \eqn{Var((1-L)^j y)} grows like
#' \eqn{\binom{2(j-1)}{j-1}}, about 3.5e10 at `j = 20`, so cost parameters
#' would have to span ten orders of magnitude to contribute comparably.
#' `cost = "geometric"` exists for exactly this reason.
#'
#' @section Growth neutrality:
#' On a balanced growth path with \eqn{\Delta y^* = g}, the exact Euler
#' solution satisfies
#' \eqn{y - y^* = g(1 - \sum_{i \geq 1} a_i - \sum_i d_i)/a_0 < 0}. The
#' optimum runs *below* a growing target: the agent optimally lags, because
#' moving is costly. This is a property of the optimum, not a bug. Setting
#' `growth` adds a correction that removes it, matching FRB/US practice --
#' a deliberate departure from the exact optimum, and it should not be
#' presented as "the solution".
#'
#' @section Instrument dating:
#' `iv_lag` trades bias against variance, and the default deliberately buys
#' bias protection at a real cost in precision.
#'
#' Monte Carlo on the reference DGP (`m = 2`, `var_lags = 2`, so the default
#' `iv_lag` is 4; RMSE for \eqn{a_0}, 60 replications at `T = 178` and 12 at
#' `T = 3000`):
#'
#' \tabular{lrrr}{
#'   \strong{measurement error in ystar} \tab \strong{NLS} \tab
#'     \strong{`iv_lag = 1`} \tab \strong{`iv_lag = 4`} \cr
#'   none, `T = 178`      \tab 0.0204 \tab 0.0212 \tab 0.0541 \cr
#'   none, `T = 3000`     \tab 0.0037 \tab 0.0035 \tab 0.0273 \cr
#'   severe, `T = 178`    \tab 0.0889 \tab 0.0880 \tab 0.0549 \cr
#'   severe, `T = 3000`   \tab 0.0871 \tab 0.0876 \tab 0.0132
#' }
#'
#' The bottom row is the one that matters. With `sd(nu)` at 0.02, NLS and
#' `iv_lag = 1` both sit at a bias of about +0.087 that **does not shrink
#' with `T`** -- they are inconsistent, and the contaminated instruments buy
#' nothing over NLS. At the default dating the bias is -0.009 and falling.
#' The two failure modes are not symmetric: dating too early gives a wrong
#' answer that more data never fixes, dating correctly gives a noisy answer
#' that more data does fix. That is why the rule is the default even though
#' it costs a factor of 2.6 in RMSE when there is no measurement error to
#' correct.
#'
#' Validity was checked directly: at `iv_lag = var_lags + 2` the residual is
#' orthogonal to every instrument column (`p > 0.05` on all of them at
#' `T = 3000`).
#'
#' The exemption matters as much as the dating. Lagging the *whole* state
#' vector, exogenous expectations variables included, destroys
#' identification -- the criterion goes numerically flat and \eqn{a_0}
#' stalls near -0.11 whatever `T` is. Only the `y`/`ystar`-derived columns
#' carry the rule; the rest stay at `t-1`, where they are strong.
#'
#' Set `iv_lag = 1` if and only if `ystar` is measured exactly.
#'
#' @section Standard errors:
#' `Z_t` is a generated regressor: it depends on the estimated auxiliary VAR
#' *and* on `ystar`, which usually comes from a first-stage cointegrating
#' regression. In Monte Carlo the conditional standard errors came in at
#' roughly half the true sampling standard deviation for the PAC parameters,
#' while the standard error on an ordinary exogenous regressor was correct.
#' Use [recm_boot()].
#'
#' @seealso [recm_profile()], [recm_boot()], [equation()], [lead_weights()];
#'   `docs/02-math-spec.md` for the underlying mathematics.
#'
#' @examples
#' set.seed(1)
#' n <- 160
#' dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5,
#'                                 method = "recursive"))
#' dat <- data.frame(ystar = cumsum(dys))
#' dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))
#'
#' fit <- recm_estimate("y", "ystar", beta = 0.995, data = dat, m = 2,
#'                      var_lags = 2, restarts = 1, quiet = TRUE)
#' print(fit)
#'
#' @export
recm_estimate <- function(y, ystar, vars = NULL, beta = 1, data,
                          m = 3, cost = c("geometric", "free"),
                          expectations = NULL, var_lags = 4,
                          diff_exp = TRUE,
                          method = c("nls", "gmm"), instruments = NULL,
                          iv_lag = NULL,
                          free_forward = FALSE, growth = NULL,
                          hac_lags = NULL, start = NULL, subset = NULL,
                          maxit = 4000, restarts = 3, quiet = FALSE) {

  cl <- match.call()
  cost <- match.arg(cost)
  method <- match.arg(method)
  stopifnot(is.character(y), is.character(ystar), m >= 1,
            beta > 0, beta <= 1)
  if (beta == 1 && !quiet) {
    warning("beta = 1: no discounting. The forward sum converges only ",
            "through the adjustment cost roots. Consider beta = 1/(1+r) ",
            "with a real quarterly rate.", call. = FALSE)
  }
  # Deliberately NOT gated on `quiet`, unlike the beta warning above: this
  # one says the specification is meaningless, not that a default was taken.
  cspec <- .cost_spec(cost)
  cwarn <- cspec$warn(m)
  if (!is.null(cwarn)) warning(cwarn, call. = FALSE)
  if (!is.null(subset)) data <- data[subset, , drop = FALSE]

  yv <- data[[y]]
  if (is.null(yv)) stop("column '", y, "' not found in `data`")
  ysv <- data[[ystar]]
  if (is.null(ysv)) stop("column '", ystar, "' not found in `data`")
  Tn <- length(yv)

  W  <- .getcols(vars, data)
  nw <- if (is.null(W)) 0L else ncol(W)
  gr <- if (is.null(growth)) rep(0, Tn) else data[[growth]]
  do_growth <- !is.null(growth)

  dy  <- c(NA, diff(yv))
  dys <- c(NA, diff(ysv))
  ecm <- .lagv(yv - ysv, 1)
  DYL <- if (m > 1) {
    sapply(seq_len(m - 1), function(j) .lagv(dy, j))
  } else {
    NULL
  }
  if (!is.null(DYL)) {
    DYL <- matrix(DYL, nrow = Tn)
    colnames(DYL) <- paste0("dy.l", seq_len(m - 1))
  }

  # ---- auxiliary VAR, estimated once and held fixed ----
  Xe <- cbind(d.ystar = dys)
  Ex <- .getcols(expectations, data)
  if (!is.null(Ex)) {
    if (diff_exp) Ex <- apply(Ex, 2, function(v) c(NA, diff(v)))
    Xe <- cbind(Xe, Ex[, setdiff(colnames(Ex), colnames(Xe)), drop = FALSE])
  }
  vc <- .var_companion(Xe, var_lags)

  # ---- expectations mechanism, built once ----
  # Built here, outside the optimiser, for the same reason the VAR itself
  # is: `zmech$support` fixes the estimation sample, and a sample that moved
  # with theta would make the criterion incomparable across trial values.
  # See the mechanism contract in expectations.R.
  zmech <- .zmech_var(vc, beta)

  # ---- GMM instruments ----
  # Only the columns built from `y` or `ystar` carry the dating rule: the
  # composite error carries measurement error in `ystar` at t-1 ..
  # t-(var_lags+1), so ecm, lagged dy and the d(ystar) states are dated
  # t - iv_lag. Everything else in the state vector is exogenous with
  # respect to that error and stays at t-1, where it is far stronger.
  #
  # Built HERE, before `ok`, so the longer lags shorten the estimation
  # sample. Built after `ok` they would carry leading NAs into the sample
  # and be dropped wholesale by the finiteness filter below -- a silent
  # loss of every automatic instrument.
  Ziv <- NULL
  if (method %in% .METHODS_IV) {
    if (is.null(iv_lag)) iv_lag <- var_lags + 2L
    iv_lag <- as.integer(iv_lag)
    if (is.na(iv_lag) || iv_lag < 1L) {
      stop("`iv_lag` must be a positive integer")
    }
    if (iv_lag < var_lags + 2L && !quiet) {
      warning("iv_lag = ", iv_lag, " < var_lags + 2 = ", var_lags + 2L,
              ": the automatic instruments are dated inside the window the ",
              "composite error contaminates when `ystar` is measured with ",
              "error, and `ystar` usually comes from a first-stage ",
              "regression. Valid only if `ystar` is measured exactly.",
              call. = FALSE)
    }
    lagmat <- function(j) {
      rbind(matrix(NA_real_, j, vc$n),
            vc$states[seq_len(Tn - j), , drop = FALSE])
    }
    S1 <- lagmat(1L)
    Sk <- if (iv_lag == 1L) S1 else lagmat(iv_lag)
    # d(ystar) is always column 1 of X, so its state columns are the ones
    # named d.ystar.l* -- see the layout note in expectations.R.
    is_ys <- grepl("^d\\.ystar\\.l", colnames(vc$states))
    Sm <- S1
    Sm[, is_ys] <- Sk[, is_ys]
    # Name by the TOTAL lag. A state column "V.lj" held at instrument lag L
    # is V at t-(j+L); calling it "V.lj.lL" invites off-by-one reading.
    vlag <- as.integer(sub("^.*\\.l", "", colnames(vc$states)[-1]))
    Sm <- Sm[, -1, drop = FALSE]
    colnames(Sm) <- paste0(sub("\\.l\\d+$", "", colnames(vc$states)[-1]),
                           ".l", vlag + ifelse(is_ys[-1], iv_lag, 1L))

    DY_iv <- if (m > 1) {
      z <- matrix(sapply(seq_len(m - 1),
                         function(j) .lagv(dy, iv_lag + j - 1L)), nrow = Tn)
      colnames(z) <- paste0("dy.l", iv_lag + seq_len(m - 1) - 1L)
      z
    }
    Ziv <- cbind(const = 1, ecm.l = .lagv(yv - ysv, iv_lag), DY_iv, Sm)
    colnames(Ziv)[2] <- paste0("ecm.l", iv_lag)
    IVx <- .getcols(instruments, data)
    if (!is.null(IVx)) Ziv <- cbind(Ziv, IVx)
  }

  # ---- sample ----
  # The mechanism decides which rows it can supply Z for; as an NA-coded
  # column that lets complete.cases() do the rest. This used to be the whole
  # `Slag` matrix, which was the VAR mechanism's answer to the same question
  # spelled out at the call site.
  zsup <- ifelse(zmech$support, 0, NA_real_)
  parts <- cbind(dy, ecm, DYL, zsup, gr, Ziv)
  if (nw) parts <- cbind(parts, W)
  ok <- which(stats::complete.cases(parts))
  if (length(ok) < 5 * (m + nw)) {
    stop("too few usable observations (", length(ok), ") for m = ", m,
         " and ", nw, " extra regressors")
  }

  # ---- residual as a function of the nonlinear parameters ----
  build <- function(th) {
    kk <- .k_from_theta(th, m, cost)
    if (any(!is.finite(kk)) || any(kk <= 0)) return(NULL)
    # A non-convergent Riccati solve is rejected, not approximated: the
    # per-step gain change understates the true error badly for high m with
    # slowly decaying costs. See the note on .TOL_RICCATI.
    lq <- .lq_alpha(kk, beta)
    if (is.null(lq) || !lq$converged || lq$maxeig >= 1) return(NULL)
    al <- lq$alpha
    s <- .scalars(al, beta)
    if (is.null(s)) return(NULL)
    # Z comes IN from the mechanism; this function no longer knows how it
    # was made. The admissibility test that used to live here --
    # rho(G)rho(H) >= 1 -- went with it, because it is the VAR mechanism's
    # condition and not every mechanism's. See expectations.R.
    Zm <- zmech$z(al, s)
    if (is.null(Zm)) return(NULL)
    aa <- .alpha_to_a(al)
    corr <- if (do_growth) (1 - sum(aa[-1]) - s$sum_d) * gr else 0
    list(alpha = al, a = aa, s = s, k = kk,
         Zraw = Zm + corr, Znorm = (Zm + corr) / s$sum_d)
  }

  # Linear block, concentrated out: [a_f on Znorm if free] + [delta on W].
  # This keeps Nelder-Mead in 2 dimensions instead of 2 + p, which matters
  # because it degrades badly above roughly 10 parameters. Concentration is
  # an optimisation device, not an inference shortcut -- the full parameter
  # vector is reassembled afterwards for the covariance.
  linX <- function(b) {
    M <- NULL
    if (free_forward) M <- cbind(M, a_f = b$Znorm[ok])
    if (nw) M <- cbind(M, W[ok, , drop = FALSE])
    M
  }
  core <- function(b) {
    v <- b$a[1] * ecm[ok]
    if (m > 1) v <- v + drop(DYL[ok, , drop = FALSE] %*% b$a[-1])
    if (!free_forward) v <- v + b$Zraw[ok]
    v
  }

  resid_full <- function(th) {
    b <- build(th)
    if (is.null(b)) return(NULL)
    u <- dy[ok] - core(b)
    Xl <- linX(b)
    if (is.null(Xl)) {
      return(list(e = u, lin = numeric(0), b = b, Xl = NULL))
    }
    f <- stats::lm.fit(Xl, u)
    list(e = f$residuals, lin = f$coefficients, b = b, Xl = Xl)
  }
  rfun <- function(th) {
    r <- resid_full(th)
    if (is.null(r)) rep(NA_real_, length(ok)) else r$e
  }

  # ---- starting values ----
  # Both the default and the fallback grid come from the parameterisation,
  # not from a branch here. They used to be written as `if geometric ... else
  # <the "free" shape>`, which handed any other parameterisation a theta of
  # length m; .k_from_theta() then read it positionally and returned numbers.
  npar_th <- cspec$npar(m)
  if (is.null(start)) {
    th0 <- cspec$start(m)
  } else {
    if (length(start) != npar_th) {
      stop("`start` has length ", length(start), " but cost = \"", cost,
           "\" at m = ", m, " takes ", npar_th, " parameter",
           if (npar_th == 1) "" else "s", ": ",
           paste(cspec$names(m), collapse = ", "))
    }
    th0 <- start
  }
  if (is.null(build(th0))) {
    grid <- cspec$grid(m)
    for (i in seq_len(nrow(grid))) {
      if (!is.null(build(as.numeric(grid[i, ])))) {
        th0 <- as.numeric(grid[i, ])
        break
      }
    }
  }
  if (is.null(build(th0))) {
    stop("no admissible starting value found: every trial theta gave ",
         "rho(G)*rho(H) >= 1 or a non-convergent Riccati solve. ",
         "Reduce m or var_lags.")
  }

  # ---- objective ----
  ssr <- function(th) {
    e <- rfun(th)
    if (any(!is.finite(e))) .PENALTY_INADMISSIBLE else sum(e^2)
  }

  if (method == "nls") {
    fit <- stats::optim(th0, ssr, method = "Nelder-Mead",
                        control = list(maxit = maxit, reltol = 1e-12))
    for (i in seq_len(restarts)) {
      fit <- stats::optim(fit$par, ssr, method = "Nelder-Mead",
                          control = list(maxit = maxit, reltol = 1e-14))
    }
    th <- fit$par
    r <- resid_full(th)
    npar <- length(th) + length(r$lin)
    gfun <- function(p) {
      thx <- p[seq_along(th)]
      rr <- resid_full(thx)
      if (is.null(rr)) return(rep(NA_real_, length(ok)))
      lin <- p[-seq_along(th)]
      rr$e + if (length(lin)) drop(rr$Xl %*% (rr$lin - lin)) else 0
    }
    par_all <- c(th, r$lin)
    J <- .jacnum(gfun, par_all)
    s2 <- sum(r$e^2) / (length(ok) - npar)
    V <- tryCatch(s2 * solve(crossprod(J)),
                  error = function(e) matrix(NA_real_, npar, npar))
    extras <- list()
    objective <- ssr
  } else if (method == "gmm") {
    n_supplied <- ncol(Ziv)
    Zi <- Ziv[ok, , drop = FALSE]
    # Column 1 is the intercept and is exempt from the variance filter --
    # that filter exists to remove degenerate columns, and a zero-variance
    # constant is the one column that is meant to be constant. It used to
    # drop it, silently deleting the E[e_t] = 0 moment.
    usable <- apply(Zi, 2, function(z) {
      all(is.finite(z)) && stats::sd(z) > 1e-12
    })
    usable[1] <- TRUE
    Zi <- Zi[, usable, drop = FALSE]
    qq <- qr(Zi)
    if (qq$rank < ncol(Zi)) {
      Zi <- Zi[, sort(qq$pivot[seq_len(qq$rank)]), drop = FALSE]
    }
    Q <- function(th, Wm) {
      e <- rfun(th)
      if (any(!is.finite(e))) return(.PENALTY_INADMISSIBLE)
      g <- colMeans(Zi * e)
      length(e) * drop(t(g) %*% Wm %*% g)
    }
    W1 <- solve(crossprod(Zi) / nrow(Zi))
    s1 <- stats::optim(th0, Q, Wm = W1, method = "Nelder-Mead",
                       control = list(maxit = maxit, reltol = 1e-12))
    W2 <- tryCatch(solve(.nw(Zi * rfun(s1$par), hac_lags)),
                   error = function(e) W1)
    fit <- s1
    for (i in seq_len(restarts + 1)) {
      fit <- stats::optim(fit$par, Q, Wm = W2, method = "Nelder-Mead",
                          control = list(maxit = maxit, reltol = 1e-14))
    }
    th <- fit$par
    r <- resid_full(th)
    par_all <- c(th, r$lin)
    gmom <- function(p) {
      thx <- p[seq_along(th)]
      rr <- resid_full(thx)
      if (is.null(rr)) return(rep(NA_real_, ncol(Zi)))
      lin <- p[-seq_along(th)]
      e <- rr$e + if (length(lin)) drop(rr$Xl %*% (rr$lin - lin)) else 0
      colMeans(Zi * e)
    }
    Gm <- .jacnum(gmom, par_all)
    V <- tryCatch(solve(t(Gm) %*% W2 %*% Gm) / nrow(Zi),
                  error = function(e) {
                    matrix(NA_real_, length(par_all), length(par_all))
                  })
    gb <- gmom(par_all)
    Jst <- nrow(Zi) * drop(t(gb) %*% W2 %*% gb)
    Jdf <- ncol(Zi) - length(par_all)
    # First-stage F for the error-correction regressor, on the EXCLUDED
    # instruments only. At iv_lag = 1 the ecm instruments itself, and
    # regressing a variable on itself reports F ~ 1e31 -- a perfect fit, not
    # a strong instrument. Drop the intercept and any column collinear with
    # the regressor; with nothing left the diagnostic does not exist.
    endog <- ecm[ok]
    excl <- Zi[, -1, drop = FALSE]
    keep <- apply(excl, 2, function(z) {
      stats::sd(z) > 1e-12 && abs(stats::cor(z, endog)) < 1 - 1e-10
    })
    excl <- excl[, keep, drop = FALSE]
    fsF <- NA_real_
    if (ncol(excl)) {
      fs <- stats::lm(endog ~ excl)
      fs0 <- stats::lm(endog ~ 1)
      fsF <- ((sum(stats::resid(fs0)^2) - sum(stats::resid(fs)^2)) /
                ncol(excl)) / (sum(stats::resid(fs)^2) / fs$df.residual)
    }
    if (is.finite(fsF) && fsF < 10 && !quiet) {
      warning("first-stage F = ", signif(fsF, 4), " < 10: weak ",
              "instruments. A shorter expectations VAR buys instrument ",
              "strength at the cost of forecast accuracy.", call. = FALSE)
    }
    extras <- list(J = Jst, Jdf = Jdf,
                   hansen_p = stats::pchisq(Jst, max(Jdf, 1),
                                            lower.tail = FALSE),
                   first_stage_F = fsF,
                   n_instruments = ncol(Zi),
                   n_instruments_supplied = n_supplied,
                   iv_lag = iv_lag)
    objective <- function(t2) Q(t2, W2)
  } else {
    # Not reachable from user code -- match.arg() above rejects an unknown
    # name first. This catches the developer who added an estimator to the
    # `method` argument and not here, which previously fell through to the
    # GMM branch. Do not delete it as dead code.
    stop("estimator '", method, "' is accepted by the `method` argument but ",
         "has no branch in recm_estimate(). It must set th, r, par_all, V, ",
         "fit, extras and objective -- see the estimator extension point in ",
         "docs/01-architecture.md.")
  }

  names(par_all) <- c(.theta_names(m, cost),
                      if (free_forward) "a_f",
                      if (nw) colnames(W))
  dimnames(V) <- list(names(par_all), names(par_all))

  b <- r$b
  a_f <- if (free_forward) unname(r$lin[1]) else b$s$sum_d

  out <- list(call = cl, method = method, cost = cost, m = m, beta = beta,
              free_forward = free_forward, growth = do_growth,
              theta = th, par = par_all, vcov = V,
              alpha = b$alpha, a = b$a, k = b$k, scalars = b$s,
              a_f = a_f,
              delta = if (nw) {
                r$lin[length(r$lin) - nw + seq_len(nw)]
              } else {
                NULL
              },
              residuals = r$e, fitted = dy[ok] - r$e, dep = dy[ok],
              index = ok, n = length(ok), npar = length(par_all),
              var = vc, extras = extras, converged = fit$convergence == 0,
              vars_names = if (nw) colnames(W) else character(0),
              objective = objective,
              y = y, ystar = ystar, growth_name = growth,
              alpha_fn = function(t2) {
                bb <- build(t2)
                if (is.null(bb)) rep(NA_real_, m) else bb$alpha
              })
  class(out) <- "recmfit"
  out
}


# ---- the estimation frame ----
#
# `objective` is a closure over recm_estimate()'s evaluation frame, so that
# frame stays alive for as long as the fit does, and it is where `data` and
# the resolved specification live. recm_boot() and predict.recmfit() both
# have to rebuild the design on resampled or new observations, and both read
# it from here rather than from a copy stored on the object: a copy could go
# stale against the closures, which are what recm_profile() and the
# delta-method standard errors actually evaluate.
#
# Reading it post-`subset` is the point, not an accident -- the design is
# rebuilt on the rows that were estimated on.
#
# The trade is the documented one (docs/01, "Object model"): a `recmfit` is
# not portable across sessions without its data. Both callers previously
# reached in by name with no check beyond `data`; a frame missing anything
# else they read would have produced predictions from a half-built design
# rather than an error.
.fit_env <- function(object, need = character(0)) {
  env <- environment(object$objective)
  if (!is.environment(env) || is.null(env$data)) {
    stop("the fit's closures no longer carry their data; a `recmfit` is ",
         "not portable across sessions (see docs/01-architecture.md)")
  }
  miss <- need[!vapply(need, exists, logical(1), envir = env,
                       inherits = FALSE)]
  if (length(miss)) {
    stop("the fit's estimation frame is missing ",
         paste(miss, collapse = ", "),
         "; it did not come from this version of recm_estimate()")
  }
  env
}
