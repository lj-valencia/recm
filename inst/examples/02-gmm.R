# -----------------------------------------------------------------------------
# 02-gmm.R
# Generalized Method of Moments (GMM) estimation of a RECM equation.
#
# Shows:
#   1. Two-step GMM estimation (method = "gmm") with Newey-West weighting matrix
#   2. Instrument dating (iv_lag): why default iv_lag = var_lags + 2 is required
#      when ystar is measured with error (e.g. from cointegrating regression) vs iv_lag = 1
#   3. Adding linear regressors (vars), auxiliary expectations (expectations),
#      and custom excluded instruments (instruments)
#   4. Inspecting GMM diagnostics (Hansen J test, first-stage F statistic, instrument counts)
# -----------------------------------------------------------------------------

# Load recm package
if (!requireNamespace("recm", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("Please install the 'recm' package to run this example.")
  }
} else {
  library(recm)
}

# 1. Generate synthetic data with measurement error in ystar
set.seed(20260727)
n <- 180

# True frictionless target & decision variable
d_ystar_true <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
ystar_true   <- cumsum(d_ystar_true)
y            <- ystar_true - 0.05 + cumsum(rnorm(n, 0, 0.003))

# Observed ystar has measurement error nu_t (e.g. estimated cointegrating vector)
nu           <- rnorm(n, 0, 0.004)
ystar_obs    <- ystar_true + nu

# Exogenous regressor (e.g., world demand or financial condition index)
x_exog       <- rnorm(n, 0, 0.01)

# Additional external instrument (e.g., interest rate spread)
z_ext        <- rnorm(n, 0, 0.02)

data_gmm <- data.frame(
  y = y,
  ystar = ystar_obs,
  x_exog = x_exog,
  z_ext = z_ext
)

cat("=== 1. Standard GMM Estimation (Default Instrument Dating) ===\n")
# In default GMM, iv_lag = var_lags + 2 (here 2 + 2 = 4).
# This dates automatic y/ystar-derived instruments at t-4 to avoid contamination
# from measurement error in ystar at t-1..t-3.
fit_gmm_default <- recm_estimate(
  y = "y",
  ystar = "ystar",
  beta = 0.995,
  data = data_gmm,
  m = 12,
  cost = "geometric",
  var_lags = 2,
  method = "gmm",
  quiet = TRUE
)

print(summary(fit_gmm_default))

cat("\n=== 2. GMM Diagnostics (Hansen J & First-Stage F) ===\n")
ex_def <- fit_gmm_default$extras
cat("Instruments used:", ex_def$n_instruments, "\n")
cat("Hansen J stat:", round(ex_def$J, 4), "(df =", ex_def$Jdf, ", p-value =", round(ex_def$hansen_p, 4), ")\n")
cat("First-stage F stat:", round(ex_def$first_stage_F, 4), "\n")

cat("\n=== 3. Comparing Instrument Dating Rules ===\n")
# If ystar were measured WITHOUT error, iv_lag = 1 could be used for extra efficiency.
# But under measurement error, iv_lag = 1 leads to inconsistency (biased estimates).
fit_gmm_lag1 <- recm_estimate(
  y = "y",
  ystar = "ystar",
  beta = 0.995,
  data = data_gmm,
  m = 2,
  cost = "geometric",
  var_lags = 2,
  method = "gmm",
  iv_lag = 1,
  quiet = TRUE
)

cat("Default dating (iv_lag = 4)   a0 =", round(fit_gmm_default$a[1], 4), "  a1 =", round(fit_gmm_default$a[2], 4), "\n")
cat("Contaminated (iv_lag = 1)     a0 =", round(fit_gmm_lag1$a[1], 4), "  a1 =", round(fit_gmm_lag1$a[2], 4), "\n")

cat("\n=== 4. GMM with Exogenous Regressors & Custom Instruments ===\n")
# - vars: linear regressor entering directly (e.g. x_exog)
# - instruments: user-supplied excluded instrument (e.g. z_ext)
fit_gmm_full <- recm_estimate(
  y = "y",
  ystar = "ystar",
  vars = "x_exog",
  instruments = "z_ext",
  beta = 0.995,
  data = data_gmm,
  m = 2,
  cost = "geometric",
  var_lags = 2,
  method = "gmm",
  quiet = TRUE
)

print(equation(fit_gmm_full, format = "compressed"))

cat("Linear coefficient (delta on x_exog):\n")
print(fit_gmm_full$delta)

cat("\n[02-gmm.R completed successfully]\n")
