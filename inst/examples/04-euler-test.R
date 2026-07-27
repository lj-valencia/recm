# -----------------------------------------------------------------------------
# 04-euler-test.R
# Testing the Euler equation restriction and applying growth neutrality.
#
# Shows:
#   1. Estimating with free_forward = TRUE to un-restrict the forward loading a_f
#   2. Comparing freely estimated a_f against the theoretical restriction sum(d_i)
#   3. Using the growth argument to add growth neutrality when the target has trend growth
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

# 1. Generate data with a trend growth rate g = 0.01 (1% per period)
set.seed(20260727)
n <- 160
g_trend <- 0.01

d_ystar <- rnorm(n, mean = g_trend, sd = 0.004)
ystar   <- cumsum(d_ystar)
y       <- ystar - 0.04 + cumsum(rnorm(n, mean = 0, sd = 0.002))

data_euler <- data.frame(
  y = y,
  ystar = ystar,
  g = rep(g_trend, n)
)

cat("=== 1. Estimating Under Structural Euler Restriction (free_forward = FALSE) ===\n")
fit_restricted <- recm_estimate(
  y = "y",
  ystar = "ystar",
  beta = 0.995,
  data = data_euler,
  m = 2,
  cost = "geometric",
  var_lags = 2,
  free_forward = FALSE,
  quiet = TRUE
)

cat("Restricted a_f (= sum(d_i)):", round(fit_restricted$a_f, 5), "\n")
print(equation(fit_restricted, format = "compressed"))

cat("\n=== 2. Estimating with Un-restricted Forward Loading (free_forward = TRUE) ===\n")
# Estimating a_f freely is the primary empirical test of the PAC specification.
# If the structural model holds, the freely estimated a_f will be close to sum(d_i).
fit_free <- recm_estimate(
  y = "y",
  ystar = "ystar",
  beta = 0.995,
  data = data_euler,
  m = 2,
  cost = "geometric",
  var_lags = 2,
  free_forward = TRUE,
  quiet = TRUE
)

cat("Theoretical sum(d_i) under implied alpha:", round(fit_free$scalars$sum_d, 5), "\n")
cat("Freely estimated a_f:                   ", round(fit_free$a_f, 5), "\n")
print(equation(fit_free, format = "compressed"))

cat("\n=== 3. Growth-Neutral Specification (growth = 'g') ===\n")
# On a balanced growth path (d.ystar = g > 0), the exact Euler optimum runs below target:
# (y - ystar) = g * (1 - sum a_i - sum d_i) / a0 < 0.
# Specifying growth = "g" adds the FRB/US growth correction term to restore steady-state neutrality.
fit_growth <- recm_estimate(
  y = "y",
  ystar = "ystar",
  growth = "g",
  beta = 0.995,
  data = data_euler,
  m = 2,
  cost = "geometric",
  var_lags = 2,
  quiet = TRUE
)

print(summary(fit_growth))

cat("\n[04-euler-test.R completed successfully]\n")
