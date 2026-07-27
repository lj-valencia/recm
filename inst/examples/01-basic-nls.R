# -----------------------------------------------------------------------------
# 01-basic-nls.R
# Minimal fit of a Rational Error Correction Model (RECM) using NLS.
#
# Shows:
#   1. Simulating or loading data for y (decision variable) and ystar (target)
#   2. Model estimation via recm_estimate() with method = "nls"
#   3. Using standard S3 methods (print, summary, coef, vcov, residuals, fitted, confint, logLik)
#   4. Reporting fitted decision rules in dual format (explicit and compressed) via equation()
#   5. Extracting lead weights and plotting the model via lead_weights() and plot()
# -----------------------------------------------------------------------------

# Load recm package (or devtools::load_all() during package development)
if (!requireNamespace("recm", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("Please install the 'recm' package to run this example.")
  }
} else {
  library(recm)
}

# 1. Simulate data driven by a PAC process
set.seed(2026)
n <- 160 # 40 years of quarterly data

# Frictionless target ystar: AR(1) in growth rates around mean quarterly growth of ~1%
d_ystar <- as.numeric(stats::filter(rnorm(n, mean = 0.0104, sd = 0.004),
                                    filter = 0.5, method = "recursive"))
ystar <- cumsum(d_ystar)

# Decision variable y: optimal adjustment towards target with noise
y <- ystar - 0.05 + cumsum(rnorm(n, mean = 0, sd = 0.002))
data_example <- data.frame(y = y, ystar = ystar)

# 2. Fit RECM via Nonlinear Least Squares (NLS)
#    - y: decision variable (in levels)
#    - ystar: frictionless target (in levels)
#    - beta: calibrated discount factor (e.g., 0.995 for 2% annual real discount rate)
#    - m: order of adjustment cost polynomial (m = 2 means 1 lag of dy)
#    - cost: "geometric" (default) or "free"
#    - var_lags: lag order of auxiliary VAR for expectations
fit_nls <- recm_estimate(
  y = "y",
  ystar = "ystar",
  beta = 0.995,
  data = data_example,
  m = 12,
  cost = "geometric",
  var_lags = 2,
  method = "nls",
  quiet = TRUE
)

# 3. Standard S3 accessors
cat("\n=== 1. Basic Print Output ===\n")
print(fit_nls)

cat("\n=== 2. Summary Output ===\n")
print(summary(fit_nls))

cat("\n=== 3. Extracted Coefficients & Covariance ===\n")
cat("Coefficients:\n")
print(coef(fit_nls))

cat("\nCovariance Matrix (vcov):\n")
print(vcov(fit_nls))

cat("\nConfidence Intervals (95%):\n")
print(confint(fit_nls))

cat("\nModel Log-Likelihood:\n")
print(logLik(fit_nls))

cat("\nFirst 6 Residuals:\n")
print(head(residuals(fit_nls)))

cat("\nFirst 6 Fitted Values (dy_t):\n")
print(head(fitted(fit_nls)))

# 4. Dual-format reporting via equation()
#    - Explicit format: d.y_t = a0*(y - ystar)_{t-1} + a1*d.y_{t-1} + Z_t + e_t (Z_t coefficient = 1)
#    - Compressed format: d.y_t = a0*(y - ystar)_{t-1} + a1*d.y_{t-1} + a_f * sum(f_i * E[d.ystar_{t+i}]) + e_t (sum f_i = 1)
cat("\n=== 4. Dual-Format Decision Rule Reporting ===\n")
eq_res <- equation(fit_nls, format = "both")

# 5. Lead weights & Plots
cat("\n=== 5. Lead Weights (d_i and normalised f_i) ===\n")
lw_unnorm <- lead_weights(fit_nls, horizon = 8, normalised = FALSE)
cat("Unnormalised lead weights (d_i):\n")
print(lw_unnorm)

lw_norm <- lead_weights(fit_nls, horizon = 8, normalised = TRUE)
cat("\nNormalised lead weights (f_i, sum to 1):\n")
print(lw_norm)

# Display 4-panel diagnostic plot if running interactively
if (interactive()) {
  plot(fit_nls, horizon = 24)
}

cat("\n[01-basic-nls.R completed successfully]\n")
