# -----------------------------------------------------------------------------
# 03-identification.R
# Diagnostics for weak identification and bootstrap inference.
#
# Shows:
#   1. Why point estimates alone are insufficient (near-flat criterion surface)
#   2. Profiling the criterion via recm_profile() to assess parameter identification
#   3. Bootstrap standard errors via recm_boot() to correct conditional standard
#      errors (which suffer from generated regressor bias in the VAR)
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

# 1. Generate data from a known PAC process
set.seed(2026)
n <- 150

d_ystar <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
ystar   <- cumsum(d_ystar)
y       <- ystar - 0.05 + cumsum(rnorm(n, 0, 0.003))

data_id <- data.frame(y = y, ystar = ystar)

# 2. Fit initial model via NLS
cat("=== 1. Fitting Base Model ===\n")
fit <- recm_estimate(
  y = "y",
  ystar = "ystar",
  beta = 0.995,
  data = data_id,
  m = 2,
  cost = "geometric",
  var_lags = 2,
  restarts = 1,
  quiet = TRUE
)

# 3. Profile the objective criterion
#    recm_profile() re-optimises remaining parameters across a grid for parameter `which`.
#    rel_range measures curvature; values < 0.02 indicate a weakly identified parameter.
cat("\n=== 2. Profiling Criterion Curvature (recm_profile) ===\n")
prof1 <- recm_profile(fit, which = 1, ngrid = 11, span = 2)
print(prof1)

cat("\nProfile Grid & Objective Values:\n")
prof_table <- data.frame(
  theta1_grid = prof1$grid,
  criterion   = prof1$value,
  implied_a0  = prof1$a0
)
print(head(prof_table, 10))

# Plot profile if running interactively
if (interactive()) {
  plot(prof1)
}

# 4. Bootstrap Inference (recm_boot)
#    Conditional SEs from summary() treat the auxiliary VAR as fixed/known and
#    underestimate standard errors by up to ~50%. recm_boot() resamples the VAR
#    and equation residuals recursively to compute honest sampling variability.
cat("\n=== 3. Running Dynamic Bootstrap (recm_boot) ===\n")
cat("Running R = 30 bootstrap replications...\n")
boot_res <- recm_boot(fit, R = 30, var_error = TRUE, eq_error = TRUE, seed = 2026)

cat("\nBootstrap Summary:\n")
cat("Successful replications:", boot_res$R_ok, "/", boot_res$R_ok + boot_res$R_failed, "\n")

cat("\nParameter Standard Errors Comparison:\n")
cond_se <- sqrt(pmax(diag(vcov(fit)), 0))
se_comp <- data.frame(
  Parameter       = names(fit$par),
  Conditional_SE  = cond_se,
  Bootstrap_SE    = boot_res$se
)
print(se_comp)

cat("\n95% Bootstrap Percentile Confidence Intervals:\n")
print(boot_res$interval)

cat("\n[03-identification.R completed successfully]\n")
