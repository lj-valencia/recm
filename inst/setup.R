# -----------------------------------------------------------------------------
# setup.R
# System setup check, dependency verification, and environment diagnostics for recm.
#
# Usage:
#   Rscript inst/setup.R
# -----------------------------------------------------------------------------

cat("=====================================================================\n")
cat("            RECM Package Setup & Environment Verification            \n")
cat("=====================================================================\n\n")

# 1. R Version Check
r_version <- getRversion()
cat("1. Checking R Version...\n")
cat("   Detected R version:", as.character(r_version), "\n")
if (r_version >= "4.1.0") {
  cat("   [OK] R version meets minimum requirement (>= 4.1.0)\n\n")
} else {
  cat("   [WARNING] R version < 4.1.0. recm requires R >= 4.1.0.\n\n")
}

# 2. Package Dependency Checks
cat("2. Checking Dependencies...\n")
required_pkgs <- c("graphics", "grDevices", "stats")
suggested_pkgs <- c("testthat", "devtools", "lintr")

req_status <- vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
cat("   Required base/stats packages:\n")
for (pkg in names(req_status)) {
  status_str <- if (req_status[[pkg]]) "[OK]" else "[MISSING]"
  cat("     -", sprintf("%-12s", pkg), status_str, "\n")
}

sug_status <- vapply(suggested_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
cat("   Suggested development packages:\n")
for (pkg in names(sug_status)) {
  status_str <- if (sug_status[[pkg]]) "[INSTALLED]" else "[NOT INSTALLED - optional]"
  cat("     -", sprintf("%-12s", pkg), status_str, "\n")
}
cat("\n")

# 3. recm Package Availability
cat("3. Loading recm Package...\n")
loaded_recm <- FALSE
if (requireNamespace("recm", quietly = TRUE)) {
  library(recm)
  cat("   [OK] Loaded installed 'recm' package (v", as.character(utils::packageVersion("recm")), ")\n", sep = "")
  loaded_recm <- TRUE
} else if (requireNamespace("devtools", quietly = TRUE)) {
  cat("   'recm' is not installed in standard library; trying devtools::load_all()...\n")
  tryCatch({
    devtools::load_all(".", quiet = TRUE)
    cat("   [OK] Successfully loaded recm from local source directory.\n")
    loaded_recm <- TRUE
  }, error = function(e) {
    cat("   [ERROR] Failed to load recm source via devtools: ", e$message, "\n")
  })
} else {
  cat("   [ERROR] 'recm' package is not installed and 'devtools' is not available.\n")
}
cat("\n")

# 4. Numerical Engine Sanity Check
if (loaded_recm) {
  cat("4. Running Numerical Engine Sanity Check...\n")
  set.seed(2026)
  n <- 120
  dys <- as.numeric(stats::filter(rnorm(n, 0.0104, 0.004), 0.5, method = "recursive"))
  dat <- data.frame(ystar = cumsum(dys))
  dat$y <- dat$ystar - 0.05 + cumsum(rnorm(n, 0, 0.002))

  fit_test <- tryCatch({
    recm_estimate(
      y = "y", ystar = "ystar", beta = 0.995,
      data = dat, m = 2, var_lags = 2, restarts = 1, quiet = TRUE
    )
  }, error = function(e) NULL)

  if (!is.null(fit_test) && isTRUE(fit_test$converged)) {
    cat("   [OK] Riccati solver & estimator converged successfully.\n")
    cat("   - Implied error correction a0:", round(fit_test$a[1], 4), "\n")
    cat("   - Implied lag coefficient a1: ", round(fit_test$a[2], 4), "\n")
    cat("   - Spectral stability rho(G)*rho(H):", round(fit_test$scalars$rhoG * fit_test$var$rho_dyn, 4), "\n\n")
  } else {
    cat("   [ERROR] Numerical engine sanity test failed or did not converge.\n\n")
  }
}

cat("=====================================================================\n")
cat("Setup verification complete. You are ready to work with 'recm'!\n")
cat("Run 'Rscript inst/workflow-guide.R' for package usage guidelines.\n")
cat("=====================================================================\n")
