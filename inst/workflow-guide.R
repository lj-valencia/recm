# -----------------------------------------------------------------------------
# workflow-guide.R
# Executable decision guide: "What to use and when" in the recm package.
#
# Usage:
#   Rscript inst/workflow-guide.R
# -----------------------------------------------------------------------------

# Load recm package
if (!requireNamespace("recm", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("Please install the 'recm' package to run this guide.")
  }
} else {
  library(recm)
}

cat("=====================================================================\n")
cat("          RECM Package Decision Guide: What to Use and When          \n")
cat("=====================================================================\n\n")

cat("---------------------------------------------------------------------\n")
cat("1. CALIBRATING THE DISCOUNT FACTOR (beta)\n")
cat("---------------------------------------------------------------------\n")
cat("RULE: Always calibrate beta outside the estimation. Do not estimate it.\n")
cat("WHY:  beta is not identified jointly with adjustment cost parameters.\n")
cat("      It enters only through the composite operator (1 - beta*F)(1 - L).\n")
cat("HOW:  Use beta = 1 / (1 + r) with a REAL quarterly discount rate r.\n")
cat("      For example, a 2% annual real discount rate gives r ~ 0.005/quarter,\n")
cat("      so set beta = 1 / 1.005 = 0.995.\n")
cat("NOTE: Default beta = 1 implies no discounting (indifference to infinite\n")
cat("      future) and triggers a warning unless quiet = TRUE.\n\n")

cat("---------------------------------------------------------------------\n")
cat("2. CHOOSING THE ESTIMATION METHOD (method = 'nls' vs 'gmm')\n")
cat("---------------------------------------------------------------------\n")
cat("USE method = 'nls' (Nonlinear Least Squares) WHEN:\n")
cat("  - Baseline modeling without measurement error in ystar.\n")
cat("  - Rapid model screening and exploratory fits.\n")
cat("USE method = 'gmm' (Two-step GMM with Newey-West matrix) WHEN:\n")
cat("  - Target variable ystar is estimated from a first-stage cointegrating\n")
cat("    regression and carries measurement/estimation error.\n")
cat("  - Overidentifying restrictions or external instruments are present.\n\n")

cat("---------------------------------------------------------------------\n")
cat("3. COST PARAMETERISATION (cost = 'geometric' vs 'free')\n")
cat("---------------------------------------------------------------------\n")
cat("USE cost = 'geometric' (DEFAULT) WHEN:\n")
cat("  - Polynomial order m > 4 (or for parsimonious 2-parameter specification).\n")
cat("  - Variance of high-order differences Var((1-L)^j y) explodes like\n")
cat("    choose(2(j-1), j-1), requiring cost weights to span ten orders\n")
cat("    of magnitude under free estimation.\n")
cat("USE cost = 'free' WHEN:\n")
cat("  - Polynomial order m <= 4 and unconstrained cost weights per lag are desired.\n")
cat("  - Warns automatically if m > 4.\n\n")

cat("---------------------------------------------------------------------\n")
cat("4. GMM INSTRUMENT DATING RULES (iv_lag)\n")
cat("---------------------------------------------------------------------\n")
cat("RULE: Automatic instruments built from y or ystar are dated t - iv_lag.\n")
cat("DEFAULT iv_lag = var_lags + 2:\n")
cat("  - Protects against measurement error in ystar at t-1..t-(var_lags+1).\n")
cat("  - Prevents inconsistent estimates when ystar comes from cointegration.\n")
cat("EXEMPTION:\n")
cat("  - Exogenous variables (W_t) and auxiliary expectation variables stay at t-1.\n")
cat("  - Lagging whole state vector destroys identification; exemption is crucial.\n")
cat("SET iv_lag = 1 ONLY IF:\n")
cat("  - ystar is known to be measured EXACTLY without error.\n\n")

cat("---------------------------------------------------------------------\n")
cat("5. EULER RESTRICTION TESTING (free_forward)\n")
cat("---------------------------------------------------------------------\n")
cat("DEFAULT free_forward = FALSE:\n")
cat("  - Restricts forward loading a_f = sum(d_i) implied by adjustment costs.\n")
cat("SET free_forward = TRUE WHEN:\n")
cat("  - Formally testing the PAC specification hypothesis.\n")
cat("  - Comparing freely estimated a_f to sum(d_i) evaluated at estimated alpha.\n\n")

cat("---------------------------------------------------------------------\n")
cat("6. TREND GROWTH neutralITY (growth)\n")
cat("---------------------------------------------------------------------\n")
cat("DEFAULT growth = NULL:\n")
cat("  - Exact Euler solution optimum runs below a growing target (y - ystar < 0).\n")
cat("SET growth = 'growth_col' WHEN:\n")
cat("  - Target follows a balanced growth path (e.g. d.ystar = g > 0).\n")
cat("  - Adds FRB/US steady-state neutrality correction term to eliminate gap bias.\n\n")

cat("---------------------------------------------------------------------\n")
cat("7. DIAGNOSTICS & INFERENCE PIPELINE\n")
cat("---------------------------------------------------------------------\n")
cat("CRITICAL INFERENCE WARNING:\n")
cat("  - Conditional standard errors in summary() treat the VAR as known and\n")
cat("    are optimistic (~half true sampling standard deviation).\n")
cat("  - The optimization criterion surface is routinely near-flat.\n")
cat("MANDATORY DIAGNOSTIC STEPS:\n")
cat("  1. Run recm_profile(fit, which = 1) to inspect criterion curvature.\n")
cat("     (rel_range < 0.02 signals weak identification).\n")
cat("  2. Run recm_boot(fit, R = 299) for honest dynamic bootstrap standard errors\n")
cat("     that propagate VAR estimation uncertainty.\n\n")

cat("---------------------------------------------------------------------\n")
cat("SUMMARY WORKFLOW CHEATSHEET\n")
cat("---------------------------------------------------------------------\n")
cat("Step 1: Fit model\n")
cat("  fit <- recm_estimate(y='y', ystar='ystar', beta=0.995, data=df, m=2,\n")
cat("                       var_lags=2, method='nls')\n\n")
cat("Step 2: Inspect decision rules & lead weights\n")
cat("  equation(fit, format='both')\n")
cat("  lead_weights(fit, horizon=12)\n")
cat("  plot(fit)\n\n")
cat("Step 3: Diagnose identification & curvature\n")
cat("  prof <- recm_profile(fit, which=1)\n")
cat("  print(prof)\n\n")
cat("Step 4: Compute bootstrap confidence intervals\n")
cat("  bs <- recm_boot(fit, R=299, seed=2026)\n")
cat("  print(bs$interval)\n")
cat("=====================================================================\n")
