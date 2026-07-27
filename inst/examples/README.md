# Runnable end-to-end examples

Scripts here are installed with the package and may use suggested packages (`testthat`, `devtools`), unlike anything under `R/`.

| Script | Shows | Key Functions |
|---|---|---|
| [`01-basic-nls.R`](file:///c:/Users/ABC/Documents/recm/inst/examples/01-basic-nls.R) | Minimal fit using NLS, S3 accessors, dual reporting formats (explicit & compressed), lead weights, plotting | `recm_estimate()`, `equation()`, `lead_weights()`, `plot()` |
| [`02-gmm.R`](file:///c:/Users/ABC/Documents/recm/inst/examples/02-gmm.R) | Two-step GMM estimation, instrument dating rules (`iv_lag`), exogenous regressors, custom instruments, GMM diagnostics | `recm_estimate(method = "gmm")`, Hansen J test, First-stage F |
| [`03-identification.R`](file:///c:/Users/ABC/Documents/recm/inst/examples/03-identification.R) | Objective criterion profiling (`recm_profile()`), near-flat identification diagnostics, dynamic bootstrap standard errors (`recm_boot()`) | `recm_profile()`, `recm_boot()` |
| [`04-euler-test.R`](file:///c:/Users/ABC/Documents/recm/inst/examples/04-euler-test.R) | Testing the PAC Euler restriction (`free_forward = TRUE`), growth-neutrality adjustment (`growth`) | `recm_estimate(free_forward = TRUE)`, `recm_estimate(growth = ...)` |

Running any script from the command line:

```bash
Rscript inst/examples/01-basic-nls.R
Rscript inst/examples/02-gmm.R
Rscript inst/examples/03-identification.R
Rscript inst/examples/04-euler-test.R
```

*Note: The worked "cointegrating regression to PAC equation" walkthrough is roadmap item R-9 and belongs in a vignette.*
