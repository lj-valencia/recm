# Runnable end-to-end examples

Scripts here are installed with the package and may use suggested packages
(`testthat`, `devtools`), unlike anything under `R/`.

Planned, blocked on the estimator landing:

| Script | Shows | Depends on |
|---|---|---|
| `01-basic-nls.R` | minimal fit, both reporting formats | `estimate.R`, `report.R` |
| `02-gmm.R` | GMM with explicit instruments, instrument dating | `estimate.R` |
| `03-identification.R` | `recm_profile()`, `recm_boot()`, the flat criterion | `diagnostics.R` |
| `04-euler-test.R` | `free_forward = TRUE` against `sum(d_i)` | `estimate.R` |

The worked "cointegrating regression to PAC equation" walkthrough is roadmap
item R-9 and belongs in a vignette, not here.
