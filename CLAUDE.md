# CLAUDE.md - Project Memory & Directives

## Project Overview
This package implements rational error correction (REC) models for time series estimation, as outlined in Tinsley (2002). The entire codebase is written in R.

## Build and Test Commands
- **Install dependencies:** `Rscript -e "devtools::install_deps(dependencies = TRUE)"`
- **Build package:** `R CMD build .`
- **Check package (CRAN-style):** `R CMD check --as-cran .`
- **Run all tests:** `Rscript -e "devtools::test()"`
- **Run a single test file:** `Rscript -e "devtools::test_active_file('tests/testthat/test-<name>.R')"`
- **Lint:** `Rscript -e "lintr::lint_package()"`
- **Style/autofix:** `Rscript -e "styler::style_pkg()"`
- **Rebuild documentation:** `Rscript -e "devtools::document()"`

## Architecture & Code Style
- **Stack:** Base R for all estimation and modeling logic. Use only `stats`, `utils`, and other base/recommended packages for core functionality — avoid tidyverse or other non-base dependencies in package code.
- **Dependencies:** External packages (e.g., `testthat`) are permitted only in `Suggests`, scoped strictly to testing. Never add a package to `Imports` unless base R cannot reasonably accomplish the task.
- **Package Structure:** Follow standard R package layout — exported functions in `R/`, documentation via `roxygen2` (`#'` tags above each function), tests in `tests/testthat/`, `DESCRIPTION` and `NAMESPACE` kept in sync via `devtools::document()`.
- **Estimation Code:** Keep model-fitting, data-transformation, and inference logic in separate files (e.g., `R/estimate.R`, `R/transform.R`, `R/inference.R`) rather than one monolithic script.
- **Error Handling:** Validate inputs at the top of exported functions and fail with informative `stop()` messages. Use `warning()` for recoverable numerical issues (e.g., convergence warnings).
- **Type Safety:** Explicitly check argument classes/dimensions (e.g., `is.numeric`, `is.matrix`) at function entry; do not rely on implicit coercion for core estimation routines.
- **Style:** Follow the tidyverse style guide for naming and formatting (enforced via `styler`/`lintr`), even though tidyverse packages are not used as dependencies.

## Workflow Rules & Guardrails
- **Pre-execution Step:** You must formulate and explain your development plan to the user before making any code modifications.
- **Post-execution Step:** Always run `lintr::lint_package()` and `devtools::check()` following file changes to guarantee clean compilation and CRAN-style compliance.
- **Placeholders:** Never insert `# TODO:` or stub out incomplete logic when fixing issues. Implement the functional code fully.
- **Testing Guardrail:** Do not mark a debugging or feature task as complete until you have successfully executed the corresponding `testthat` suite.
- **Dependency Guardrail:** Before adding any package to `DESCRIPTION`, confirm it cannot be reasonably replaced with base R; testing-only packages belong in `Suggests`, never `Imports`.
- **Secrets Protection:** Never store, embed, or expose API credentials, tokens, or encryption keys in the codebase.
