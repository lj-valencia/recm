# 06 — Conventions

## Dependencies

`R/` uses **base R and `stats` only**. No `dplyr`, `data.table`, `ggplot2`,
`Matrix`, or `Rcpp`. `graphics` and `grDevices` for plotting.

Rationale: this package is used inside institutional environments where
installing a dependency tree is a procurement exercise. The numerical work is
small dense linear algebra that base R handles fine.

Suggested packages (`testthat`, `devtools`, `lintr`) are permitted in
`tests/` and `inst/examples/` only.

Never write `library()` inside `R/`. Use `stats::lm.fit()`, not `lm.fit()`.

## Naming

| Kind | Convention | Example |
|---|---|---|
| exported function | `snake_case`, `recm_` prefix if ambiguous | `recm_estimate()`, `lead_weights()` |
| internal helper | leading dot | `.hvec()`, `.alpha_to_a()` |
| S3 method | `generic.class` | `summary.recmfit()` |
| class | lowercase, no separator | `recmfit`, `recm_profile` |
| local variable | short, matching the maths | `alpha`, `beta`, `G`, `H` |

**Match the mathematical notation.** A variable holding `A(L)` coefficients is
`alpha`, not `poly_coefs`. Someone reading the code with docs/02 open should be
able to follow line by line. This overrides any general preference for
descriptive names.

Exception: single-letter names that shadow base R (`c`, `t`, `F`, `T`, `diff`)
are forbidden even when the maths uses them. Use `cc`, `tt`, `Fg`, `dd`.

## Style

- 2-space indent, no tabs
- 80 columns
- `<-` for assignment; `=` only for arguments
- spaces around binary operators, after commas
- `TRUE`/`FALSE` spelled out
- one blank line between functions, two between sections
- section banners as `# ---- name ----` at the top level

## Error handling

Three tiers, and the distinction matters:

1. **Internal builders return `NULL`** on inadmissible parameters. The
   optimiser must be able to reject a trial `theta` cheaply without an
   exception. Never `stop()` inside a function the optimiser calls.

2. **Exported functions `stop()`** with a message naming the violated
   condition and the value that violated it:
   ```r
   stop("rho(G)*rho(H) = ", signif(rg * rh, 4),
        " >= 1; the forward sum does not converge. ",
        "Reduce var_lags or m.")
   ```

3. **`warning()`** for choices that are legal but likely unintended:
   `beta = 1`, `cost = "free"` with `m > 4`, first-stage F below 10,
   spectral factorisation above m = 5.

Never `message()` from `R/`. Never `print()` outside a `print.*` method.

## Numerical practice

- **Tolerances are named constants, not literals.** `.TOL_ROUNDTRIP <- 1e-12`
  at the top of the file that uses it.
- **Never widen a tolerance without a same-line comment** giving the reason.
- Prefer `solve(A, b)` to `solve(A) %*% b`.
- Prefer `crossprod(X)` to `t(X) %*% X`.
- Use `chol2inv(qr.R(qr(X)))` for cross-product inverses in covariance code.
- Guard every `solve()` that could be singular with `tryCatch(..., error =
  function(e) NULL)` and propagate the `NULL`.
- Check `is.finite()` on anything that will be handed to an optimiser.

## Documentation

Roxygen on every exported function. Required tags: `@param` for all arguments,
`@return` describing the class and its meaningful fields, `@details` where the
mathematics needs it, `@examples` that run in under five seconds.

Where an argument has a non-obvious interaction, say so in its own `@param`
rather than only in `@details`. The `var_lags` / GMM instrument dating
interaction is the canonical example — a user reading only `@param var_lags`
must learn that it constrains valid instruments.

Cross-reference the specification: `@seealso docs/02-math-spec.md section 5`.

## Commits

```
<area>: <imperative summary>          # 60 chars

<body: what changed and why. Reference the roadmap ID if applicable.>

Refs: R-3
```

Areas: `algebra`, `factorise`, `expectations`, `estimate`, `methods`,
`report`, `diagnostics`, `tests`, `docs`, `build`.

Any commit that changes a frozen regression value must say which invariant or
route changed and why the new number is more correct. Commits that touch `R/`
must have `devtools::check()` passing.

## Pull request checklist

- [ ] `devtools::check()` clean, no new NOTEs
- [ ] invariant tests untouched, or the change explicitly justified
- [ ] new numerical route has an independent-route test
- [ ] roxygen regenerated (`devtools::document()`)
- [ ] NEWS.md entry for anything user-visible
- [ ] no new dependency in `R/`
- [ ] generated-regressor caveat still printed by `summary()`
