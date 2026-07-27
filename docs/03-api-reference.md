# 03 — API reference

Exported surface. Changing any signature here is a breaking change and needs a
NEWS.md entry.

---

## `recm_estimate()`

```r
recm_estimate(y, ystar, vars = NULL, beta = 1, data,
             m = 3, cost = c("geometric", "free"),
             expectations = NULL, var_lags = 4, diff_exp = TRUE,
             method = c("nls", "gmm"), instruments = NULL, iv_lag = NULL,
             free_forward = FALSE, growth = NULL,
             hac_lags = NULL, start = NULL, subset = NULL,
             maxit = 4000, restarts = 3, quiet = FALSE)
```

### Core arguments

| Argument | Type | Notes |
|---|---|---|
| `y` | character(1) | column name, **in levels** (usually logs). Estimated in first differences. |
| `ystar` | character(1) | frictionless target, same units as `y`. Supplied by the user, typically from a cointegrating regression. Its estimation error is **not** propagated. |
| `vars` | formula or character | extra regressors, entering linearly with free coefficients. Concentrated out of the optimiser. |
| `beta` | numeric(1) | discount factor. Default `1`; warns. See below. |
| `data` | data.frame | ordered in time, **no gaps**. Not checked — gaps will silently corrupt lags. |

### Model structure

| Argument | Default | Notes |
|---|---|---|
| `m` | `3` | cost polynomial order; equation carries `m-1` lags of `dy`. |
| `cost` | `"geometric"` | `k_j = kappa*psi^(j-1)`, 2 params. `"free"` gives `m` params — only usable for `m <= 4`, and warns above it regardless of `quiet`. Adding a third is one `.COST` entry plus one word here; see the extension point in `docs/01`. |
| `free_forward` | `FALSE` | if `TRUE`, `a_f` estimated freely instead of restricted to `sum(d_i)`. This is the Euler test. |
| `growth` | `NULL` | column name of trend growth; adds the growth-neutrality correction. |

### Expectations

| Argument | Default | Notes |
|---|---|---|
| `expectations` | `NULL` | extra VAR variables. `d(ystar)` always included and always occupies position 1. |
| `var_lags` | `4` | **also sets GMM instrument dating** through `iv_lag`: instruments from `y`/`ystar` must be dated `t-(var_lags+2)` or earlier. |
| `diff_exp` | `TRUE` | difference the `expectations` variables before entering the VAR. |

### Estimation

| Argument | Default | Notes |
|---|---|---|
| `method` | `"nls"` | or `"gmm"` (two-step, Newey-West). |
| `instruments` | `NULL` | *excluded* instruments for GMM, added to the automatic set. Taken as supplied — `iv_lag` is not applied to them. Collinear columns dropped by QR pivot. |
| `iv_lag` | `NULL` → `var_lags + 2` | lag for the `y`/`ystar`-derived automatic instruments (ecm, lagged `dy`, `d(ystar)` states). Other state variables stay at `t-1` regardless — that exemption is what keeps the model identified. `iv_lag = 1` is more efficient and valid only if `ystar` is measured exactly; it warns. |
| `hac_lags` | `NULL` | Newey-West bandwidth; `NULL` uses `floor(4*(T/100)^(2/9))`. |
| `start`, `subset`, `maxit`, `restarts`, `quiet` | | optimiser and sample control. `restarts` re-runs Nelder-Mead from its own solution, which matters — a single pass routinely stops short. `start` must have one element per free cost parameter and errors otherwise. |

### On `beta = 1`

Admissible: the forward sum still converges provided `rho(G) < 1`. But it
means indifference between the present and the arbitrarily distant future, and
it lengthens the effective lead horizon. `beta = 1/(1+r)` with a **real**
quarterly rate is the defensible choice — the loss is quadratic in log
deviations, which are unit-free, so the real rate is the relevant one even
when `y` and `ystar` are nominal. A warning fires on the default path unless
`quiet = TRUE`.

### Value

Object of class `recmfit`. Fields: `theta`, `par`, `vcov`, `alpha`, `a`, `k`,
`scalars`, `a_f`, `delta`, `residuals`, `fitted`, `dep`, `index`, `n`, `npar`,
`var`, `extras`, `converged`, plus two closures:

- `alpha_fn(theta)` — any `theta` to `alpha`; powers delta-method SEs
- `objective(theta)` — the criterion; powers `recm_profile()`

**The closures capture `data` by reference.** A `recmfit` is not portable
across sessions without its data. Deliberate.

---

## `summary()` / `print()`

`summary(fit)` prints both reporting formats in one call:

1. estimated parameters on the optimisation scale, `printCoefmat` style
2. cost parameters on the natural scale, delta method
3. **Format 1 — explicit**: `a_0`, `a_1..a_{m-1}`, `delta`, forward coefficient fixed at 1
4. **Format 2 — compressed**: `a_0`, `sum a_i`, `a_f`, mean lead, half-life, plus the `f_i` table
5. fit and diagnostics: residual SE, R-squared, Ljung-Box, Durbin-Watson, `rho(G)*rho(H)`, cost admissibility, Hansen J, first-stage F

The generated-regressor caveat is printed unconditionally. **Do not add a
switch to suppress it.**

`print(fit)` gives a one-paragraph digest.

---

## `equation()`

```r
equation(object, digits = 4, format = c("both","explicit","compressed"),
         max_lags = 12)
```

Writes the fitted equation with coefficients substituted, in either or both
formats. Returns `list(a, a_f, d, f)` invisibly. `max_lags` truncates the
printed lag list; the count of suppressed lags is shown.

---

## `predict()`

```r
predict(object, newdata = NULL,
        interval = c("none", "confidence", "prediction"),
        level = 0.95, type = c("diff", "level"),
        boot_object = NULL, quiet = FALSE, ...)
```

Obtains predictions and optional confidence or prediction intervals from a fitted `recmfit` object, matching standard `predict.lm()` conventions.

- `newdata`: optional `data.frame` containing out-of-sample observations.
- `type`: `"diff"` (predicts \eqn{\Delta y_t}) or `"level"` (predicts \eqn{y_t = y_{t-1} + \Delta y_t}).
- `interval`: `"none"`, `"confidence"`, or `"prediction"`.
- `boot_object`: optional `recm_boot` object from `recm_boot()`. When supplied, standard errors and intervals are computed across bootstrap replicates. When `NULL`, delta-method standard errors are used.

---

## `lead_weights()`

```r
lead_weights(object, horizon = 24, level = 0.95, normalised = FALSE)
```

Data frame `(i, weight, se, lo, hi)`. `normalised = TRUE` returns `f_i`
instead of `d_i`. Standard errors are delta method: numerically differentiate
`d_i(theta)` and sandwich the estimated covariance. No re-estimation.

**Known limitation:** the cumulative band drawn by `plot()` sums pointwise
standard errors rather than accumulating the Jacobian first, so it is
conservative. Fixing this is roadmap item R-3.

---

## `recm_profile()`

```r
recm_profile(object, which = 1, span = 3, ngrid = 25, level = 0.95)
```

Profiles the criterion in one parameter, re-optimising the others at each grid
point. Returns `grid`, `value`, `a0`, `cutoff`, `interval`, `rel_range`.

`rel_range < 0.02` prints a weak-identification warning. Given how routinely
that fires, **`recm_profile()` should be treated as part of the standard output,
not an optional diagnostic.**

---

## `recm_boot()`

```r
recm_boot(object, R = 299, var_error = TRUE, eq_error = TRUE, seed = NULL)
```

Bootstrap standard errors that propagate VAR estimation error:

- recursive residual bootstrap of the auxiliary VAR, giving `H*`
- fixed-design wild bootstrap of the equation error

Does **not** propagate uncertainty in `ystar`. Superconsistency makes that
second-order asymptotically; at T around 170 it is not zero. Roadmap item R-1.

In Monte Carlo, `recm_boot()` reproduced the true sampling standard deviations
closely, while the conditional standard errors were roughly half of them.

---

## Internal helpers

Not exported; stable within the package but may change without notice.

| Function | File | Purpose |
|---|---|---|
| `.alpha_to_a`, `.a_to_alpha` | algebra.R | exact inverse maps |
| `.G_mat`, `.scalars`, `.dweights` | algebra.R | companion, totals, scalar leads |
| `.lq_alpha` | factorise.R | Riccati; production route |
| `factorise_spectral` | factorise.R | polynomial route; cross-check only, m <= 5 |
| `.var_companion` | expectations.R | VAR, constant-augmented companion, state matrix |
| `.hvec` | expectations.R | closed-form `h` so `Z_t = h'z_{t-1}` |
| `.jacnum`, `.nw`, `.lagv`, `.getcols` | aaa-utils.R | numerical Jacobian, Newey-West, lags, model matrix |
| `.halflife` | algebra.R | interpolated gap half-life (must stay differentiable) |

`.halflife()` interpolates the 0.5 crossing linearly. An integer step would
have zero numerical derivative and the delta-method SE would come back as
exactly 0 — which is what happened before the fix. Keep it continuous.
