# recm

> [!WARNING]
> **This package is experimental.** It is at version 0.1.0, it is not on CRAN,
> and it has not been validated against a reference implementation on real
> data. The API is unstable: argument names, the contents of the fitted object,
> and the printed output may all change without deprecation. The reported
> standard errors are known to understate sampling uncertainty (see
> [Caveats](#caveats)). Do not use it for published work without checking the
> numbers yourself.

Estimation of **rational error correction (REC)** models in R — the
econometric representation of polynomial adjustment cost (PAC) behaviour
described in Tinsley (2002) and in the FRB/US *PAC basics* note.

A rational error correction model is not written down, it is *derived* from an
intertemporal quadratic optimisation problem under rational expectations. The
derivation buys a set of cross-equation restrictions linking the speed of
adjustment, the weights on expected future changes in the target, and the
structural adjustment costs. `recm()` imposes those restrictions rather than
assuming them:

$$\Delta y_t = a_0 (y^*_{t-1} - y_{t-1}) + \sum_{i=1}^{m-1} a_i \Delta y_{t-i} + Z_t + \delta' W_t + \varepsilon_t$$

where $Z_t = \sum_{i \ge 0} d_i E[\Delta y^*_{t+i}]$ is the forward term. The
lead weights $d_i$, the speed of adjustment $a_0$ and the autoregressive
coefficients $a_i$ are all functions of the same lag polynomial $A(L)$, and
imposing that dependence is what separates this from an error correction model
written down by hand. The coefficient on $Z_t$ is fixed at 1 by construction.

## Installation

The package is written in base R and depends only on packages that ship with
R (`stats`, `utils`, `graphics`, `grDevices`). There is nothing to compile.

```r
# install.packages("remotes")
remotes::install_github("lj-valencia/recm")
```

Requires R >= 3.5.0.

### From a local clone

The project uses [renv](https://rstudio.github.io/renv/) for its development
library, which activates automatically from `.Rprofile`.

```sh
git clone https://github.com/lj-valencia/recm.git
cd recm
```

```r
renv::restore()                                  # development library
devtools::install_deps(dependencies = TRUE)
devtools::load_all()                             # load without installing
```

## Quick start

The example below builds a mock data generating process *from the decision
rule itself*, so there is a known truth to compare the estimates against.

### 1. Simulate

```r
library(recm)

set.seed(42)

n   <- 240
a0  <- 0.3   # true speed of adjustment
phi <- 0.6   # persistence of the target's growth rate
mu  <- 0.5   # mean growth rate of the target

# The frictionless target: growth is a persistent AR(1) around mu.
dystar <- numeric(n)
dystar[1] <- mu
for (t in 2:n) {
  dystar[t] <- mu + phi * (dystar[t - 1] - mu) + rnorm(1, 0, 0.4)
}
ystar <- 100 + cumsum(dystar)

# The forward term. At m = 1 and discount = 1 the lead weights are
# d_i = a0 (1 - a0)^i, which sum to one, and expectations dated t-1 of an
# AR(1) target collapse the infinite forward sum to a closed form.
z <- mu + a0 * phi * (c(NA, dystar[-n]) - mu) / (1 - (1 - a0) * phi)

# The decision variable.
y <- numeric(n)
y[1] <- ystar[1]
for (t in 2:n) {
  y[t] <- y[t - 1] + a0 * (ystar[t - 1] - y[t - 1]) + z[t] + rnorm(1, 0, 0.15)
}

dat <- ts(cbind(y = y, ystar = ystar), start = c(1965, 1), frequency = 4)
```

### 2. Estimate

`y` and `y_star` are bare column names of `data`. Every *other* numeric column
of `data` is taken as an exogenous regressor $W_t$ — see
[Gotchas](#gotchas).

```r
fit <- recm(y, ystar, data = dat)
fit
```

```
Rational error correction model

Call:
  recm(y = y, y_star = ystar, data = dat)

dy[t] = 0.3099 * (ystar[t-1] - y[t-1]) + Z[t]

  m = 1 (0 lags of the dependent variable)   discount = 1
  expectations: var   estimator: iterative OLS (17 iterations)
  observations: 238   residual std. error: 0.151
```

The estimated speed of adjustment is 0.3099 against a true 0.30, and the
implied lead weights track $a_0(1-a_0)^i$:

```r
round(head(fit$lead_weights, 6), 4)   # 0.3099 0.2139 0.1476 0.1018 0.0703 0.0485
round(a0 * (1 - a0)^(0:5), 4)         # 0.3000 0.2100 0.1470 0.1029 0.0720 0.0504
```

### 3. Extract the summary

`summary()` prints the coefficient table and then everything the
cross-equation restrictions imply — the roots of $A(L)$, the total forward
loading, the mean lead and half-life, the structural adjustment costs, the
growth neutrality status, and the auxiliary autoregression.

```r
summary(fit)
```

```
Coefficients:
   Estimate Std. Error t value Pr(>|t|)
ec 0.309905   0.008471   36.59   <2e-16
Standard errors: just-identified GMM, HAC lag 4, including dZ/da'.
...

Fit, dy against its fitted value:
  R-squared (uncentered)  0.9384   adjusted  0.9381
...

Implied dynamics:
  A(L) roots (modulus)  1.449
  total forward loading sum(d_i)  1
  mean lead  2.227 quarters   half-life  1.869 periods
  adjustment costs k  0.09604, 0.69009   (k0 = A(1)A(beta), residual 2.2e-16)

Growth neutrality:
  slack at discount = 1: R(a) = 0.00e+00 for any coefficients
  target drift 0.4595, implied permanent level gap 0

Forward loading, restricted against free (levels 2 vs 3):
  restricted 1 (s.e. 4.702e-13)   free 1.06 (s.e. 0.02022)
  The two standard errors omit different terms, so no test statistic
  is formed; compare the intervals by eye.

Auxiliary model: AR(1) in dystar with a constant (order by AIC)
  constant 0.196   ar 0.5735   rho(Phi) 1
```

The usual extractors work as expected: `coef()`, `vcov()`, `fitted()`,
`residuals()`, `nobs()`, `logLik()`.

### 4. Forecast

`predict()` is a **dynamic simulation** of the decision rule, not a one step
ahead conditional prediction. `y` is iterated forward from the end of the
estimation sample, so `newdata` supplies the future path of the *target* and
of any exogenous regressors only. A `y` column in `newdata` is ignored.

With `n_ahead` alone, the target is projected from the fitted auxiliary
autoregression:

```r
fc <- predict(fit, n_ahead = 12)
fc
```

```
Forecast from a rational error correction model

  horizon: 12   origin: 2024.75   expectations: var
  target ystar: projected from the auxiliary autoregression

Contributions to dy (the term columns sum to fit):
 horizon time       ec forward    fit
       1 2025 0.148689  0.5343 0.6830
       2 2025 0.124633  0.5024 0.6270
       3 2026 0.098639  0.4841 0.5827
       4 2026 0.075313  0.4736 0.5489
...
```

With `newdata`, you supply the path yourself:

```r
future <- ts(cbind(ystar = ystar[n] + cumsum(rep(0.5, 8))),
             start = c(2025, 1), frequency = 4)

fc2 <- predict(fit, newdata = future)
head(fc2$levels, 4)
```

```
  horizon    time    base        ec   forward      fit
1       1 2025.00 210.874 0.1486888 0.5343011 211.5570
2       2 2025.25 210.874 0.2406681 1.0056970 212.1203
3       3 2025.50 210.874 0.3130071 1.4770930 212.6641
4       4 2025.75 210.874 0.3717924 1.9484890 213.1943
```

The forecast is decomposed term by term, exactly:

```r
all.equal(rowSums(fc$contributions[, fc$terms]), unname(fc$fit))  # TRUE
```

`$contributions` decomposes $\Delta y$; `$levels` cumulates the same
decomposition and carries the forecast origin $y_T$ in its own `base` column.
These are contributions to the *simulated path*, not partial derivatives — the
error correction column carries the feedback from every term before it.

### 5. Plot

```r
op <- par(mfrow = c(2, 2))
plot(fit)                          # panels 1:4 by default, 1:9 available
par(op)

plot(fit, which = c(3, 4, 9), ask = FALSE)   # lead weights and the ec response

plot(fc)                           # stacked contributions to the forecast
plot(fc, type = "level")           # the same, cumulated
```

## Key arguments

| Argument | Default | What it does |
| --- | --- | --- |
| `m` | `1` | Order of the adjustment cost polynomial. Fixes the number of lags of $\Delta y$ at `m - 1`. `m = 1` is the canonical REC model. Unrelated to the length of the forward sum, which is never truncated. |
| `discount` | `1` | The discount factor $\beta$, calibrated and never estimated. Must lie in $(0, 1]$. |
| `expectations` | `"var"` | `"var"` uses an auxiliary univariate autoregression with information dated $t-1$, giving $Z_t$ in closed form. `"mce"` uses model consistent expectations on the realised target path. |
| `method` | `"ols"` | Only the FRB/US iterative OLS zig-zag is implemented. `"nls"` and `"gmm"` are accepted by the signature and error informatively. |
| `tr_exog` | `TRUE` | Enters exogenous regressors as $\Delta W_t$. Set `FALSE` only when the regressor is genuinely stationary in levels (a dummy, a spread, a gap). It applies to all of them at once. |

### Growth neutrality

On a balanced growth path the equation delivers $y = y^*$ only if
$R(a) = 1 - \sum_i a_i - \sum_i d_i = 0$. FRB/US and Dynare satisfy this by
*adding* $R(a) g$ to the equation. `recm()` instead imposes $R(a) = 0$ as a
nonlinear restriction on the estimated coefficients, so nothing is added and
neutrality holds by construction.

- At `discount = 1` (the default) the restriction holds identically and is
  not imposed — it is reported as slack.
- At `discount < 1` with `m = 1` it degenerates ($a_0 = 1$, instantaneous
  adjustment), so that combination is refused with an error. Use `m >= 2`.

## Gotchas

- **Every remaining numeric column of `data` becomes an exogenous regressor.**
  A `data.frame` with a numeric `time = 1:n` column will silently estimate a
  coefficient on `d_time`. Use a `ts` object, or a `Date`/`character` column
  as the time index — at most one non-numeric column is allowed and it is
  taken as the index rather than as a regressor.
- **No intercept is fitted.** A free constant is inconsistent with growth
  neutrality; the auxiliary autoregression carries the drift instead. This is
  why `r.squared` is **uncentered** — the total sum of squares is taken about
  zero rather than about $\overline{\Delta y}$, so the drift in the target
  counts as explained variation. It is not comparable with an R-squared from
  an equation carrying a constant.
- **`y_star` is given, not estimated.** The frictionless target is an input.
- `y` and `y_star` must be complete: the forward term and the auxiliary
  autoregression are built from contiguous differences.

## Caveats

Beyond the experimental status of the package as a whole:

- **The standard errors understate sampling uncertainty.** They are the
  just-identified GMM sandwich, including the numerical $\partial Z/\partial a'$
  that the final zig-zag sweep's OLS standard errors would omit — but they
  condition on the *estimated* auxiliary model. In simulation, holding the
  target path fixed across replications, they matched the sampling standard
  deviation to within a tenth; resampling the target process as well, they
  came in at roughly a **third** of it. Bootstrap the auxiliary model before
  believing an interval.
- **No forecast intervals are reported**, deliberately. An interval built on
  the standard errors above would inherit the problem and be too narrow by
  about as much.
- **Only iterative OLS is implemented.** `"nls"` and `"gmm"` are planned.
- **No test statistic is formed** for the restricted-versus-free forward
  loading comparison. The two standard errors omit different terms and their
  covariance is not available in closed form, so a t-ratio on the difference
  would be arbitrary. Compare the intervals by eye.

## Development

```sh
Rscript -e "devtools::install_deps(dependencies = TRUE)"
Rscript -e "devtools::test()"                # run the testthat suite
Rscript -e "devtools::document()"            # rebuild roxygen docs + NAMESPACE
Rscript -e "devtools::load_all(); lintr::lint_package()"
Rscript -e "styler::style_pkg()"
R CMD build .
R CMD check --as-cran .
```

`lint_package()` needs the namespace loaded first, otherwise it reports a
large number of phantom `object_usage` warnings.

Testing-only packages belong in `Suggests`, never `Imports`. Core estimation
logic stays in base R.

## References

- Tinsley, P. A. (2002). Rational error correction. *Computational Economics*,
  19(2), 197–225.
- Brayton, F., Davis, M., and Tulip, P. (2000). *Polynomial adjustment costs in
  FRB/US*. Federal Reserve Board.

## License

GPL-3.
