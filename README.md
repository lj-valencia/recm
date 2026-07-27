# recm

Estimation of **Rational Error Correction (REC)** models in R — the
forward-looking error-correction specification used in the Federal Reserve
Board's FRB/US model and the Bank of Canada's LENS.

Base R and `stats` only.

## What it does

A REC equation is the solution to

```
min  E sum_i beta^i [ (y - y*)^2 + sum_j k_j ((1-L)^j y)^2 ]
```

Supply a decision variable, a frictionless target, and a discount factor;
`recm` recovers the adjustment cost parameters and reports the implied
decision rule in both the explicit and the compressed (LENS-style) formats.

```r
fit <- recm_estimate(
  y     = "cons_nd",           # decision variable, in log levels
  ystar = "cons_nd_star",      # frictionless target, from your first stage
  vars  = ~ rel_price,         # extra regressors
  beta  = 0.995,               # calibrate this; do not estimate it
  data  = quarterly,
  m     = 3,
  expectations = ~ hh_income,
  var_lags = 4,
  method = "gmm"
)

summary(fit)       # both reporting formats, fit, diagnostics
equation(fit)      # the fitted equation, coefficients substituted
plot(fit)          # lead weights, cumulative, lag distribution, gap response
recm_profile(fit)  # run this. the criterion is usually flatter than you expect
```

## Installation

```r
# install.packages("devtools")
devtools::install_github("lj-valencia/recm")
```

## Documentation

| | |
|---|---|
| `docs/01-architecture.md` | pipeline, module boundaries, object model |
| `docs/02-math-spec.md` | the mathematics the code implements |
| `docs/03-api-reference.md` | exported surface |
| `docs/04-testing.md` | invariant tests and statistical testing rules |
| `docs/05-roadmap.md` | open work |
| `docs/06-conventions.md` | style, errors, numerics, commits |

## Three things to know before using it

**Calibrate `beta`.** It is not identified jointly with the cost parameters —
it enters only through `(1 - beta F)(1 - L)`. Use `1/(1+r)` with a real
quarterly rate. The default of `1` is legal and warns.

**The criterion is often nearly flat.** This is inverse optimal control: many
cost configurations rationalise nearly identical behaviour. Point estimates
alone are not a defensible output. Use `recm_profile()` and `recm_boot()`.

**Reported standard errors are optimistic.** `Z_t` is a generated regressor
depending on the auxiliary VAR and on `ystar`. In Monte Carlo the conditional
standard errors came in at roughly half the true sampling standard deviation.

## Status

Pre-release, but working: the estimator, both reporting formats and the
diagnostics are implemented, with 435 passing expectations and a clean
`devtools::check()`. Consistency is verified against a synthetic REC DGP at
`T = 4000`; at `T = 178` there is a genuine finite-sample downward bias in
`|a_0|` of roughly half a standard deviation.

Two things are **not** claimed, and are worth reading `NEWS.md` for before
relying on them: Hansen J rejects a correctly specified fixture at
`p = 0.0032`, so treat the overidentification test as uncalibrated, and
`recm_boot()` is not calibrated either.

GMM instruments are dated `t-(var_lags+2)` by default, which is what makes
the estimator robust to measurement error in `ystar` — the usual case, since
`ystar` normally comes from a first-stage regression. It costs precision when
`ystar` is measured exactly; `iv_lag = 1` buys that back and is valid only
then.

The API is not stable. See `docs/05-roadmap.md`.
