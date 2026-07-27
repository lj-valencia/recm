# 01 — Architecture

## Design principle

The package is a pipeline from **cost parameters** to **observable residuals**,
wrapped in an optimiser. Every stage is a pure function of the stage before it.
That structure is deliberate: it means any stage can be tested in isolation
against an independent implementation, and it means the optimiser only ever
sees a low-dimensional parameter vector.

```
theta  (2 or m free parameters, unconstrained scale)
  |
  |  .k_from_theta()                       [estimate.R]
  v
k      (m+1 cost weights, k_0 normalised to 1)
  |
  |  .lq_alpha()   Riccati                 [factorise.R]
  |  factorise_spectral()  cross-check, m <= 5 only
  v
alpha  (A(L) coefficients, roots outside unit circle)
  |
  +---> .alpha_to_a()      -> a_0, a_1..a_{m-1}      [algebra.R]
  |
  +---> .scalars()         -> A(1), A(beta), G, sum_d
  |
  +---> .dweights()        -> d_0..d_H scalar leads
  |
  +---> zmech$z(alpha, scalars)                      [expectations.R]
        |     built ONCE, before the optimiser
        v
      Z  (length T), the forward sum. Under .zmech_var() this is
         h' z_{t-1} with h from .hvec(); the residual function does
         not know which mechanism produced it.
        |
        v
  residual  e_t = dy_t - [a_0*ecm + sum a_i dy_{t-i} + a_f*Z + delta'W]
        |
        v
  objective  SSR  or  GMM criterion                  [estimate.R]
```

## Why the auxiliary VAR is estimated once, outside the loop

`H` is a nuisance parameter. Re-estimating it inside the optimiser would make
the criterion jagged (the VAR coefficients would jitter with each trial
`theta`) and `which.min` would pick numerical noise rather than a minimum.
It is estimated once in `recm_estimate()` and held fixed.

The cost is that reported standard errors condition on `H`. That is the
generated-regressor problem, and it is handled by `recm_boot()` rather than by
moving the VAR inside the loop.

## Why linear parameters are concentrated out

`delta` (extra regressors) and, when `free_forward = TRUE`, `a_f` enter the
residual linearly given `alpha`. They are recovered by `lm.fit()` inside the
residual function rather than passed to the optimiser. This keeps Nelder-Mead
in 2 dimensions (geometric costs) instead of `2 + p` dimensions, which matters
because Nelder-Mead degrades badly above roughly 10 parameters.

The full parameter vector is reassembled afterwards for the covariance
calculation, which differentiates the residual with respect to *all*
parameters including the concentrated ones. Concentration is an optimisation
device, not an inference shortcut.

## Two solution routes, deliberately kept

| Route | File | Use |
|---|---|---|
| Spectral factorisation | `factorise.R::factorise_spectral()` | m <= 5, cross-check, teaching |
| Riccati / discounted LQ | `factorise.R::.lq_alpha()` | production, all m |

They solve the same problem. The spectral route roots `z^m P(q(z))`, a
degree-2m polynomial; the Riccati route solves the equivalent LQ regulator and
reads `alpha` off the feedback gain. Keeping both is what makes invariant #6
in `CLAUDE.md` testable.

Do not delete the spectral route because it is "slower" or "redundant". It is
the independent implementation that certifies the fast one.

## Module boundaries

**`algebra.R`** knows about `alpha`, `beta`, `m`. It knows nothing about data,
VARs, or estimation. Pure functions of small numeric vectors. Everything here
should be testable without loading a dataset.

**`factorise.R`** maps `k -> alpha`. Also pure. Does not know what `alpha` will
be used for.

**`expectations.R`** owns the auxiliary VAR and the companion form. This is
the only module that knows the state vector layout — that the leading element
is a constant, that variable `j` at lag `l` sits at position `1 + j + l*k`.
If that layout ever changes, it changes here and nowhere else. `.var_states()`
is that layout, split out of `.var_companion()` so a caller holding an
already-estimated `H` can build states on new observations without
re-estimating anything.

**`estimate.R`** is the only module that touches `data` at estimation time.
It assembles the design, calls the pipeline, and runs the optimiser.

**`predict.R`** touches `data` too, and is the one place that is allowed to:
it rebuilds the same design on new observations. It builds nothing itself —
the state layout comes from `.var_states()`, the forward sum from
`.zmech_var()` at the **fitted** `H`, and the specification from
`recm_estimate()`'s frame via `.fit_env()`. It had its own copy of all three
and disagreed with `.var_companion()` on the first.

**`methods.R` / `report.R`** are presentation only. They must not recompute
anything that changes the answer; they may recompute derived quantities
(delta-method Jacobians, lead weights at new horizons) from stored inputs.

## State vector layout

Fixed by `.var_companion()`:

```
z_t = ( 1, X_t, X_{t-1}, ..., X_{t-p+1} )'      length n_z = 1 + k*p
```

The leading constant carries the VAR intercept, so **no demeaning is
required** anywhere in the package. `H` has `H[1,1] = 1`, which puts a unit
eigenvalue in the companion matrix. This is expected. It means:

- `rho(H)` as reported by the raw eigenvalue is always 1; report `rho_dyn`
  (excluding the constant block) for display.
- The convergence condition `rho(G) * rho(H) < 1` therefore reduces to
  `rho(G) < 1` in practice, which is the correct conservative check.

`d(ystar)` is **always** column 1 of `X`, hence position 2 of `z`. Several
call sites pass `sel = 2L` on that assumption. If you ever allow the target
to move, introduce a named lookup rather than changing the constant.

## Object model

`recm_estimate()` returns an object of class `recmfit` containing:

- the optimisation-scale `theta` and its covariance
- `alpha`, `a`, `k`, `scalars` at the optimum
- `alpha_fn`, a **closure** mapping any `theta` to `alpha`
- `objective`, a closure giving the criterion at any `theta`
- residuals, fitted values, the retained sample index
- the VAR object and GMM diagnostics

The two closures are what make delta-method standard errors and
`recm_profile()` possible without re-specifying the model. They capture the
data by reference, so a `recmfit` object is not portable across sessions
without the data. That is a deliberate trade; document it, do not "fix" it by
serialising the design matrices.

The same closures keep `recm_estimate()`'s evaluation frame alive, and
`recm_boot()` and `predict()` read `data` and the resolved specification back
out of it through `.fit_env()` rather than from a copy on the object — a copy
could go stale against the closures, which are what actually gets evaluated.
Reading it *after* `subset` is the point. `.fit_env()` names what it needs
and fails loudly when the frame does not have it; both callers used to reach
in by name having checked only `data`.

## Extension points

Adding a new **cost parameterisation**: one entry in the `.COST` registry in
`estimate.R`, plus the new name in the `cost` argument of `recm_estimate()`.
Nothing else in the package branches on `cost`.

An entry supplies seven fields, all functions of `m`:

| Field | Returns |
|---|---|
| `npar(m)` | number of free parameters at order `m` |
| `k(th, m)` | cost vector, length `m+1`, with `k[1] == 1` |
| `names(m)` | labels on the **optimisation** scale (these reach `summary()`) |
| `natural(th, m)` | named vector on the interpretable scale |
| `start(m)` | default starting value |
| `grid(m)` | `npar(m)`-column matrix of fallback starts, tried in order when `start(m)` is inadmissible |
| `warn(m)` | `character(1)` if the order is a bad idea, else `NULL` |

The first four are reached through the thin wrappers `.k_from_theta()`,
`.theta_names()`, `.theta_natural()` and `.theta_npar()`; those names are
kept because they read better at the call sites, not because they hold any
logic.

`npar(m)` is the contract: `names()`, `natural()` and `start()` must all have
that length and `grid()` must have that many columns. `test-cost.R` loops the
registry and asserts it at several `m`, so **a new entry is checked without
anyone writing a new test** — if you add one and that file starts failing,
the entry is inconsistent.

Both edit sites fail loudly if you do only one. A name in the signature but
not the registry reaches `.cost_spec()` and stops with the registered set in
the message; a name in the registry but not the signature is rejected by
`match.arg()`. Neither can silently estimate the wrong model — which is what
the previous structure did, three functions each shaped
`if (cost == "geometric") ... else ...`, where any unrecognised name fell
into the `else` and was quietly treated as `"free"`.

Verified end to end by adding a third parameterisation with `npar = 1`
(constant costs, a count equal to neither `2` nor `m`) and touching only
those two sites: NLS and GMM both fit, and `summary()`, `coef()`, `vcov()`,
`confint()`, `equation()`, `lead_weights()` and `recm_profile()` all worked
unchanged.

One caveat found doing that: at `npar = 1`, `stats::optim()` warns that
one-dimensional Nelder-Mead is unreliable. The fit is fine on this fixture,
but a one-parameter parameterisation should probably route to `"Brent"`, and
nothing does that today.

Adding a new **estimator**: the residual function is shared, and the rest is
a documented contract rather than a registry. Three edits.

1. The name in `recm_estimate()`'s `method` argument.
2. If it consumes the automatic instruments, add it to `.METHODS_IV`. That
   gate used to be a bare `method == "gmm"` at the instrument block, which
   made it an invisible second edit site — a new estimator reached the block,
   got `Ziv = NULL`, and died inside `apply()` with `dim(X) must have a
   positive length`.
3. A branch in the estimator block that sets **all seven** of:

| Name | What |
|---|---|
| `th` | `theta` at the optimum |
| `r` | `resid_full(th)`; carries `$b`, `$lin`, `$e` |
| `par_all` | `c(th, r$lin)`, the full parameter vector |
| `V` | covariance, `length(par_all)` square |
| `fit` | optimiser result; only `$convergence` is read |
| `extras` | estimator-specific diagnostics, possibly an empty list |
| `objective` | closure giving the criterion at any `theta` |

This section previously named only `theta`, `V` and `extras`. The other four
are not optional — the object assembly at the end of `recm_estimate()` reads
every one of them, so a branch written to the old description produced a
half-built `recmfit`.

A branch that is declared but not written now hits an explicit `stop()`
naming the seven fields. That guard is unreachable from user code, because
`match.arg()` rejects an unknown name first; it exists for the developer
mid-edit and should not be deleted as dead code.

`extras` is duck-typed on purpose: `print.summary.recmfit()` tests each field
for `NULL` rather than branching on `method`, so a new estimator can print
its own diagnostics without touching `methods.R`.

**Why a contract and not a `.COST`-style registry.** The two estimators share
the residual function and almost nothing else, and the asymmetry is
structural rather than incidental: GMM must build its instruments *before*
sample selection, because the longer lags shorten the estimation sample —
which is upstream of anywhere a plug-in estimator function could run. A
registry would have to hoist that, and there is no third estimator on the
roadmap to justify the restructure. Revisit if one lands; R-6
(Anderson-Rubin) inverts the existing GMM criterion and does not count.

Adding a new **expectations mechanism** (e.g. model-consistent rather than
VAR): one constructor in `expectations.R` returning a mechanism. The residual
function takes `Z` as an input and no longer builds it, so nothing in
`estimate.R` needs to change to add one.

A mechanism is three things:

| Field | What |
|---|---|
| `name` | `character(1)`, reported by `summary()` |
| `support` | logical, length `T`. Rows it can supply `Z` for at **any** admissible `alpha` |
| `z(alpha, s)` | numeric length `T`, or `NULL` if this `alpha` is inadmissible **for this mechanism**. `s` is `.scalars(alpha, beta)`, passed in rather than recomputed |

**`support` must not depend on `theta`.** This is the load-bearing part and
the reason the seam is a list rather than a bare function. The estimation
sample is fixed once, before the optimiser runs, so that every trial `theta`
is scored on the same observations — a criterion compared across moving
samples is not a criterion. `recm_estimate()` reads `support` into the
`complete.cases()` that builds `ok`, where the whole `Slag` matrix used to
sit.

The VAR mechanism satisfies this for free: its `NA` pattern comes from the
lag structure, not from `alpha`. **Perfect foresight does not**, because its
horizon follows `rho(G)` — which is why `.zmech_pf()` takes a *fixed*
horizon and refuses any `alpha` needing more, rather than quietly shortening
the sum. That trade is the open question in R-4: a fixed horizon buys a
fixed sample and costs admissible parameter space.

The admissibility condition belongs to the mechanism, not to the residual
function. The VAR route needs `rho(G) rho(H) < 1` (INVARIANT 8); perfect
foresight carries no `H` and needs `rho(G) < 1` plus a horizon inside
budget. The `rho(G) rho(H) >= 1` test used to sit in the residual builder,
which is why it read as universal when it never was.

Two mechanisms exist. `.zmech_var()` is the production one;
`.zmech_pf()` is not reachable from `recm_estimate()` yet — there is no
`expectations_backend` argument, because how the horizon should be chosen is
undecided. It is exercised through the contract test, which is what keeps
this seam from being notional.
