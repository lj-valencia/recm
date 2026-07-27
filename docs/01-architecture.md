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
  +---> .hvec(alpha, beta, H_var, sel)               [expectations.R]
        |
        v
      h  (length n_z), so Z_t = h' z_{t-1}
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
If that layout ever changes, it changes here and nowhere else.

**`estimate.R`** is the only module that touches `data`. It assembles the
design, calls the pipeline, and runs the optimiser.

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

## Extension points

Adding a new **cost parameterisation**: implement in `.k_from_theta()`,
`.theta_names()`, `.theta_natural()`. Nothing downstream changes.

Adding a new **estimator**: add a branch in `recm_estimate()` that produces
`theta`, `V`, and `extras`. The residual function is shared.

Adding a new **expectations mechanism** (e.g. model-consistent rather than
VAR): implement an alternative to `.hvec()` returning a `Z_t` series. The
residual function should take `Z` as an input rather than building it, once a
second mechanism exists. Do this refactor when the second mechanism lands, not
before.
