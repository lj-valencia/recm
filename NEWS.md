# recm (development version)

## 0.0.0.9000

First working implementation. Every routine named in
`docs/03-api-reference.md` is now implemented and tested; the previous
release was scaffolding only.

### Estimation

* `recm_estimate()` fits the equation by NLS or two-step Newey-West GMM,
  with `cost = "geometric"` or `"free"`, optional `free_forward`, `growth`,
  `vars`, `subset` and Nelder-Mead `restarts`. `delta` and, when free,
  `a_f` are concentrated out via `lm.fit()`.
* `recm_profile()`, `recm_boot()`, `equation()` and `lead_weights()`, plus
  `print`/`summary`/`coef`/`vcov`/`confint`/`logLik`/`residuals`/`fitted`/
  `plot` methods for `recmfit` and `print`/`plot` for `recm_profile`.
* `confint()` and `logLik()` are additions beyond `docs/03-api-reference.md`.

### GMM instruments

* **The automatic instrument set now obeys the package's own dating rule.**
  It was `cbind(1, ecm, DYL, Slag[, -1])`, with `ecm` at `t-1` and
  `d(ystar)` at `t-1 .. t-p` — all inside the `t-(p+2)` window the composite
  error contaminates when `ystar` carries measurement error. New `iv_lag`
  argument, defaulting to `var_lags + 2`, dates the `y`/`ystar`-derived
  columns correctly. `iv_lag = 1` restores the old set, is more efficient,
  and warns.
* **Only the `y`/`ystar`-derived columns are lagged; exogenous expectations
  variables stay at `t-1`.** This exemption is what makes the fix work.
  Lagging the whole state vector leaves the moments valid but
  underidentified — the criterion goes numerically flat (about 4e-05 at the
  true parameters against 0.0 at a `psi` boundary point) and `a_0` stalls
  near -0.11 at any `T`.
* Measured on the reference DGP, RMSE for `a_0`. With `ystar` measured
  exactly the dating costs a factor of 2.6 (0.0541 against 0.0212 at
  `T = 178`). With `ystar` badly mismeasured it is the difference between a
  consistent estimator and an inconsistent one: at `T = 3000` NLS and
  `iv_lag = 1` both sit at a bias of +0.087 that does not shrink, while the
  new default gives -0.009. The full table is in `?recm_estimate`.
* Validity checked directly: at the default dating the residual is
  orthogonal to every instrument column, `p > 0.05` on all of them at
  `T = 3000`.
* **The intercept is no longer dropped from the instrument matrix.** The
  zero-variance filter, which exists to remove degenerate columns, was
  deleting the one column that is meant to be constant, silently removing
  the `E[e_t] = 0` moment.
* **The first-stage F no longer regresses the error-correction term on
  itself.** With the intercept restored, `ecm` appeared on both sides at
  `iv_lag = 1` and F came back as 1e31. It is now computed on the excluded
  instruments only, and is `NA` when there are none.
* `extras` gained `iv_lag`; `summary()` reports the dating and how many
  instrument columns survived pruning.
* Hansen J on the correctly specified fixture improved from `p = 2.8e-05` to
  `p = 0.0032`, but still rejects. Roadmap R-11.

### Cost parameterisations

* **The extension point in `docs/01-architecture.md` now works as
  advertised.** It claimed a new cost parameterisation was three functions
  and that "nothing downstream changes"; both halves were wrong.
* Dispatch was `if (cost == "geometric") ... else ...` in `.k_from_theta()`,
  `.theta_names()` and `.theta_natural()`, so **any name that was not
  `"geometric"` was silently treated as `"free"`** — `.k_from_theta()`
  returned `c(1, exp(th))`, a perfectly plausible cost vector, and
  estimation ran to completion on the wrong model. Dispatch is now a lookup
  in a `.COST` registry that `stop()`s on an unregistered name.
* Starting values and the fallback restart grid were written the same way at
  their call sites, and the `cost = "free"` order warning was a third such
  branch, so *three* things downstream did change. All three now come from
  the registry entry. A parameterisation whose parameter count is neither
  `2` nor `m` previously got a starting vector of length `m`, which
  `.k_from_theta()` then read positionally without error.
* `npar(m)` makes the parameter count explicit; `test-cost.R` loops the
  registry and asserts `names()`, `natural()`, `start()` and `grid()` all
  agree with it, so a parameterisation added later is checked without a new
  test being written.
* **A wrong-length `start` now errors** instead of being read positionally.
* `cost = "free"` is now exercised end to end. The suite previously tested
  only that it *warns* above `m = 4`, never that it fits.
* No behaviour change for `"geometric"` or `"free"`: every frozen value in
  `test-regression.R` is unmoved.

### Numerics

* **The Riccati iteration now converges on the feedback gain, not on `P`.**
  At `m = 20` with slowly decaying costs `max|P|` reaches 1e8, so a
  `P`-scaled threshold is really a scale-dependent one and the iteration
  stalls short of it — `.lq_alpha()` returned `converged = FALSE`, which
  made every `m = 20` fit fail outright. The gain is the quantity wanted
  (it *is* `alpha`) and stays O(1).
* `.TOL_RICCATI` is 1e-10 **on the gain**. It must not be loosened to make
  slow cases certify: value iteration converges linearly at a rate
  approaching 1 for high `m`, so the per-step change understates the
  remaining error by roughly `1/(1-rate)`. Measured at `m = 20`,
  `psi = 0.79`, a per-step change of 1.1e-07 sat 4.5e-05 from the converged
  gain — outside INVARIANT 6. Such cases are rejected, not approximated.
* `factorise_spectral()` written from scratch (it was absent from the
  source file) and returns the normalisation constant `cc` and the rescaled
  `k_norm` alongside `alpha`, `roots` and a `reciprocity` conditioning
  diagnostic. It warns above `m = 5`.
* `.a_to_alpha()` written from scratch; round-trips against `.alpha_to_a()`
  to 5.6e-17.
* `.hvec()` now returns `NULL` when `rho(G) rho(H) >= 1` rather than
  relying on its caller to check (INVARIANT 8).
* `.scalars()` returns `NULL` on non-finite `alpha` instead of erroring
  inside `eigen()`.
* `.dweights()` returns both `h` and `d`; `.scalars()` gained `sum_h`.
* `.polymul()` replaces `convolve()` for polynomial multiplication, to keep
  FFT noise out of the already ill-conditioned spectral route.
* `.getcols()` no longer silently drops incomplete rows. `model.matrix()`
  defaults to `na.omit`, so a supplied `vars` or `instruments` column with
  any `NA` returned a matrix shorter than `data` and misaligned against
  every other series. Any lagged column has leading `NA`s — and the dating
  rule tells users to supply exactly that — so this made the documented
  workflow impossible. It now passes `NA` through and errors if the row
  count still does not match.

### INVARIANT 3 restated

**INVARIANT 3 did not hold as written, and the documentation has been
corrected.** `k_0 = A(1)A(beta)` requires the factorisation constant
`cc = 1`, which is a *different* normalisation from the `k_0 = 1` the
package estimates in; the two cannot both be imposed. Matching leading
coefficients gives `cc = k_m * prod(zeta_i)`; evaluating `z^m P(q(z))` at
`z = 1`, where `q(1) = 0`, gives `k_0 = cc * A(1)A(beta)`. For the
documented fixture `k = c(1, 27.6486, 3.2239)`, `beta = 0.98`: `k[1] = 1`
but `A(1)A(beta) = 0.02481524` and `cc = 40.29781`.

The invariant is now stated as `k_norm[1] = A(1)A(beta)`, with
`k_norm = k / cc`, verified to 4.5e-16. The same rescaling applies to the
`m = 2` identity, now stated as `k_2 / cc = a_1`. Amended in
`docs/02-math-spec.md` sections 3 and 9, `docs/04-testing.md`, and the
invariant list and `m = 2` trap in `CLAUDE.md`. No code changed —
`factorise_spectral()` already returned `cc` and `k_norm`, and the tests
already asserted the corrected form.

### Known discrepancy with the documentation

Flagged rather than silently reconciled.

* **`recm_boot()` is a recursive bootstrap, not a fixed-design one.**
  `docs/03-api-reference.md` describes the equation step as fixed-design; a
  fixed design is not available, because the error-correction regressor is
  built from `y`, so perturbing `y` necessarily moves the design.

### Testing

* 435 passing expectations, no skips. All eight invariant tests from
  `docs/04-testing.md` pass, INVARIANT 3 in the restated form above. The
  count is up from 196 largely because `test-cost.R` loops the `.COST`
  registry across several `m` rather than naming each parameterisation.
* Fixtures `simulate_dgp_pac()` / `simulate_dgp_ecm()` in
  `tests/testthat/helper-recm.R`, with the target's forward sum in closed
  form so the fixture satisfies the Euler equation exactly rather than
  approximately. Mirrored by `data-raw/dgp.R`; R-8 will ship them as data.
* Growth neutrality verified: predicted balanced-growth gap −0.001438
  against −0.001636 realised, negative as `docs/02` section 8 requires.
* Consistency verified at `T = 4000` (`a_0` −0.1535 against a true
  −0.1550). At `T = 178` there is a genuine finite-sample downward bias in
  `|a_0|` of 0.51–0.67 sd, measured at `var_lags` 1, 2 and 4 alike, so it
  is not VAR over-parameterisation. `test-estimate.R` therefore asserts
  consistency and a 1.0 sd bias bound rather than the 0.4 sd used as an
  illustration in `docs/04-testing.md`.

### Not claimed

* That the GMM instrument set is *efficient*. It is now valid, which it was
  not, and the cost is a factor of 2.6 in RMSE when there is no measurement
  error to correct. Roadmap R-10.
* That Hansen J is calibrated. On the correctly specified `dgp_pac` it
  returns `p = 0.0032`, which needs explaining before the statistic is
  relied on. Roadmap R-11.
* That `recm_boot()` is calibrated (roadmap R-1). Its test asserts shape
  and that it runs, nothing more.
