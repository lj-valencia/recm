# 05 — Roadmap

Work items in rough priority order. Each has an ID used in commit messages and
in `# TODO(R-n)` comments.

Status: `open` / `in progress` / `blocked` / `done`.

---

## R-1 — Propagate first-stage uncertainty in `ystar`  *(open, high)*

`ystar` normally comes from a cointegrating regression estimated outside the
package. Nothing downstream accounts for its sampling error. Superconsistency
makes this second-order asymptotically, but at T around 170 quarters it is not
zero, and it feeds *both* the error-correction term and the forward term
through the VAR.

**Approach.** Accept an optional `ystar_fit` argument (an `lm` or a
user-supplied covariance for the cointegrating vector). Inside `recm_boot()`,
redraw the cointegrating coefficients, rebuild `ystar`, re-estimate the VAR,
and re-run the whole pipeline per replication.

**Acceptance.** On `dgp_pac` with a deliberately noisy first stage, bootstrap
standard errors must cover the Monte Carlo sampling standard deviation, where
the current ones do not.

**Watch for.** Cost. This makes every replication a full re-estimation. Cache
the VAR lag matrices; do not recompute the design.

---

## R-2 — Optimiser bounds from admissibility  *(open, high)*

`m = 2` has closed-form admissibility constraints (docs/02 section 9) that are
currently checked only after the fact, by returning `NULL` from the builder.
The optimiser therefore wanders into inadmissible regions and wastes
evaluations, and Nelder-Mead handles the resulting `1e12` penalty badly.

**Approach.** Map the constraints into the unconstrained parameterisation and
supply `lower`/`upper` to `nlminb`, keeping Nelder-Mead as a fallback. For
`m > 2` there is no closed form — use a penalty that is continuous at the
boundary rather than a cliff.

**Acceptance.** No `1e12` returns in a normal fit; convergence in fewer
function evaluations on the fixtures; identical optima.

---

## R-3 — Correct cumulative bands on lead weights  *(open, medium)*

`plot.recmfit()` draws the cumulative-lead band by summing pointwise standard
errors. That is conservative and wrong in the same way that adding standard
deviations instead of variances is wrong.

**Approach.** Accumulate the Jacobian rows before sandwiching:
`J_cum <- apply(J, 2, cumsum)` then `sqrt(diag(J_cum V J_cum'))`. The Jacobian
is already computed in `lead_weights()`; expose it as an attribute.

**Acceptance.** Band width strictly narrower than the current one at every
horizon, and matching a bootstrap of the cumulative sum.

---

## R-4 — Model-consistent expectations  *(open, medium)*

Only VAR-based expectations are supported. FRB/US also runs a
model-consistent (MCE) mode where expectations come from solving the full
model forward.

**Approach.** This is the refactor flagged at the end of docs/01: the residual
function should take a `Z` series as an input rather than building it, once a
second mechanism exists. Introduce an `expectations_backend` argument with
`"var"` and `"perfect"` first — perfect foresight is a one-line implementation
and a genuine cross-check, since certainty equivalence means the analytic rule
under perfect foresight must reproduce the exact quadratic-program solution
(verified to 3.4e-07 over 40 periods during development).

**Do this refactor when the second backend lands, not before.**

**Status.** The perfect-foresight mechanism has landed as `.zpf()` in
`expectations.R`, standalone — nothing in `estimate.R` calls it yet, so the
refactor above is still to do, but its precondition is now met.

It was not a one-line implementation, for one reason worth knowing before
the wiring is written: **the sum cannot be collapsed.** The VAR route
telescopes the infinite forward sum with the Kronecker identity and needs no
horizon; perfect foresight sums over realised data, so it must truncate, and
the last `H` observations have no future to sum over. `.zpf()` returns those
as `NA` rather than padding them.

`H` is not small. At the reference calibration (`.ref_k`, `beta = 0.98`,
`rho(G) = 0.817`) it is **114** at a relative remainder of `1e-10`, and 223
for a sluggish `m = 2`, 716 for `m = 3` at `beta = 0.999`. On 178 quarterly
observations that leaves 64. **Perfect foresight is a cross-check to run on
long simulated samples, not a backend to offer on typical macro data**, and
`expectations_backend = "perfect"` should say so when it is exposed.

Verified against three routes that share nothing with it but `.dweights()`:
a constant target (`Z = c * sum_d`), a geometric target (`Z_t = dys_t *
sum_i d_i rho^i`, the closed form the DGP fixture uses), and the VAR closure
itself on a path a VAR(2) forecasts exactly — agreement 1.1e-10, which is
the truncation tolerance and nothing else. That last one also pins the `t-1`
dating, the easiest thing to get wrong when the second call site is written.

Still to do here: the residual-function refactor, the
`expectations_backend` argument, and the certainty-equivalence check against
a direct quadratic-program solution that the note above describes.

---

## R-5 — Weak-identification reporting as a first-class output  *(open, medium)*

`recm_profile()` is documented as a diagnostic, but given how routinely the
criterion is flat, it should be part of standard output.

**Approach.** Add `profile = TRUE` to `summary.recmfit()`, printing the profile
interval alongside each delta-method interval, and flagging when the two
disagree materially. Consider making a flat profile a `warning()` rather than
a printed note.

**Acceptance.** On the `m = 20` fixture, `summary()` surfaces the flatness
without the user having to know to ask.

---

## R-6 — Anderson-Rubin intervals for GMM  *(open, low)*

Weak instruments are the norm here, not the exception — first-stage F was 3.3
at `var_lags = 4` on the reference dataset. Wald intervals are unreliable
under weak identification.

**Approach.** Implement identification-robust Anderson-Rubin confidence sets by
inverting the GMM criterion over a grid of `theta`. Two dimensions under
`cost = "geometric"` makes this tractable.

**Blocked on** R-2, since it needs a reliable admissible region to grid over.

---

## R-7 — Support for data with gaps  *(open, low)*

Currently `data` must be regularly spaced with no gaps, and **this is not
checked**. A gap silently corrupts every lag.

**Approach.** Minimum viable fix is detection: require a date column or an
explicit `frequency`, verify regular spacing, and `stop()` otherwise. Actually
handling gaps is a much larger change and probably not worth it.

**Acceptance.** A test that gapped data errors rather than returning numbers.

---

## R-8 — Package the reference DGP as exported data  *(open, low)*

`data-raw/dgp.R` is currently a script. Ship `dgp_ecm` and `dgp_pac` as
documented datasets so examples and vignettes do not regenerate them.

---

## R-9 — Vignette: "From cointegrating regression to PAC equation"  *(open, low)*

End-to-end worked example: Engle-Granger first stage, VAR specification, both
estimators, both reporting formats, the `free_forward` restriction test, and
an honest section on what the sample can and cannot identify.

**Depends on** R-8.

---

## R-10 — Recover the efficiency lost to instrument dating  *(open, medium)*

`iv_lag` now defaults to `var_lags + 2`, so the automatic instrument set
obeys the dating rule and GMM corrects measurement error in `ystar` as it is
supposed to. It is paid for in variance: on the reference DGP with no
measurement error, RMSE on `a_0` is 0.0541 at the default against 0.0212 at
`iv_lag = 1`, a factor of 2.6. The default is still right — `iv_lag = 1`
under measurement error is *inconsistent*, bias +0.087 at `T = 3000` and not
shrinking — but the trade should not be this expensive.

**Approach.** Three candidates, in order of expected payoff. (a) Add the
exogenous expectations variables at more lags, since they are exempt from
the rule and are where the surviving identification comes from. (b) Use a
continuously-updated or iterated GMM weight matrix; the two-step `W2` is
estimated on 5-7 moments at `T ~ 170` and is probably the binding
constraint. (c) Drop the weakest instrument columns by first-stage
contribution rather than keeping all of them.

**Acceptance.** RMSE at the default dating within 1.5x of `iv_lag = 1` on
the clean DGP, with the measurement-error correction of the current default
preserved — both arms of the table in `?recm_estimate` must be re-run.

**Watch for.** The exemption is load-bearing. Lagging the whole state vector
leaves the moments valid but underidentified: the criterion goes numerically
flat, about 4e-05 at the true parameters against 0.0 at a `psi` boundary
point, and `a_0` stalls near -0.11 at any `T`. Any reshuffle of the
instrument set has to keep the exogenous columns at `t-1`.

---

## R-11 — Calibrate Hansen J  *(open, medium)*

On the correctly specified `dgp_pac`, Hansen J returns `p = 0.0032`. A
correctly specified model should not be rejected at that rate, so either the
statistic, the weight matrix it is built from, or the fixture is wrong, and
until that is known J should not be relied on. The R-10 dating fix moved it
from 2.8e-05 to 0.0032 without curing it, which points at the weight matrix
rather than at instrument validity.

**Approach.** Monte Carlo the null distribution of `extras$J` on `dgp_pac`
and compare against `chi2(Jdf)`. Suspect the HAC weight matrix first:
`.nw()` is estimated on few moments and short samples, and `solve()` of a
near-singular estimate inflates the quadratic form. Compare against an
unweighted J and against a bootstrap null.

**Acceptance.** Empirical size within Monte Carlo error of nominal at 5% on
the correctly specified fixture, and power against a fixture with a
deliberately omitted lead term.

---

## Explicitly out of scope

Recorded so they do not get re-proposed.

- **Estimating `beta`.** Not identified jointly with the costs. Calibrate.
- **Free `k_j` for large `m`.** Variance amplification makes it meaningless.
  `cost = "geometric"` exists for this.
- **Replacing the spectral factorisation with the Riccati route entirely.**
  The slow route is the independent check that certifies the fast one.
- **Suppressing the generated-regressor caveat.** It is always true.
- **A `dplyr`/`data.table` dependency.** Base R and `stats` only in `R/`.
