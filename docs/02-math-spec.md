# 02 — Mathematical specification

Everything in `algebra.R` and `factorise.R` implements this document. If code
and this document disagree, one of them is a bug — stop and resolve it rather
than picking whichever seems more likely.

Notation: `L` lag operator, `F` forward operator, `D = (1-L)`.

---

## 1. The optimisation problem

The agent chooses a path for `y` to minimise

```
    E_{t-1} sum_{i>=0} beta^i [ k_0 (y_{t+i} - ystar_{t+i})^2
                                + sum_{j=1}^{m} k_j ( D^j y_{t+i} )^2 ]
```

taking `y_{t-1}, ..., y_{t-m}` as given and `{ystar}` as an exogenous
stochastic process. All `k_j > 0`. Only ratios are identified, so `k_0 = 1`
by normalisation throughout the package.

`j = 1` penalises velocity, `j = 2` acceleration, and so on. The order `m`
fixes the number of lags of `dy` in the decision rule at `m - 1`.

---

## 2. Euler equation

Differentiating with respect to `y_s` and collecting own-period, one-ahead and
further-ahead terms, every penalty collapses onto powers of a single composite
operator:

```
    Q = (1 - beta F)(1 - L)

    [ k_0 + k_1 Q + k_2 Q^2 + ... + k_m Q^m ] y_t  =  k_0 ystar_t
```

Write `P(q) = sum_j k_j q^j`. This form is why the penalties are written as
powers of `(1-L)` rather than as arbitrary lag polynomials: each `k_j` term
differentiates into exactly `k_j * ((1-beta F)(1-L))^j`.

---

## 3. Factorisation

Setting `L = z`, `F = 1/z`, so `q(z) = (1 - beta/z)(1 - z)`, there is a unique
factorisation

```
    P(q(z)) = cc * A(beta/z) A(z),     A(z) = 1 + alpha_1 z + ... + alpha_m z^m
```

with all roots of `A` outside the unit circle. Equivalently

```
    A(beta F) A(L) y_t = A(1) A(beta) ystar_t
```

`cc` is not free. `A` is pinned by its unit leading coefficient and by the
root-selection rule above, so matching the leading coefficients of `z^m P(q(z))`
against `cc * z^m A(beta/z)A(z)` fixes

```
    cc = k_m * prod_i zeta_i
```

**`cc = 1` and `k_0 = 1` are different normalisations and cannot both hold.**
Section 1 fixes `k_0 = 1`, and that is what the package estimates in, so `cc`
is whatever the factorisation returns. Evaluating `z^m P(q(z))` at `z = 1`,
where `q(1) = 0` and only the `k_0` term survives, gives the relation that
holds under *either* normalisation:

```
    k_0 = cc * A(1) A(beta)

    equivalently, with k_norm = k / cc,

    k_norm[1] = A(1) A(beta)                                    [INVARIANT 3]
```

`factorise_spectral()` returns `cc` and `k_norm` for this reason; it is
`k_norm[1]`, not `k[1]`, that the invariant test compares against
`A(1)A(beta)`. For the fixture `k = c(1, 27.6486, 3.2239)`, `beta = 0.98`:
`cc = 40.29781` and `A(1)A(beta) = 0.02481524`, against `k[1] = 1`.

This remains the cheapest available check that a factorisation succeeded — it
is one polynomial evaluation — but it checks a ratio, not `k_0` itself.

### 3a. Spectral route (m <= 5)

Multiply by `z^m`. Then `z^m P(q(z))` is an ordinary polynomial of degree `2m`
whose roots come in pairs `(zeta, beta/zeta)` — exactly `m` outside the unit
circle and `m` inside. Keep the outside set and build

```
    A(z) = prod_i (1 - z / zeta_i)
```

Implementation note: `z*q(z) = -beta + (1+beta) z - z^2`, so the coefficient
vector of `z^m P(q(z))` is assembled by convolving that quadratic `j` times and
shifting by `m - j`.

**This route is ill-conditioned for larger m.** At m = 20 the coefficients
span roughly 5e12 and the reciprocity check `|zeta_in * zeta_out / beta - 1|`
degrades to about 1e-3. Warn above m = 5; do not use it in production.

### 3b. Riccati route (all m)

The same problem is a discounted linear-quadratic regulator. With state
`x_t = (y_{t-1}, ..., y_{t-m})'` and control `u_t = y_t`:

```
    x_{t+1} = A x_t + B u_t          A = shift-down, B = e_1
```

Each penalty `D^j y_t = sum_{l=0}^{j} (-1)^l C(j,l) y_{t-l}` is a linear
function of `(u_t, x_t)`, so the stage cost is
`v' M v` with `v = (u, x)` and `M = sum_j k_j g_j g_j'`. Partition `M` into
`R` (on `u`), `N` (cross), `Q` (on `x`), and iterate

```
    S  = R + beta B' P B
    K  = N + beta A' P B
    P <- Q + beta A' P A - K K' / S
```

to convergence. The feedback gain

```
    F = ( N' + beta B' P A ) / S
```

satisfies `u_t = -F x_t`, i.e. `y_t = -sum_i F_i y_{t-i}`, so

```
    alpha = F                                                   [INVARIANT 6]
```

Closed-loop eigenvalues of `A - B F` are the reciprocals of the roots of `A(z)`
and must all lie strictly inside the unit circle.

---

## 4. Reduced form

From `A(L) y_t = A(1) y_t - sum_{k=0}^{m-1} c_k D y_{t-k}` with
`c_k = sum_{i>k} alpha_i`:

```
    a_0 = -A(1) = -(1 + sum_i alpha_i)                          [INVARIANT 2]
    a_i = sum_{j>i} alpha_j          i = 1..m-1
```

The inverse map, used to go from a target decision rule back to costs:

```
    c_0 = -a_0 - 1
    c_i = a_i                                    i = 1..m-1
    alpha_i = c_{i-1} - c_i                      i = 1..m-1
    alpha_m = c_{m-1}
```

These are exact inverses.                                       [INVARIANT 1]

Note `a_0 < 0` for a stable equation, hence `A(1) > 0`. The FRB/US "PAC
Basics" note prints `a_0 = A(1)` without the minus sign; that is inconsistent
with its own equations (4) and (5).

---

## 5. Forward weights

With `G` the companion matrix of `A(beta F)` (bottom row
`(-alpha_m beta^m, ..., -alpha_1 beta)`) and `iota = e_m`:

```
    h_i = A(1) A(beta) * iota' G^i iota                    weights on ystar
    d_i = A(1) A(beta) * iota' (I-G)^{-1} G^i iota         weights on d.ystar
```

The scalar `h_i` equal `A(1)A(beta)` times the `i`-th coefficient of the power
series `1 / A(beta F)`, obtainable by polynomial long division. [INVARIANT 5]

Totals in closed form:

```
    sum_i h_i = A(1) A(beta) * iota' (I-G)^{-1}   iota
    sum_i d_i = A(1) A(beta) * iota' (I-G)^{-2}   iota
```

---

## 6. VAR closure

With the auxiliary VAR in companion form `z_t = H z_{t-1}` and `e_v` selecting
the target growth rate from the state,

```
    Z_t = sum_{i>=0} d_i E_{t-1}[ D ystar_{t+i} ]
        = A(1) A(beta) ( iota'(I-G)^{-1} (x) e_v' ) [I - G (x) H]^{-1}
          ( iota (x) H ) z_{t-1}
        = h' z_{t-1}
```

using `sum_i G^i (x) H^i = [I - G (x) H]^{-1}`. **No horizon truncation is
required.**                                                     [INVARIANT 4]

Convergence requires `rho(G) rho(H) < 1`.                       [INVARIANT 8]

---

## 7. Decision rule, both reporting formats

**Explicit:**

```
    D y_t = a_0 (y - ystar)_{t-1} + sum_{i=1}^{m-1} a_i D y_{t-i}
            + sum_{i>=0} d_i E_{t-1}[D ystar_{t+i}]  + delta' W_t + e_t
```

The coefficient on the forward sum is 1 by construction.

**Compressed (LENS / FRB-US presentation):**

```
    D y_t = a_0 (y - ystar)_{t-1} + sum_{i=1}^{m-1} a_i D y_{t-i}
            + a_f sum_{i>=0} f_i E_{t-1}[D ystar_{t+i}] + delta' W_t + e_t

    f_i = d_i / sum_j d_j,   sum_i f_i = 1                      [INVARIANT 7]
    a_f = sum_j d_j          under the Euler restriction
```

Setting `free_forward = TRUE` estimates `a_f` freely. Comparing it to
`sum d_i` evaluated at the estimated `alpha` is the direct test of the Euler
restriction, and is the single most informative diagnostic the package offers.

---

## 8. Growth neutrality

On a balanced growth path with `D ystar = g`, the exact Euler solution
satisfies

```
    y - ystar = g ( 1 - sum_{i>=1} a_i - sum_i d_i ) / a_0   <  0
```

The optimum runs *below* a growing target — the agent optimally lags, because
moving is costly. This is correct behaviour.

`growth = "<column>"` adds `(1 - sum a_i - sum d_i) * g_t` to the forward
term, which forces `E[y - ystar] = 0`. This is a deliberate departure from the
exact optimum and matches FRB/US practice. Document it as such wherever it
appears; do not present the corrected rule as "the solution".

---

## 9. Admissibility

Beyond root stability, cost non-negativity restricts the reduced form. At
`m = 2` the constraints are sharp:

```
    k_2 / cc = a_1                exactly  ->  a_1 > 0 required
    k_1 > 0  <=>  (1 + a_0 + a_1)(1 + beta a_1) > 2 (1 + beta) a_1
             <=>  a_1 < ( B - sqrt(B^2 - 4 beta (1+a_0)) ) / (2 beta),
                  B = 1 + beta - beta a_0
    k_0 / cc = A(1)A(beta) > 0    automatic when a_0 < 0
```

The two `cc` divisions are the section 3 normalisation and matter only to the
first line as an *equality*: `k_2 = a_1` is false under `k_0 = 1`, `k_2 / cc =
a_1` is exact. The sign conditions are unaffected, `cc > 0`, and the `k_1`
inequality is stated in the reduced-form `a` which are normalisation-free.

A negative coefficient on `D y_{t-1}` is not rationalisable by any quadratic
adjustment cost structure at `m = 2`. These belong in optimiser bounds.

For general `m` there is no closed form; check `all(k > 0)` and report it.

---

## 10. Identification — known and expected

Documented so nobody spends a week rediscovering it.

- **beta is not identified jointly with the costs.** It enters only through
  `q(z)`. Moving beta from 0.95 to 0.999 changed a GMM criterion by less than
  a factor of two while `kappa` moved 10 -> 162. Calibrate it externally.

- **The criterion is often nearly flat.** At `m = 2`, a fourfold change in
  `k_1/k_0` cost 2.2% of SSR; the true parameters sat 1.7% above the optimum.
  Monte Carlo showed the estimator unbiased but with a standard deviation on
  `a_1` roughly equal to `a_1` itself.

- **This is inverse optimal control.** Many cost configurations rationalise
  nearly identical observed behaviour. Report profile regions, not points.

- **`m` is not free.** At `m = 20` on 178 quarterly observations nothing is
  identified even after reducing to two free cost parameters. Consistency was
  confirmed at T = 4000. The limitation is sample size, not the estimator.
