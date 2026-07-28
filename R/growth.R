## Growth neutrality.
##
## On a balanced growth path the target grows at a constant rate g and a
## growth neutral equation must deliver y = ystar rather than a permanent
## level gap. Substituting dy = g and dystar = g into
##
##   dy[t] = a0 (ystar[t-1] - y[t-1]) + sum_i a_i dy[t-i] + Z[t]
##
## gives a0 * gap = g * (1 - sum_i a_i - sum_i d_i), so the condition is
##
##   R(a) = 1 - sum_{i=1}^{m-1} a_i - sum_i d_i(a; beta) = 0.
##
## FRB/US and Dynare satisfy this by adding the correction R(a) * g to the
## equation. This package instead imposes R(a) = 0 directly as a nonlinear
## cross-coefficient restriction on the estimated {a}, so nothing is added to
## the equation and the balanced growth path is neutral by construction.
##
## Note that at beta = 1 the operator A(beta F) A(L) is symmetric in L, a path
## linear in t passes through it exactly, and R(a) vanishes identically for
## every a and every m. The restriction therefore only binds when beta < 1.

# The growth neutrality gap R(a).
gn_gap <- function(a, beta) {
  alg <- pac_algebra(a, beta)
  gn_gap_alg(alg)
}

gn_gap_alg <- function(alg) {
  ar_sum <- if (alg$m > 1L) sum(alg$a[-1L]) else 0
  1 - ar_sum - alg$d_sum
}

# Gradient of R with respect to a, by central differences. R is smooth in a
# wherever the model is admissible, so a fixed step is adequate.
gn_grad <- function(a, beta, step = 1e-6) {
  vapply(
    seq_along(a),
    function(j) {
      up <- a
      dn <- a
      up[j] <- up[j] + step
      dn[j] <- dn[j] - step
      (gn_gap(up, beta) - gn_gap(dn, beta)) / (2 * step)
    },
    numeric(1)
  )
}

# The permanent level gap y - ystar implied by a non-neutral equation on a
# balanced growth path with drift g.
gn_implied_gap <- function(a, beta, g) {
  -gn_gap(a, beta) * g / a[1L]
}
