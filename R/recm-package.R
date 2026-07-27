#' @keywords internal
#' @aliases recm-package
#'
#' @details
#' `recm` estimates polynomial adjustment cost (PAC) equations under a
#' rational error correction (REC) framework. A PAC equation is the solution
#' to the intertemporal problem
#'
#' \deqn{\min E_{t-1} \sum_i \beta^i \left[ k_0 (y - y^*)^2 +
#'   \sum_{j=1}^{m} k_j ((1-L)^j y)^2 \right]}
#'
#' The user supplies a decision variable `y`, a frictionless target `ystar`
#' and a **calibrated** discount factor `beta`; the package recovers the
#' adjustment cost weights `k` and reports the implied decision rule.
#'
#' The entry point is [recm_estimate()]. Because the criterion is routinely
#' near-flat, [recm_profile()] and [recm_boot()] should be treated as part of
#' the standard output rather than as optional diagnostics.
#'
#' @section Three things to know before using this package:
#'
#' **Calibrate `beta`.** It is not identified jointly with the cost
#' parameters; it enters only through the composite operator
#' \eqn{(1 - \beta F)(1 - L)}. Use `1/(1+r)` with a real quarterly rate. The
#' default of `1` is admissible and warns.
#'
#' **The criterion is often nearly flat.** This is inverse optimal control:
#' many cost configurations rationalise nearly identical behaviour. Point
#' estimates alone are not a defensible output.
#'
#' **Reported standard errors are optimistic.** `Z_t` is a generated
#' regressor depending on the auxiliary VAR and on `ystar`. In Monte Carlo
#' the conditional standard errors came in at roughly half the true sampling
#' standard deviation. [recm_boot()] is the fix.
#'
#' @seealso `docs/02-math-spec.md` for the mathematics this package
#'   implements.
"_PACKAGE"
