## Both imports below are load-bearing and neither is decorative. Registering
## an S3 method requires the generic to be resolvable as the namespace loads,
## and R only finds it without an import when the generic is one of base R's
## .knownS3Generics - which print, plot, predict, summary, coef and the rest
## are, but nobs and simulate are not. Dropping either name here does not fail
## at lint or at build; it fails at load, as three R CMD check WARNINGs
## reading "object 'simulate' not found whilst loading namespace 'recm'".

#' recm: Rational Error Correction Models
#'
#' Estimation of rational error correction (REC) models, the econometric
#' representation of polynomial adjustment cost (PAC) behaviour described in
#' Tinsley (2002) and in the FRB/US "PAC basics" note.
#'
#' A rational error correction model is not written down, it is derived from an
#' intertemporal quadratic optimisation problem under rational expectations.
#' The derivation buys a set of cross-equation restrictions linking the speed of
#' adjustment, the weights on expected future changes in the target, and the
#' structural adjustment costs. `recm()` imposes those restrictions.
#'
#' @section Entry point:
#' [recm()] is the only exported estimator. Everything else in the package is
#' internal machinery reachable through the fitted object.
#'
#' @references
#' Tinsley, P. A. (2002). Rational error correction. *Computational
#' Economics*, 19(2), 197--225.
#'
#' Brayton, F., Davis, M., and Tulip, P. (2000). Polynomial adjustment costs in
#' FRB/US. Federal Reserve Board.
#'
#' @keywords internal
#'
#' @importFrom stats nobs simulate
"_PACKAGE"
