## Graphical diagnostics for fitted rational error correction models.
##
## The panels fall into two groups. One is the usual regression diagnostic set,
## adapted to a time series equation: residuals against fitted values, their
## distribution, their path through the sample, their autocorrelation. The
## other shows the estimated structure itself - the lead weights d_i and the
## deterministic adjustment path - which is where a rational error correction
## model differs from an error correction model written down by hand. The lead
## weights are not free parameters; they are a nonlinear function of {a}
## through A(L), so their shape is the whole content of the cross-equation
## restrictions and is worth looking at directly.

# A restrained palette, kept in one place so the panels agree with one another.
recm_pal <- function() {
  list(pt = "grey30", ref = "grey55", hi = "#B03A2E", alt = "#21618C",
       fill = "grey85")
}

# The reduced-form coefficients a = (a0, a1, ..., a_{m-1}). They are the first
# m entries of the coefficient vector by construction of the design.
recm_a <- function(x) unname(x$coefficients[seq_len(x$m)])

# Where to stop drawing a decaying sequence: the last element that is not
# negligible beside the largest one, floored so a fast decay still gets a
# readable axis.
recm_trim <- function(v, tol = 1e-4, min_h = 8L) {
  scale <- max(abs(v))
  if (!is.finite(scale) || scale <= 0) {
    return(min_h)
  }
  big <- which(abs(v) > tol * scale)
  # v[1] is horizon 0, so the horizon of element i is i - 1.
  max(min_h, if (length(big)) max(big) - 1L else 0L)
}

# Lead weights at an arbitrary horizon, recomputed from the fit rather than
# read off `x$lead_weights`, which stores only the first 41.
recm_weights <- function(x, horizon = NULL, cap = 400L) {
  alg <- pac_algebra(recm_a(x), x$discount)
  if (is.null(horizon)) {
    d <- lead_weights(alg, cap)
    horizon <- min(cap, recm_trim(d))
    d <- d[seq_len(horizon + 1L)]
  } else {
    d <- lead_weights(alg, horizon)
  }
  list(d = d, horizon = horizon, d_sum = alg$d_sum)
}

# The horizontal axis of the time series panels. `index` carries whatever the
# user's `data` supplied: a numeric time, a Date, or period labels. The first
# two plot directly and keep their own axis method; labels are plotted against
# position with the ticks written in by hand.
recm_time <- function(index, n) {
  usable <- length(index) == n &&
    (is.numeric(index) || inherits(index, "Date") || inherits(index, "POSIXt"))
  if (usable) {
    list(x = index, labels = NULL)
  } else if (length(index) == n) {
    list(x = seq_len(n), labels = as.character(index))
  } else {
    list(x = seq_len(n), labels = NULL)
  }
}

# Draw the tick labels for the character-index case.
recm_time_axis <- function(tm, pal) {
  if (is.null(tm$labels)) {
    return(invisible(NULL))
  }
  n <- length(tm$x)
  at <- unique(round(seq(1, n, length.out = min(6L, n))))
  graphics::axis(1, at = at, labels = tm$labels[at], col = pal$ref)
  invisible(NULL)
}

# Name the id_n most extreme points, so an outlier can be traced back to a
# date rather than to a row number.
#
# A label centred above its point runs off the panel when the point is near
# either end of the horizontal axis, which is exactly where the extremes of a
# quantile plot sit, so labels near an edge are pushed inwards instead.
recm_label_extremes <- function(px, py, labels, id_n, pal) {
  if (id_n < 1L || !length(py) || is.null(labels)) {
    return(invisible(NULL))
  }
  ord <- order(abs(py), decreasing = TRUE)[seq_len(min(id_n, length(py)))]
  span <- range(as.numeric(px), finite = TRUE)
  where <- if (diff(span) > 0) {
    (as.numeric(px[ord]) - span[1L]) / diff(span)
  } else {
    rep(0.5, length(ord))
  }
  pos <- ifelse(where > 0.85, 2L, ifelse(where < 0.15, 4L, 3L))
  graphics::text(px[ord], py[ord], labels = labels[ord], pos = pos,
                 cex = 0.7, col = pal$hi, xpd = TRUE)
  invisible(NULL)
}

# Panel 1. Residuals against fitted values, with a lowess smooth. A visible
# slope or bend in the smooth is the sign that the equation is missing
# something the fitted value already knows about.
recm_panel_resid_fit <- function(x, main, id_n, pal, ...) {
  fv <- as.numeric(stats::fitted(x))
  r <- as.numeric(stats::residuals(x))
  plot(fv, r, main = main, pch = 16L, cex = 0.6, col = pal$pt,
       xlab = paste0("Fitted d", x$variables$y), ylab = "Residual", ...)
  graphics::abline(h = 0, lty = 3L, col = pal$ref)
  graphics::lines(stats::lowess(fv, r), col = pal$hi, lwd = 2)
  recm_label_extremes(fv, r, names(stats::residuals(x)), id_n, pal)
  invisible(NULL)
}

# Panel 2. The distribution of the residuals against the normal the reported
# standard errors would imply. The reference curve is centred at zero rather
# than at the residual mean: no intercept is fitted, so the residuals are not
# mechanically centred and a visible shift is information rather than an
# artefact of the plot.
recm_panel_hist <- function(x, main, pal, ...) {
  r <- as.numeric(stats::residuals(x))
  h <- graphics::hist(r, breaks = "FD", plot = FALSE)
  dens <- stats::density(r)
  sigma <- sqrt(x$sigma2)
  grid <- seq(min(h$breaks), max(h$breaks), length.out = 256L)
  norm <- stats::dnorm(grid, mean = 0, sd = sigma)
  plot(h, freq = FALSE, main = main, xlab = "Residual", col = pal$fill,
       border = "white", ylim = c(0, max(h$density, dens$y, norm)), ...)
  graphics::lines(grid, norm, col = pal$hi, lwd = 2)
  graphics::lines(dens, col = pal$alt, lwd = 2, lty = 2L)
  graphics::rug(r, col = pal$ref)
  graphics::abline(v = 0, lty = 3L, col = pal$ref)
  graphics::legend(
    "topright", bty = "n", cex = 0.75, lwd = 2, lty = c(1L, 2L),
    col = c(pal$hi, pal$alt),
    legend = c(paste0("N(0, ", format(sigma, digits = 3), "^2)"),
               "kernel density")
  )
  invisible(NULL)
}

# Panel 3. The lead weights themselves. The centre of mass is the mean lead
# reported by summary(), drawn in so the two agree visibly.
recm_panel_leads <- function(x, main, w, pal, ...) {
  i <- seq_along(w$d) - 1L
  plot(i, w$d, type = "h", lwd = 2, col = pal$pt, main = main,
       xlab = "Horizon i (periods ahead)",
       ylab = expression(paste("lead weight ", d[i])), ...)
  graphics::points(i, w$d, pch = 16L, cex = 0.5, col = pal$pt)
  graphics::abline(h = 0, col = pal$ref)
  if (is.finite(x$mean_lead) && x$mean_lead >= 0 &&
        x$mean_lead <= max(i)) {
    graphics::abline(v = x$mean_lead, lty = 2L, col = pal$hi)
    graphics::text(x$mean_lead, max(w$d), pos = 4L, cex = 0.7, col = pal$hi,
                   labels = paste0("mean lead ", format(x$mean_lead,
                                                        digits = 3)))
  }
  # Inside the panel rather than under the title: the weights decay, so the
  # top right corner is empty, and a subtitle crowds the title in a multi
  # panel layout.
  graphics::legend("topright", bty = "n", cex = 0.75, text.col = pal$ref,
                   legend = paste0("sum(d_i) = ",
                                   format(w$d_sum, digits = 4)))
  invisible(NULL)
}

# Panel 4. The same weights accumulated, normalised by their total so the
# vertical axis reads as the share of the forward response already accrued.
# The normalisation is meaningless if the total is zero, so that case falls
# back to the raw cumulative sum and says so on the axis.
recm_panel_cum_leads <- function(x, main, w, pal, ...) {
  i <- seq_along(w$d) - 1L
  cum <- cumsum(w$d)
  degenerate <- abs(w$d_sum) < sqrt(.Machine$double.eps)
  yy <- if (degenerate) cum else cum / w$d_sum
  ylab <- if (degenerate) {
    "cumulative sum of d_i"
  } else {
    "share of total forward loading"
  }
  plot(i, yy, type = "s", lwd = 2, col = pal$alt, main = main,
       xlab = "Horizon i (periods ahead)", ylab = ylab, ...)
  graphics::abline(h = 0, col = pal$ref)
  if (!degenerate) {
    graphics::abline(h = 1, lty = 3L, col = pal$ref)
    for (share in c(0.5, 0.9)) {
      hit <- which(yy >= share)
      if (!length(hit)) {
        next
      }
      at <- hit[1L]
      graphics::points(i[at], yy[at], pch = 16L, col = pal$hi)
      graphics::text(i[at], yy[at], pos = 4L, cex = 0.7, col = pal$hi,
                     labels = paste0(round(100 * share), "% by i = ", i[at]))
    }
  }
  invisible(NULL)
}

# Panel 5. The residuals through the sample, with the two standard error band
# the reported sigma implies.
recm_panel_resid_time <- function(x, main, id_n, pal, ...) {
  r <- as.numeric(stats::residuals(x))
  tm <- recm_time(x$index, length(r))
  plot(tm$x, r, type = "h", col = pal$pt, main = main,
       xlab = "Time", ylab = "Residual",
       xaxt = if (is.null(tm$labels)) "s" else "n", ...)
  recm_time_axis(tm, pal)
  graphics::abline(h = 0, col = pal$ref)
  graphics::abline(h = c(-2, 2) * sqrt(x$sigma2), lty = 3L, col = pal$hi)
  recm_label_extremes(tm$x, r, names(stats::residuals(x)), id_n, pal)
  invisible(NULL)
}

# Panel 6. Normal quantile plot of the standardised residuals.
recm_panel_qq <- function(x, main, id_n, pal, ...) {
  z <- as.numeric(stats::residuals(x)) / sqrt(x$sigma2)
  qq <- stats::qqnorm(z, plot.it = FALSE)
  plot(qq$x, qq$y, main = main, pch = 16L, cex = 0.6, col = pal$pt,
       xlab = "Theoretical quantile", ylab = "Standardised residual", ...)
  stats::qqline(z, col = pal$hi, lwd = 2)
  recm_label_extremes(qq$x, qq$y, names(stats::residuals(x)), id_n, pal)
  invisible(NULL)
}

# Panel 7. Residual autocorrelation, with the HAC bandwidth the reported
# standard errors actually used marked on the lag axis. Autocorrelation
# reaching past that line is autocorrelation the Bartlett kernel gave no
# weight to.
recm_panel_acf <- function(x, main, lag_max, pal, ...) {
  r <- as.numeric(stats::residuals(x))
  n <- length(r)
  if (is.null(lag_max)) {
    lag_max <- max(1L, min(as.integer(10 * log10(n)), n - 2L))
  }
  ac <- stats::acf(r, lag.max = lag_max, plot = FALSE)
  lags <- as.numeric(ac$lag)[-1L]
  vals <- as.numeric(ac$acf)[-1L]
  band <- stats::qnorm(0.975) / sqrt(n)
  plot(lags, vals, type = "h", lwd = 2, col = pal$pt, main = main,
       xlab = "Lag", ylab = "Residual autocorrelation",
       ylim = 1.2 * range(c(vals, band, -band)), ...)
  graphics::abline(h = 0, col = pal$ref)
  graphics::abline(h = c(-band, band), lty = 2L, col = pal$alt)
  if (!is.na(x$hac_lag) && x$hac_lag >= 1L && x$hac_lag <= max(lags)) {
    graphics::abline(v = x$hac_lag + 0.5, lty = 3L, col = pal$hi)
    graphics::text(x$hac_lag + 0.5, 1.1 * max(abs(vals)), pos = 4L,
                   cex = 0.7, col = pal$hi,
                   labels = paste0("HAC lag ", x$hac_lag))
  }
  invisible(NULL)
}

# Panel 8. The differenced dependent variable against its fitted value. The
# uncentered R-squared reported by summary() is flattered by the drift in the
# target; this panel shows how much of what is left the equation tracks.
recm_panel_actual_fitted <- function(x, main, pal, ...) {
  fv <- as.numeric(stats::fitted(x))
  r <- as.numeric(stats::residuals(x))
  dy <- fv + r
  tm <- recm_time(x$index, length(dy))
  plot(tm$x, dy, type = "l", col = pal$pt, main = main,
       xlab = "Time", ylab = paste0("d", x$variables$y),
       xaxt = if (is.null(tm$labels)) "s" else "n", ...)
  recm_time_axis(tm, pal)
  graphics::lines(tm$x, fv, col = pal$hi, lwd = 1.5)
  graphics::legend("topleft", bty = "n", cex = 0.75, lwd = c(1, 1.5),
                   col = c(pal$pt, pal$hi),
                   legend = c(paste0("d", x$variables$y), "fitted"))
  invisible(NULL)
}

# Panel 9. The deterministic response of the error correction gap to a one-off
# unit displacement, with the target held flat. This is the half-life reported
# by summary(), drawn out in full: at m > 1 the path need not be monotone, and
# a single number cannot say so.
recm_panel_adjustment <- function(x, main, horizon, pal, cap = 200L, ...) {
  a <- recm_a(x)
  if (is.null(horizon)) {
    full <- gap_path(a, cap)
    horizon <- min(cap, recm_trim(full, tol = 1e-2))
    path <- full[seq_len(horizon + 1L)]
  } else {
    path <- gap_path(a, horizon)
  }
  h <- seq_along(path) - 1L
  plot(h, path, type = "l", lwd = 2, col = pal$alt, main = main,
       xlab = "Periods since the displacement",
       ylab = "Remaining gap (unit shock)", ...)
  graphics::abline(h = 0, col = pal$ref)
  graphics::abline(h = 0.5, lty = 3L, col = pal$ref)
  if (is.finite(x$half_life) && x$half_life <= max(h)) {
    graphics::abline(v = x$half_life, lty = 2L, col = pal$hi)
    graphics::text(x$half_life, max(path), pos = 4L, cex = 0.7, col = pal$hi,
                   labels = paste0("half-life ",
                                   format(x$half_life, digits = 3)))
  }
  invisible(NULL)
}

#' Diagnostic and structural plots for a rational error correction model
#'
#' Nine panels, selected by `which`, in the manner of [plot.lm()]. The first
#' four are drawn by default: two residual diagnostics and the two views of
#' the estimated lead weights.
#'
#' The lead weights are the reason the last group of panels exists. They are
#' not estimated freely - they are a nonlinear function of the reduced-form
#' coefficients through \eqn{A(L)} - so their shape is what the cross-equation
#' restrictions actually deliver, and `sum(d_i)` together with the mean lead
#' does not convey it.
#'
#' @param x An object of class `"recm"` returned by [recm()].
#' @param which Integer vector selecting panels from `1:9`; see below.
#'   Defaults to `1:4`.
#' @param caption Titles for the nine panels, recycled to length nine.
#' @param main Title actually used, recycled over the selected panels.
#'   `NULL`, the default, takes the titles from `caption`.
#' @param horizon Horizon for the lead weight panels (3 and 4) and the
#'   adjustment path (9). `NULL`, the default, stops each of them where the
#'   sequence has decayed to nothing, so a fast adjustment is not drawn as a
#'   spike followed by a long flat stretch.
#' @param lag_max Highest lag shown in the residual autocorrelation panel.
#'   `NULL` uses `stats::acf()`'s own rule.
#' @param id_n Number of extreme residuals labelled with their time index in
#'   panels 1, 5 and 6. Zero labels none.
#' @param ask Pause between panels. Defaults to `TRUE` when more panels were
#'   asked for than the current layout can hold and the device is interactive,
#'   which is [plot.lm()]'s rule.
#' @param ... Further graphical parameters passed to the underlying plot.
#'
#' @section Panels:
#' \describe{
#'   \item{1}{Residuals against fitted \eqn{\Delta y}, with a lowess smooth.}
#'   \item{2}{Histogram of the residuals, against the normal density implied
#'     by the reported residual standard error and a kernel estimate. The
#'     reference normal is centred at zero, not at the residual mean: no
#'     intercept is fitted, so a visible shift is worth seeing.}
#'   \item{3}{Lead weights \eqn{d_i} against the horizon, with the mean lead
#'     marked.}
#'   \item{4}{Cumulative lead weights as a share of \eqn{\sum_i d_i}, with the
#'     horizons carrying half and nine tenths of the forward response marked.
#'     If the total loading is zero the raw cumulative sum is drawn instead.}
#'   \item{5}{Residuals through the sample, with a two standard error band.}
#'   \item{6}{Normal quantile plot of the standardised residuals.}
#'   \item{7}{Residual autocorrelation with its \eqn{\pm 1.96/\sqrt{T}} bands
#'     and the HAC bandwidth used by the reported standard errors.}
#'   \item{8}{\eqn{\Delta y} and its fitted value through the sample.}
#'   \item{9}{Deterministic response of the error correction gap to a unit
#'     displacement with the target held flat, with the half-life marked. At
#'     `m > 1` this path need not be monotone, which the half-life alone
#'     cannot say.}
#' }
#'
#' @return `x`, invisibly.
#'
#' @seealso [recm()], [summary.recm()][recm-methods].
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' dystar <- as.numeric(stats::filter(rnorm(n, 0.5, 0.4), 0.6, "recursive"))
#' ystar <- 100 + cumsum(dystar)
#' y <- numeric(n)
#' y[1] <- ystar[1]
#' for (t in 2:n) {
#'   y[t] <- y[t - 1] + 0.3 * (ystar[t - 1] - y[t - 1]) +
#'     0.7 * dystar[t] + rnorm(1, 0, 0.1)
#' }
#' fit <- recm(y, ystar, data.frame(y = y, ystar = ystar))
#'
#' op <- par(mfrow = c(2, 2))
#' plot(fit)
#' par(op)
#'
#' # The structural panels on their own.
#' plot(fit, which = c(3, 4, 9), ask = FALSE)
#'
#' @export
plot.recm <- function(x,
                      which = 1:4,
                      caption = c("Residuals vs fitted",
                                  "Histogram of residuals",
                                  "Lead weights by horizon",
                                  "Cumulative lead weights",
                                  "Residuals over time",
                                  "Normal Q-Q",
                                  "Residual autocorrelation",
                                  "Actual and fitted",
                                  "Adjustment to a unit gap"),
                      main = NULL,
                      horizon = NULL,
                      lag_max = NULL,
                      id_n = 3L,
                      ask = prod(graphics::par("mfcol")) < length(which) &&
                        grDevices::dev.interactive(),
                      ...) {
  n_panel <- 9L
  if (!is.numeric(which) || !length(which) || anyNA(which) ||
        any(which != round(which)) || any(which < 1L | which > n_panel)) {
    stop(
      "`which` must be a non-empty vector of whole numbers in 1:", n_panel,
      ". See ?plot.recm for what each panel shows.",
      call. = FALSE
    )
  }
  panels <- as.integer(which)

  check_count <- function(v, arg) {
    if (is.null(v)) {
      return(invisible(NULL))
    }
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != round(v) ||
          v < 1) {
      stop("`", arg, "` must be a single positive whole number or NULL.",
           call. = FALSE)
    }
    invisible(NULL)
  }
  check_count(horizon, "horizon")
  check_count(lag_max, "lag_max")
  if (!is.numeric(id_n) || length(id_n) != 1L || is.na(id_n) ||
        id_n != round(id_n) || id_n < 0) {
    stop("`id_n` must be a single non-negative whole number.", call. = FALSE)
  }
  id_n <- as.integer(id_n)
  horizon <- if (is.null(horizon)) NULL else as.integer(horizon)
  lag_max <- if (is.null(lag_max)) NULL else as.integer(lag_max)

  caption <- rep_len(as.character(caption), n_panel)
  titles <- if (is.null(main)) {
    caption[panels]
  } else {
    rep_len(as.character(main), length(panels))
  }

  if (isTRUE(ask)) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask), add = TRUE)
  }

  pal <- recm_pal()
  # Computed once: both weight panels share it, and it is the only part of
  # the drawing that costs anything.
  w <- if (any(panels %in% c(3L, 4L))) recm_weights(x, horizon) else NULL

  for (k in seq_along(panels)) {
    tt <- titles[k]
    switch(
      panels[k],
      recm_panel_resid_fit(x, tt, id_n, pal, ...),
      recm_panel_hist(x, tt, pal, ...),
      recm_panel_leads(x, tt, w, pal, ...),
      recm_panel_cum_leads(x, tt, w, pal, ...),
      recm_panel_resid_time(x, tt, id_n, pal, ...),
      recm_panel_qq(x, tt, id_n, pal, ...),
      recm_panel_acf(x, tt, lag_max, pal, ...),
      recm_panel_actual_fitted(x, tt, pal, ...),
      recm_panel_adjustment(x, tt, horizon, pal, ...)
    )
  }

  invisible(x)
}
