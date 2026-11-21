#' Bubble dates with optional volatility correction
#'
#' Estimate the emergence, collapse and recovery of one bubble using
#' sequential OLS or adaptive WLS sample splitting.
#'
#' @param y A finite numeric vector, univariate time series or one-column
#'   matrix containing the initial level and at least eight observations.
#' @param trim Fraction of the full sample excluded at each end of a search
#'   window, strictly between zero and one quarter.
#' @param method Use adaptive weighted least squares (`"wls"`) or ordinary
#'   least squares (`"ols"`).
#' @param bandwidth Positive Gaussian bandwidth as a fraction of the full
#'   sample size. NULL selects it by leave-one-out cross-validation.
#' @param bandwidth.grid Positive bandwidth candidates, used only when
#'   `bandwidth` is NULL. The default is 25 logarithmically spaced values
#'   from `1 / (length(y) - 1)` to `0.5`. Bandwidth arguments require WLS.
#'
#' @details
#' With input \eqn{y_0,\ldots,y_T}, `y[1]` is the initial level. A date
#' \eqn{k} is the last equation before a regime change. `Te_est`, `Tc_est`
#' and `Tr_est` use this convention; `break.points` adds one to give row
#' indices in `y`, suitable for [plot_bubble].
#'
#' Each split minimises the weighted residual sum of squares of two AR(1)
#' regressions without intercepts. Collapse is estimated on the full sample;
#' emergence on equations `1:Tc_est`; recovery on `(Tc_est + 1):T`, keeping
#' the lag at `Tc_est`. For a window \eqn{(a,b]}, candidates run from
#' \eqn{a+\max(h+1,2)} to \eqn{b-\max(h,2)}, where
#' \eqn{h=\lfloor T\,trim\rfloor}. This follows the paper's recovery
#' search convention and requires at least two equations in each regime.
#' Exact ties select the earliest date.
#'
#' WLS first obtains OLS dates. In the pilot residual regression,
#' `diff(y)` has an intercept and a lag coefficient in each of the explosive
#' and collapsing regimes, and fitted value zero outside them (equation 4
#' of Kurozumi and Skrobotov). Gaussian leave-one-out smoothing of squared
#' residuals supplies inverse-variance weights for all three searches.
#' CV minimises the squared prediction error of the squared residuals;
#' the explicit criterion and default grid are implementation choices since
#' the paper does not specify them. Exact CV ties select the smallest
#' bandwidth. Extremely small bandwidths use the nearest-neighbour kernel
#' limit if every Gaussian weight underflows. No variance floor is added.
#'
#' If a split leaves too little data or cannot identify emergence or recovery,
#' that date is NA and the other dates are retained (`status = "partial"`).
#' If the OLS pilot is incomplete, WLS cannot estimate its variance:
#' all final dates are NA (`status = "incomplete.pilot"`) and the available
#' OLS dates remain in `pilot`. An unidentified collapse, a rank-deficient
#' pilot regression or undefined variances cause an error.
#'
#' These are conditional point estimates for one four-regime bubble, not a
#' bubble-existence test or confidence sets. The model assumes independent
#' innovations with a positive deterministic variance function. WLS need not
#' outperform OLS for every volatility path, and its dates do not validate
#' weighted confidence sets from [cs.bubble.emerge] and related functions.
#'
#' @return An object of class `bt_bubbleDates` with:
#' * `method`, `status`: requested method and `"complete"`, `"partial"` or
#'   `"incomplete.pilot"`; `reason` explains missing dates.
#' * `Te_est`, `Tc_est`, `Tr_est`: equation indices; `break.points`: corresponding
#'   row indices named `emerge`, `collapse`, `recovery`. Missing dates are NA.
#' * `coefficients`: four AR coefficients at the selected dates, or NA where
#'   a regime boundary is missing; `phi_a`, `phi_b`: the two middle coefficients.
#'   `regime.check` checks `phi_a > 1` and `0 < phi_b < 1`, not significance.
#' * `pilot`: preliminary OLS dates, row indices, status and reasons; for WLS,
#'   also residuals and coefficients if its residual regression was fitted.
#' * `variance`, `bandwidth`, `bandwidth.grid`, `cv`: WLS smoothing results,
#'   NULL when unused. CV scores use residuals divided by `scale`.
#' * `search`: final candidate-date/objective profiles, NULL if WLS was not run.
#'   Objectives use `y / scale` and inverse variances divided by their maximum.
#' * `y`, `trim`, `n.obs`, `scale`: data, trimming fraction, number of equations
#'   and numerical scaling factor.
#'
#' @references
#' Kurozumi, E. and A. Skrobotov (2026). Improving the accuracy of bubble date
#' estimators under time-varying volatility. arXiv:2306.02977v2, Sections 3--6.
#' https://arxiv.org/abs/2306.02977v2.
#'
#' @examples
#' set.seed(42)
#' y <- numeric(201)
#' y[1] <- 100
#' phi <- rep(c(1, 1.04, 0.94, 1), c(70, 50, 40, 40))
#' for (t in seq_along(phi)) y[t + 1] <- phi[t] * y[t] + rnorm(1)
#' bubble.dates.NBCN(y, bandwidth = 0.1)$break.points
#'
#' @export
bubble.dates.NBCN <- function(y, trim = 0.05, method = c("wls", "ols"),
                             bandwidth = NULL, bandwidth.grid = NULL) {
  method <- match.arg(method)
  if (!is.numeric(y) || is.complex(y) ||
      (!is.null(dim(y)) && (length(dim(y)) != 2L || ncol(y) != 1L)) ||
      length(y) < 9L || any(!is.finite(y))) {
    stop("y must contain at least nine finite numeric values in one series.")
  }
  y <- as.numeric(y)
  if (!is.numeric(trim) || is.complex(trim) || length(trim) != 1L ||
      !is.finite(trim) || trim <= 0 || trim >= 0.25) {
    stop("trim must be a single number in (0, 0.25).")
  }
  if (all(y == y[1L])) stop("A constant series does not identify bubble dates.")
  if (method == "ols" && (!is.null(bandwidth) || !is.null(bandwidth.grid))) {
    stop("Bandwidth arguments require method = 'wls'.")
  }
  if (!is.null(bandwidth) && !is.null(bandwidth.grid)) {
    stop("Supply bandwidth or bandwidth.grid, not both.")
  }
  check.bandwidth <- function(x, scalar = FALSE) {
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || !length(x) ||
        (scalar && length(x) != 1L) || any(!is.finite(x)) || any(x <= 0)) {
      stop("Bandwidths must be finite positive numbers.")
    }
  }
  if (!is.null(bandwidth)) check.bandwidth(bandwidth, scalar = TRUE)
  if (!is.null(bandwidth.grid)) check.bandwidth(bandwidth.grid)

  N <- length(y) - 1L
  scale <- max(abs(y))
  z <- y / scale
  x <- z[-length(z)]
  dy <- diff(z)
  precision <- rep(1, N)
  fit <- .bubble.nbcn.split(x, dy, precision, trim)
  pilot <- c(as.list(fit$dates), list(
    break.points = .bubble.nbcn.rows(fit$dates),
    status = fit$status, reason = fit$reason
  ))
  variance <- cv <- NULL

  if (method == "wls" && fit$status == "partial") {
    fit$dates[] <- NA_integer_
    fit$status <- "incomplete.pilot"
    fit$reason <- c(pilot = "The OLS pilot is incomplete; WLS was not estimated.")
    fit$search <- NULL
    bandwidth <- bandwidth.grid <- NULL
  } else if (method == "wls") {
    residual.fit <- .bubble.nbcn.residuals(x, dy, fit$dates)
    residuals <- residual.fit$residuals
    pilot$residuals <- residuals * scale
    pilot$coefficients <- residual.fit$coefficients
    pilot$coefficients[1:2] <- pilot$coefficients[1:2] * scale

    if (is.null(bandwidth)) {
      if (is.null(bandwidth.grid)) {
        bandwidth.grid <- exp(seq(log(1 / N), log(0.5), length.out = 25L))
      }
      bandwidth.grid <- sort(unique(bandwidth.grid))
      cv <- vapply(bandwidth.grid, function(b) {
        v <- .bubble.nbcn.variance(residuals, b)
        if (any(!is.finite(v)) || any(v <= 0)) return(Inf)
        mean((residuals^2 - v)^2)
      }, numeric(1))
      if (!any(is.finite(cv))) {
        stop("No bandwidth gives strictly positive finite variances.")
      }
      bandwidth <- bandwidth.grid[which.min(cv)]
    }
    variance.scaled <- .bubble.nbcn.variance(residuals, bandwidth)
    if (any(!is.finite(variance.scaled)) || any(variance.scaled <= 0)) {
      stop("Estimated variances must be strictly positive and finite.")
    }
    precision <- min(variance.scaled) / variance.scaled
    if (any(!is.finite(precision)) || any(precision <= 0)) {
      stop("Estimated variance ratios exceed numerical precision.")
    }
    variance <- (variance.scaled * scale) * scale
    if (any(!is.finite(variance)) || any(variance <= 0)) {
      stop("Variance is outside numeric range; rescale y before estimation.")
    }
    fit <- .bubble.nbcn.split(x, dy, precision, trim)
  }

  bounds <- c(0L, fit$dates, N)
  coefficients <- vapply(seq_len(4L), function(j) {
    if (anyNA(bounds[c(j, j + 1L)])) return(NA_real_)
    ii <- seq.int(bounds[j] + 1L, bounds[j + 1L])
    .bubble.nbcn.ar(x[ii], dy[ii], precision[ii])$coefficient
  }, numeric(1))
  names(coefficients) <- c("normal.before", "explosive", "collapsing", "normal.after")
  result <- list(
    y = y, method = method, estimator = "NBCN", trim = trim, n.obs = N,
    status = fit$status, reason = fit$reason,
    Te_est = unname(fit$dates["Te_est"]),
    Tc_est = unname(fit$dates["Tc_est"]),
    Tr_est = unname(fit$dates["Tr_est"]),
    break.points = .bubble.nbcn.rows(fit$dates),
    coefficients = coefficients,
    phi_a = unname(coefficients[2L]), phi_b = unname(coefficients[3L]),
    regime.check = c(explosive = unname(coefficients[2L]) > 1,
                     collapsing = unname(coefficients[3L]) > 0 &
                       unname(coefficients[3L]) < 1),
    pilot = pilot, variance = variance, bandwidth = bandwidth,
    bandwidth.grid = bandwidth.grid, cv = cv, scale = scale, search = fit$search
  )
  class(result) <- "bt_bubbleDates"
  result
}


.bubble.nbcn.rows <- function(dates) {
  stats::setNames(as.integer(dates + 1L), c("emerge", "collapse", "recovery"))
}


.bubble.nbcn.ar <- function(x, dy, precision) {
  denominator <- sum(precision * x^2)
  if (!is.finite(denominator) || denominator <= 0) {
    return(list(coefficient = NA_real_, SSR = Inf))
  }
  delta <- sum(precision * x * dy) / denominator
  list(coefficient = 1 + delta, SSR = sum(precision * (dy - delta * x)^2))
}


.bubble.nbcn.split <- function(x, dy, precision, trim) {
  N <- length(x)
  h <- floor(N * trim)
  search <- function(a, b, name) {
    lower <- a + max(h + 1L, 2L)
    upper <- b - max(h, 2L)
    if (lower > upper) {
      return(list(
        date = NA_integer_,
        profile = data.frame(date = integer(), objective = numeric()),
        reason = paste("Too few observations for", name, "after trimming.")
      ))
    }
    candidates <- seq.int(lower, upper)
    objective <- vapply(candidates, function(k) {
      left <- seq.int(a + 1L, k)
      right <- seq.int(k + 1L, b)
      .bubble.nbcn.ar(x[left], dy[left], precision[left])$SSR +
        .bubble.nbcn.ar(x[right], dy[right], precision[right])$SSR
    }, numeric(1))
    profile <- data.frame(date = candidates, objective = objective)
    if (!any(is.finite(objective))) {
      return(list(date = NA_integer_, profile = profile,
                  reason = paste("No identified AR fit for", name, "in the search window.")))
    }
    list(date = candidates[which.min(objective)], profile = profile, reason = NULL)
  }
  collapse <- search(0L, N, "collapse")
  if (is.na(collapse$date)) stop("No identified collapse date in the search window.")
  emerge <- search(0L, collapse$date, "emergence")
  recovery <- search(collapse$date, N, "recovery")
  dates <- c(Te_est = emerge$date, Tc_est = collapse$date, Tr_est = recovery$date)
  reason <- c(emerge = emerge$reason, recovery = recovery$reason)
  list(dates = dates, status = if (anyNA(dates)) "partial" else "complete",
       reason = reason,
       search = list(emerge = emerge$profile, collapse = collapse$profile,
                     recovery = recovery$profile))
}


.bubble.nbcn.residuals <- function(x, dy, dates) {
  t <- seq_along(x)
  explosive <- as.numeric(t > dates["Te_est"] & t <= dates["Tc_est"])
  collapsing <- as.numeric(t > dates["Tc_est"] & t <= dates["Tr_est"])
  design <- cbind(mu.a = explosive, mu.b = collapsing,
                  delta.a = explosive * x, delta.b = collapsing * x)
  fit <- stats::lm.fit(design, dy)
  if (fit$rank < ncol(design)) stop("The pilot residual regression is rank deficient.")
  list(residuals = fit$residuals, coefficients = fit$coefficients)
}


.bubble.nbcn.variance <- function(residuals, bandwidth) {
  N <- length(residuals)
  squared <- residuals^2
  time <- seq_len(N)
  log.weights <- -0.5 * ((seq_len(N - 1L) / N) / bandwidth)^2
  kernel <- numeric(N)
  # Kernel[1] is the excluded own observation. Normalise at the nearest lag.
  if (!is.finite(log.weights[1L])) {
    kernel[2L] <- 1
  } else {
    kernel[-1L] <- exp(log.weights - log.weights[1L])
  }
  vapply(time, function(t) {
    weights <- kernel[abs(time - t) + 1L]
    sum(weights * squared) / sum(weights)
  }, numeric(1))
}
