#' Sign-based bubble date estimation
#'
#' Estimate the emergence and collapse of one explosive episode using the
#' modified sign-based criterion in Section 6 of Harvey, Leybourne and Zu (2020).
#'
#' @param y A finite numeric series or one-column numeric matrix. Observations
#'   must be consecutive; missing values are not removed.
#' @param min.window Minimum window length as a fraction of `length(y)`.
#' @param epsilon Power in the modified residual variance, in \eqn{(0,1]}.
#' @param trim Fraction excluded at each outer end of the search. The default
#'   allows both sample boundaries; `trim = 0.05` excludes 5 percent at each end.
#'
#' @details
#' For \eqn{N=length(y)}, form the cumulative signs
#' \eqn{C_1=0}, \eqn{C_t=\sum_{i=2}^t[1-2I(y_i-y_{i-1}\leq0)]}.
#' Zero increments have sign -1. On each window \eqn{(a,b]}, regress
#' \eqn{\Delta C_t} on \eqn{C_{t-1}} without an intercept, using the same
#' cumulative series throughout. The dates maximize
#' \deqn{\hat\rho(a,b)\sqrt{\sum_{t=a+1}^b C_{t-1}^2}/
#'        [\tilde s^2(a,b)]^{\epsilon/2},}
#' where \eqn{\tilde s^2(a,b)=[b\hat s^2(0,b)-a\hat s^2(0,a)]/(b-a-1)}
#' and \eqn{\hat s^2(0,j)=SSR(0,j)/(j-1)}.
#'
#' Candidate boundaries satisfy \eqn{a\geq ceiling(N\,trim)},
#' \eqn{b\leq floor(N(1-trim))}, and
#' \eqn{b-a\geq max(3,ceiling(N\,min.window))}.
#' The initial row has \eqn{C_0=\Delta C_1=0}; prefix terms for
#' \eqn{a=0,1} are zero. A prefix with a zero regressor uses coefficient zero.
#' Windows with zero regressor sum of squares or nonpositive modified variance
#' are excluded and counted in `diagnostics`. If none remain, the function
#' errors. Exact ties select the earliest end, then the earliest start.
#'
#' `break.points` gives boundaries in the rows of `y`: emergence is the last
#' row before the explosive window, and collapse is its last row. An emergence
#' of zero means before the first observed row; it has no observation or calendar
#' date. No index offset is needed. The `window` is open at its start and closed
#' at its end.
#'
#' These are conditional point estimates for one bubble, with no recovery date
#' or confidence sets. The function does not test whether a bubble exists.
#' The theory assumes independent innovations with zero mean and median,
#' positive deterministic time-varying volatility, and the other regularity
#' conditions of HLZ. The minimum window must not exceed the true episode length.
#'
#' @return An object of class `bt_bubbleDates`: a list with `y`, `method`, the
#'   tuning parameters, `n.obs`, `min.length`, `break.points` (emergence and
#'   collapse), `window`, the maximum `statistic`, `C.t`, and a `profile` of the
#'   best window ending at each candidate endpoint. `diagnostics` gives the
#'   numbers of candidate, singular, nonpositive-variance and valid windows,
#'   the number of zero increments, and whether emergence precedes the sample.
#'
#' @references
#' Harvey, D. I., S. J. Leybourne, and Y. Zu (2020).
#' Sign-Based Unit Root Tests for Explosive Financial Bubbles in the Presence
#' of Deterministically Time-Varying Volatility. Econometric Theory 36, 122--169.
#' \url{https://doi.org/10.1017/S0266466619000057}.
#'
#' @examples
#' set.seed(42)
#' y <- cumsum(rnorm(100))
#' fit <- bubble.dates.HLZ(y)
#' fit$break.points
#'
#' @export
bubble.dates.HLZ <- function(
  y,
  min.window = 0.1,
  epsilon = 0.01,
  trim = 0
) {
  if (!is.numeric(y) || is.complex(y) ||
      (!is.null(dim(y)) && (length(dim(y)) != 2L || ncol(y) != 1L)) ||
      length(y) < 5L || any(!is.finite(y))) {
    stop("y must contain at least five finite numeric observations in one series.")
  }
  y <- as.numeric(y)
  if (all(y == y[1L])) {
    stop("A constant series contains no dating information.")
  }
  is.scalar <- function(z) {
    is.numeric(z) && !is.complex(z) && length(z) == 1L && is.finite(z)
  }
  if (!is.scalar(min.window) || min.window <= 0 || min.window > 1) {
    stop("min.window must be a single number in (0, 1].")
  }
  if (!is.scalar(epsilon) || epsilon <= 0 || epsilon > 1) {
    stop("epsilon must be a single number in (0, 1].")
  }
  if (!is.scalar(trim) || trim < 0 || trim >= 0.5) {
    stop("trim must be a single number in [0, 0.5).")
  }

  N <- length(y)
  width <- max(3L, ceiling(N * min.window))
  first <- ceiling(N * trim)
  last <- floor(N * (1 - trim))
  if (last - first < width) {
    stop("No window satisfies trim and min.window.")
  }

  # Compare levels directly to avoid overflow in their differences.
  dC <- c(0, ifelse(y[-1L] > y[-N], 1, -1))
  C.t <- cumsum(dC)
  lag.C <- c(0, C.t[-N])
  prefix.xx <- c(0, cumsum(lag.C^2))
  prefix.xy <- c(0, cumsum(lag.C * dC))
  prefix.yy <- c(0, cumsum(dC^2))

  # Entry j + 1 contains j rows, including the initial zero row.
  prefix.SSR <- prefix.yy
  identified <- prefix.xx > 0
  prefix.SSR[identified] <- prefix.yy[identified] -
    prefix.xy[identified]^2 / prefix.xx[identified]
  prefix.SSR <- pmax(prefix.SSR, 0)
  scaled.variance <- numeric(N + 1L)
  j <- 2:N
  scaled.variance[j + 1L] <- j * prefix.SSR[j + 1L] / (j - 1L)

  ends <- seq.int(first + width, last)
  profile <- data.frame(
    start = rep(NA_integer_, length(ends)),
    end = as.integer(ends),
    statistic = rep(NA_real_, length(ends)),
    variance = rep(NA_real_, length(ends)),
    valid.windows = integer(length(ends))
  )
  candidates <- singular <- nonpositive <- 0
  for (i in seq_along(ends)) {
    b <- ends[i]
    a <- seq.int(first, b - width)
    Sxx <- prefix.xx[b + 1L] - prefix.xx[a + 1L]
    Sxy <- prefix.xy[b + 1L] - prefix.xy[a + 1L]
    variance <- (scaled.variance[b + 1L] - scaled.variance[a + 1L]) /
      (b - a - 1L)
    candidates <- candidates + length(a)
    singular <- singular + sum(Sxx <= 0)
    nonpositive <- nonpositive + sum(Sxx > 0 &
                                      (!is.finite(variance) | variance <= 0))
    valid <- Sxx > 0 & is.finite(variance) & variance > 0
    profile$valid.windows[i] <- sum(valid)
    if (!any(valid)) next

    statistic <- Sxy[valid] / sqrt(Sxx[valid]) / variance[valid]^(epsilon / 2)
    pick <- which.max(statistic)
    profile$start[i] <- a[valid][pick]
    profile$statistic[i] <- statistic[pick]
    profile$variance[i] <- variance[valid][pick]
  }
  if (!any(is.finite(profile$statistic))) {
    stop("No identified window has a positive modified residual variance.")
  }
  pick <- which.max(replace(profile$statistic, is.na(profile$statistic), -Inf))
  points <- c(emerge = profile$start[pick], collapse = profile$end[pick])
  result <- list(
    y = y,
    method = "HLZ",
    min.window = min.window,
    epsilon = epsilon,
    trim = trim,
    n.obs = N,
    min.length = width,
    break.points = points,
    window = c(start = unname(points[1L]), end = unname(points[2L])),
    statistic = profile$statistic[pick],
    C.t = C.t,
    profile = profile,
    diagnostics = list(
      candidate.windows = candidates,
      singular.windows = singular,
      nonpositive.variance.windows = nonpositive,
      valid.windows = sum(profile$valid.windows),
      zero.increments = sum(y[-1L] == y[-N]),
      boundary.before.sample = unname(points[1L] == 0L)
    )
  )
  class(result) <- "bt_bubbleDates"
  result
}
