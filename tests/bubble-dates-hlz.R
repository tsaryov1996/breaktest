# Independent, direct-regression checks for the HLZ dating criterion.
library(breaktest)

checks <- 0L
check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  checks <<- checks + 1L
  invisible(TRUE)
}
near <- function(x, y, tolerance = 1e-9) {
  isTRUE(all.equal(unname(x), unname(y), tolerance = tolerance,
                   check.attributes = FALSE))
}
error <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  check(!is.null(msg) && grepl(pattern, msg), paste("Expected error:", pattern))
}

# No prefix sums are used here. Each prefix and candidate window is independently
# fitted with lm.fit on the actual lag/current rows, including the initial zero.
oracle <- function(y, min.window = 0.1, epsilon = 0.01, trim = 0,
                   ordinary.variance = FALSE) {
  N <- length(y)
  delta <- numeric(N)
  C <- numeric(N)
  for (t in 2:N) {
    delta[t] <- if (y[t] > y[t - 1L]) 1 else -1
    C[t] <- C[t - 1L] + delta[t]
  }
  lagged <- c(0, C[-N])
  regress <- function(a, b) {
    rows <- seq.int(a + 1L, b)
    x <- lagged[rows]
    d <- delta[rows]
    q <- sum(x^2)
    xy <- sum(x * d)
    yy <- sum(d^2)
    if (q == 0) return(list(q = 0, beta = 0, SSR = yy, xy = xy, yy = yy))
    fit <- stats::lm.fit(matrix(x, ncol = 1L), d)
    list(q = q, beta = if (xy == 0) 0 else unname(fit$coefficients),
         SSR = sum(fit$residuals^2), xy = xy, yy = yy)
  }
  prefix <- numeric(N + 1L)
  prefix.num <- numeric(N + 1L)
  prefix.den <- rep(1, N + 1L)
  for (j in 2:N) {
    fit <- regress(0, j)
    prefix[j + 1L] <- j * fit$SSR / (j - 1L)
    # The sign data give integer cross-products. Exact rational numerators
    # identify mathematical zero differences that QR residual roundoff can
    # turn into a spurious, tiny positive variance. No positive value is floored.
    prefix.num[j + 1L] <- if (fit$q) j * (fit$yy * fit$q - fit$xy^2) else j * fit$yy
    prefix.den[j + 1L] <- (j - 1L) * if (fit$q) fit$q else 1
  }
  first <- ceiling(N * trim)
  last <- floor(N * (1 - trim))
  width <- max(3L, ceiling(N * min.window))
  rows <- list()
  for (b in seq.int(first + width, last)) {
    for (a in seq.int(first, b - width)) {
      fit <- regress(a, b)
      variance <- if (ordinary.variance) fit$SSR / (b - a - 1L) else
        (prefix[b + 1L] - prefix[a + 1L]) / (b - a - 1L)
      if (!ordinary.variance) {
        terms <- c(prefix.num[b + 1L] * prefix.den[a + 1L],
                   prefix.num[a + 1L] * prefix.den[b + 1L])
        stopifnot(all(abs(terms) < 2^53))
        if (terms[1L] == terms[2L]) variance <- 0
      }
      valid <- fit$q > 0 && is.finite(variance) && variance > 0
      score <- if (valid) fit$beta * sqrt(fit$q) / variance^(epsilon / 2) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        start = a, end = b, q = fit$q, variance = variance,
        statistic = score, valid = valid
      )
    }
  }
  candidates <- do.call(rbind, rows)
  valid <- candidates[candidates$valid, ]
  list(C = C, candidates = candidates, valid = valid,
       best = if (nrow(valid)) valid[which.max(valid$statistic), ] else NULL)
}

compare <- function(y, min.window = 0.1, epsilon = 0.01, trim = 0) {
  expected <- oracle(y, min.window, epsilon, trim)
  fit <- bubble.dates.HLZ(y, min.window, epsilon, trim)
  check(inherits(fit, "bt_bubbleDates"), "Result class")
  check(near(fit$C.t, expected$C), "Global cumulative signs")
  check(near(fit$statistic, expected$best$statistic),
        paste("Maximum criterion:", fit$statistic, "versus", expected$best$statistic,
              "for y =", paste(y, collapse = ","), "epsilon =", epsilon))
  check(identical(names(fit$break.points), c("emerge", "collapse")), "Only two HLZ dates")
  check(near(fit$window, fit$break.points), "Window coordinates are row boundaries")
  check(fit$n.obs == length(y), "HLZ uses N observed rows")
  check(fit$min.length == max(3L, ceiling(length(y) * min.window)), "Minimum length rounding")
  check(fit$diagnostics$candidate.windows == nrow(expected$candidates), "Candidate count")
  check(fit$diagnostics$singular.windows == sum(expected$candidates$q == 0), "Singular count")
  check(fit$diagnostics$valid.windows == nrow(expected$valid), "Valid-window count")
  check(fit$diagnostics$nonpositive.variance.windows ==
          sum(expected$candidates$q > 0 & expected$candidates$variance <= 0),
        "Nonpositive modified-variance count")
  for (i in seq_len(nrow(fit$profile))) {
    row <- fit$profile[i, ]
    candidates <- expected$valid[expected$valid$end == row$end, ]
    check(row$valid.windows == nrow(candidates), "Endpoint valid-window count")
    if (!nrow(candidates)) {
      check(is.na(row$start) && is.na(row$statistic) && is.na(row$variance),
            "Invalid profile row is missing")
    } else {
      check(near(row$statistic, max(candidates$statistic)), "Endpoint maximum")
      chosen <- candidates[candidates$start == row$start, ]
      check(nrow(chosen) == 1L && near(row$variance, chosen$variance),
            "Modified variance of chosen window")
    }
  }
  invisible(fit)
}

set.seed(2019)
for (epsilon in c(0.01, 0.5, 1)) {
  compare(cumsum(stats::rnorm(36)), min.window = 0.15, epsilon = epsilon,
          trim = 0.05)
  compare(cumsum(stats::rnorm(28)), min.window = 0.2, epsilon = epsilon, trim = 0)
}

# The public default includes the final sample row. A 5% outer trim would
# exclude it and select a different window on this path.
y <- 0:6
default <- bubble.dates.HLZ(y)
expected <- oracle(y)
check(identical(default, bubble.dates.HLZ(y, trim = 0)), "Default has no outer trim")
check(near(default$break.points, c(2, length(y))), "Default can reach the final row")
check(near(default$statistic, expected$best$statistic), "Default agrees with direct oracle")
check(bubble.dates.HLZ(y, trim = 0.05)$break.points["collapse"] < length(y),
      "Optional outer trim excludes the final row")

# Small exhaustive sign paths supply repeated regressors and exact zero cross
# products. A lm.fit oracle may perturb an exact tie, so compare the value and
# the variance at the implementation's chosen maximizing candidate.
for (mask in 0:31) {
  signs <- ifelse(as.integer(intToBits(mask))[1:5] == 0L, -1, 1)
  compare(c(0, cumsum(signs)), min.window = 0.1, epsilon = 0.5, trim = 0)
}

# A fixture that separates the modified prefix variance from ordinary window
# residual variance. The two criteria select different, unique maxima.
signs <- c(-1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1,
           -1, -1, 1, 1, 1, -1, 1, 1, -1, -1, 1, -1)
y <- c(0, cumsum(signs))
fit <- compare(y, min.window = 0.16, epsilon = 1, trim = 0)
ordinary <- oracle(y, min.window = 0.16, epsilon = 1, trim = 0,
                   ordinary.variance = TRUE)
check(near(fit$break.points, c(8, 14)), "HLZ fixture dates")
check(near(unlist(ordinary$best[c("start", "end")]), c(11, 15)), "Ordinary-variance fixture dates")
check(!near(fit$statistic, ordinary$best$statistic), "Criterion is not ordinary window t")

# Zero increments map to -1, rather than sign(0)=0. The initial row is 0.
z <- c(0, 0, -1, 1, 1, 2, 0, 0, 3)
zfit <- compare(z, trim = 0)
check(near(zfit$C.t, c(0, -1, -2, -1, -2, -1, -2, -3, -2)), "HLZ sign(0) convention")
set.seed(92)
y <- cumsum(stats::rnorm(40))
reference <- bubble.dates.HLZ(y)
for (converted in list(7 + 2 * y, -y, as.matrix(y), stats::ts(y))) {
  transformed <- bubble.dates.HLZ(converted)
  check(near(transformed$break.points, reference$break.points) &&
          near(transformed$statistic, reference$statistic), "Scale/location/shape invariance")
}
huge <- rep(c(-1e308, 1e308), 5)
check(near(bubble.dates.HLZ(huge, trim = 0)$C.t,
           bubble.dates.HLZ(sign(huge), trim = 0)$C.t), "No overflow in sign differences")

# Exactly tied maxima with cross-product zero: earliest end, then start.
tied <- bubble.dates.HLZ(c(0, -1, 0, -1, -2), trim = 0, epsilon = 1)
check(near(tied$break.points, c(0, 5)) && tied$statistic == 0, "Exact-tie policy")
check(isTRUE(unname(tied$diagnostics$boundary.before.sample)), "Boundary zero has no row/date")
tied.ends <- bubble.dates.HLZ(c(0, -1, 0, -1, -2, -1, -2), trim = 0, epsilon = 1)
check(near(tied.ends$break.points, c(0, 5)) && tied.ends$statistic == 0,
      "Exact ties across endpoints choose the earliest end")
set.seed(729)
rng.before <- .Random.seed
invisible(bubble.dates.HLZ(c(0, -1, 0, -1, -2, -1, -2), trim = 0))
check(identical(.Random.seed, rng.before), "Deterministic dating does not consume random numbers")

# Monotone data can yield negative modified variances despite positive ordinary
# residual sums of squares. They must be excluded, not floored to positive.
monotone <- compare(0:4, trim = 0, epsilon = 1)
check(monotone$diagnostics$nonpositive.variance.windows == 1L, "Skip negative modified variance")
error(bubble.dates.HLZ(0:6, trim = 0.2), "positive modified residual variance")

for (bad in list(NULL, numeric(), 1:4, rep(1, 20), c(1:7, NA), c(1:7, Inf),
                 c(1:7, NaN), letters[1:8], TRUE, matrix(1:12, ncol = 2),
                 array(1:12, dim = c(12, 1, 1)), (1:8) + 1i)) {
  error(bubble.dates.HLZ(bad), "finite numeric|constant")
}
y <- seq_len(20)
for (bad in list(0, -0.1, 1.1, NA, NaN, Inf, numeric(), c(0.1, 0.2), TRUE, "0.1", 0.1 + 1i)) {
  error(bubble.dates.HLZ(y, min.window = bad), "min.window")
  error(bubble.dates.HLZ(y, epsilon = bad), "epsilon")
}
for (bad in list(-0.1, 0.5, NA, NaN, Inf, numeric(), c(0.1, 0.2), TRUE, "0.1", 0.1 + 1i)) {
  error(bubble.dates.HLZ(y, trim = bad), "trim")
}
error(bubble.dates.HLZ(y, min.window = 0.95, trim = 0.1), "No window")
compare(y, min.window = 1, trim = 0)
cat("HLZ date-estimator checks passed:", checks, "\n")
