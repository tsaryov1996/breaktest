# Independent numerical references and legacy regression tests; base R only.
library(breaktest)

checks <- character()
check <- function(condition, label) {
  if (!isTRUE(condition)) stop(label, call. = FALSE)
  checks <<- c(checks, label)
}
equal <- function(actual, expected, label, tolerance = 1e-8) {
  difference <- all.equal(actual, expected, tolerance = tolerance,
                          check.attributes = FALSE)
  check(isTRUE(difference), paste0(label, ": ",
                                   paste(difference, collapse = "; ")))
}
fails <- function(expression, label) {
  result <- tryCatch(force(expression), error = identity)
  check(inherits(result, "error"), label)
}

# These references use lm.fit / QR, independently of all package ADF,
# detrending, break-grid and bootstrap helpers.
adf_reference <- function(z, lag = 0L) {
  n <- length(z)
  times <- seq.int(lag + 2L, n)
  dz <- c(NA_real_, diff(z))
  x <- matrix(z[times - 1L], ncol = 1L)
  for (j in seq_len(lag)) x <- cbind(x, dz[times - j])
  fit <- lm.fit(x, dz[times])
  stopifnot(fit$rank == ncol(x))
  sigma2 <- sum(fit$residuals^2) / (nrow(x) - fit$rank)
  covariance <- chol2inv(qr.R(fit$qr))
  position <- match(1L, fit$qr$pivot)
  unname(fit$coefficients[1L] / sqrt(sigma2 * covariance[position, position]))
}

select_lag_reference <- function(z, max.lag) {
  times <- seq.int(max.lag + 2L, length(z))
  dz <- c(NA_real_, diff(z))
  x <- matrix(z[times - 1L], ncol = 1L)
  criterion <- numeric(max.lag + 1L)
  # All MAIC candidates use the same sample, with the variance divisor N.
  # The final DF t-ratio is subsequently refitted using the selected lag.
  for (lag in 0:max.lag) {
    if (lag > 0L) x <- cbind(x, dz[times - lag])
    fit <- lm.fit(x, dz[times])
    sigma2 <- sum(fit$residuals^2) / length(times)
    tau <- unname(fit$coefficients[1L]^2 * sum(z[times - 1L]^2) / sigma2)
    criterion[lag + 1L] <- log(sigma2) + 2 * (ncol(x) + tau) / length(times)
  }
  which.min(criterion) - 1L
}

minimum_reference <- function(y, breaks, trim, max.lag = 0L) {
  y <- as.numeric(y)
  n <- length(y)
  time <- seq_len(n)
  tuples <- t(combn(seq.int(floor(trim * n), floor((1 - trim) * n)), breaks))
  tuples <- tuples[apply(tuples, 1L, function(b) all(diff(b) >= floor(trim * n))),
                   , drop = FALSE]
  cbar <- c(`2` = 21.5, `3` = 25.5)[as.character(breaks)]
  rho <- 1 - cbar / n
  selected_lags <- integer(nrow(tuples))
  statistics <- vapply(seq_len(nrow(tuples)), function(i) {
    x <- cbind(1, time,
               vapply(tuples[i, ], function(b) pmax(time - b, 0), numeric(n)))
    quasi_x <- rbind(x[1L, , drop = FALSE],
                     x[-1L, , drop = FALSE] - rho * x[-n, , drop = FALSE])
    quasi_y <- c(y[1L], y[-1L] - rho * y[-n])
    fit <- lm.fit(quasi_x, quasi_y)
    stopifnot(fit$rank == ncol(x))
    z <- drop(y - x %*% fit$coefficients)
    lag <- select_lag_reference(lm.fit(x, y)$residuals, max.lag)
    selected_lags[i] <<- lag
    adf_reference(z, lag)
  }, numeric(1L))
  best <- which.min(statistics)
  list(statistic = statistics[best], dates = unname(tuples[best, ]),
       lag = selected_lags[best], grid = tuples, cbar = unname(cbar),
       statistics = statistics, lags = selected_lags)
}

set.seed(9271)
n <- 28L
time <- seq_len(n)
y <- cumsum(rnorm(n)) + 0.035 * time +
  0.35 * pmax(time - 8L, 0) - 0.55 * pmax(time - 19L, 0)
set.seed(2319)
n3 <- 30L
time3 <- seq_len(n3)
y3 <- cumsum(rnorm(n3)) + 0.3 * pmax(time3 - 7L, 0) -
  0.6 * pmax(time3 - 15L, 0) + 0.4 * pmax(time3 - 23L, 0)
trim <- 0.2
bt <- asNamespace("breaktest")
mdf <- get("MDF.mlt", bt)

# Compare every reported estimate with a separately enumerated full grid.
for (m in 2:3) {
  series <- if (m == 2L) y else y3
  for (max_lag in 0:2) {
    ref <- minimum_reference(series, m, trim, max_lag)
    obj <- uroot.CHLT(series, max.lag = max_lag, trim = trim, breaks = m)
    label <- paste("breaks", m, "max.lag", max_lag)
    equal(obj$ADF$statistic, ref$statistic, paste(label, "DF minimum"))
    equal(obj$params$break.points, ref$dates, paste(label, "break estimates"), 0)
    equal(obj$params$lag, ref$lag, paste(label, "OLS-MAIC lag"), 0)
    equal(obj$params$configurations, nrow(ref$grid), paste(label, "grid size"), 0)
    equal(obj$params$cbar, ref$cbar, paste(label, "GLS localizing constant"), 0)
    equal(obj$params$separation, floor(trim * length(series)),
          paste(label, "minimum separation"), 0)
    equal(obj$params$n.obs, length(series), paste(label, "sample size"), 0)
    equal(obj$params$breaks, m, paste(label, "requested break count"), 0)
    equal(obj$params$max.lag, max_lag, paste(label, "maximum lag"), 0)
    check(inherits(obj, "bt_mdfCHLT"), paste(label, "S3 class"))
    check(identical(obj$params$method, "CHLT2014-MDF"), paste(label, "method"))
    check(length(obj$ADF$cv) == 1L && is.na(obj$ADF$cv),
          paste(label, "no single-break critical value"))

    core <- mdf(series, breaks = m, trim = trim, max.lag = max_lag, gls.only = TRUE)
    equal(core$statistic, ref$statistic, paste(label, "MDF.mlt statistic"))
    equal(as.numeric(core$break.points), ref$dates, paste(label, "MDF.mlt dates"), 0)
    check(is.matrix(core$break.points), paste(label, "MDF.mlt date matrix"))
    equal(unname(core$grid), ref$grid, paste(label, "MDF.mlt full grid"), 0)
    equal(core$lag, ref$lag, paste(label, "MDF.mlt lag"), 0)
    equal(core$cbar, ref$cbar, paste(label, "MDF.mlt localizing constant"), 0)
    equal(core$max.lag, max_lag, paste(label, "MDF.mlt maximum lag"), 0)
  }
}

reference <- minimum_reference(y, 2L, trim)
plain <- uroot.CHLT(y, max.lag = 0L, trim = trim, breaks = 2L)
for (transformed in list(y + 80, y - 0.8 * time, 3.5 * y, -2.5 * y,
                        matrix(y, ncol = 1L), ts(y))) {
  obj <- uroot.CHLT(transformed, max.lag = 0L, trim = trim, breaks = 2L)
  equal(obj$ADF$statistic, plain$ADF$statistic, "affine/shape invariant statistic")
  equal(obj$params$break.points, plain$params$break.points,
        "affine/shape invariant minimizer", 0)
}

# Algorithm 1: raw second differences padded with two zeros, Gaussian
# multipliers, one integration, and a fresh exhaustive minimum with lag zero.
bootstrap_reference <- function(y, breaks, trim, iterations, seed) {
  innovations <- c(0, 0, diff(y, differences = 2L))
  set.seed(seed)
  multipliers <- matrix(rnorm(length(y) * iterations), nrow = length(y))
  statistics <- vapply(seq_len(iterations), function(i) {
    minimum_reference(cumsum(innovations * multipliers[, i]), breaks, trim)$statistic
  }, numeric(1L))
  list(innovations = innovations, statistics = statistics, rng = .Random.seed)
}

b <- 9L
seed <- 19304L
boot_ref <- bootstrap_reference(y, 2L, trim, b, seed)
set.seed(seed)
inline <- uroot.CHLT(y, max.lag = 0L, trim = trim, breaks = 2L,
                     boot.cv = TRUE, boot.iter = b)
equal(inline$params$bootstrap$innovations, boot_ref$innovations,
      "bootstrap raw second differences and padding", 0)
equal(inline$params$bootstrap$statistics, boot_ref$statistics,
      "seeded bootstrap exhaustive statistics")
equal(.Random.seed, boot_ref$rng, "bootstrap consumes exactly N times B multipliers", 0)
equal(inline$params$bootstrap$difference.order, 2L, "double-difference diagnostic", 0)
equal(inline$params$bootstrap$lag, 0L, "bootstrap unaugmented DF diagnostic", 0)
equal(inline$params$bootstrap$iterations, b, "bootstrap iteration diagnostic", 0)
equal(inline$ADF$boot.cv, quantile(boot_ref$statistics, 0.05, type = 1, names = FALSE),
      "bootstrap empirical lower-tail critical value")
equal(inline$ADF$p.value, mean(boot_ref$statistics <= reference$statistic),
      "paper p-value without plus-one correction", 0)
equal(inline$ADF$statistic, plain$ADF$statistic, "bootstrap keeps observed statistic", 0)

set.seed(seed)
repeated <- uroot.CHLT(y, max.lag = 0L, trim = trim, breaks = 2L,
                       boot.cv = TRUE, boot.iter = b)
check(identical(inline, repeated), "bootstrap repeats exactly with the same seed")
set.seed(seed)
with_lags <- uroot.CHLT(y, max.lag = 2L, trim = trim, breaks = 2L,
                        boot.cv = TRUE, boot.iter = b)
equal(with_lags$params$bootstrap$statistics, boot_ref$statistics,
      "bootstrap lag zero despite positive original maximum lag")
equal(with_lags$ADF$statistic, minimum_reference(y, 2L, trim, 2L)$statistic,
      "original statistic retains OLS-MAIC during bootstrap")

boot_ref3 <- bootstrap_reference(y3, 3L, trim, 5L, seed)
set.seed(seed)
inline3 <- uroot.CHLT(y3, max.lag = 0L, trim = trim, breaks = 3L,
                      boot.cv = TRUE, boot.iter = 5L)
equal(inline3$params$bootstrap$statistics, boot_ref3$statistics,
      "three-break bootstrap scans all dates with cbar 25.5")

set.seed(seed)
deferred <- bootstrap(plain, iter = b, y = y)
equal(deferred$ADF, inline$ADF$boot.cv, "standalone/integrated bootstrap critical value", 0)
equal(deferred$p.value, inline$ADF$p.value, "standalone/integrated bootstrap p-value", 0)
equal(deferred$statistics, inline$params$bootstrap$statistics,
      "standalone/integrated bootstrap statistics", 0)
equal(deferred$innovations, inline$params$bootstrap$innovations,
      "standalone/integrated bootstrap innovations", 0)
for (obj in list(plain, inline, inline3)) {
  printed <- capture.output(print(obj))
  check(length(printed) > 0L, "multiple-break print dispatch")
}
fails(bootstrap(plain, iter = b, y = y[-1L]), "bootstrap rejects wrong data length")
fails(bootstrap(plain, iter = 0L, y = y), "bootstrap rejects zero iterations")

baseline <- list(y = y, max.lag = 0L, trim = trim, breaks = 2L)
invalid <- list(
  missing_observations = list(y = replace(y, 3L, NA_real_)),
  infinite_observations = list(y = replace(y, 3L, Inf)),
  character_series = list(y = as.character(y)),
  multivariate_series = list(y = cbind(y, y)),
  short_series = list(y = y[seq_len(10L)]),
  zero_breaks = list(breaks = 0),
  four_breaks = list(breaks = 4),
  fractional_breaks = list(breaks = 2.5),
  missing_breaks = list(breaks = NA_real_),
  vector_breaks = list(breaks = c(2, 3)),
  zero_trim = list(trim = 0),
  impossible_grid = list(trim = 0.4, breaks = 3L),
  negative_lag = list(max.lag = -1),
  fractional_lag = list(max.lag = 1.5),
  missing_lag = list(max.lag = NA_real_),
  excessive_lag = list(max.lag = length(y)),
  zero_iterations = list(boot.cv = TRUE, boot.iter = 0),
  fractional_iterations = list(boot.cv = TRUE, boot.iter = 2.5),
  missing_boot_flag = list(boot.cv = NA)
)
for (name in names(invalid)) {
  arguments <- baseline
  arguments[names(invalid[[name]])] <- invalid[[name]]
  fails(do.call(uroot.CHLT, arguments), paste("rejects", name))
}
for (z in list(rep(1, 30), seq_len(30), 1e6 + seq_len(30),
               0.2 * seq_len(30) + pmax(seq_len(30) - 9L, 0) -
                 0.5 * pmax(seq_len(30) - 20L, 0))) {
  fails(uroot.CHLT(z, max.lag = 0L, breaks = 2L), "rejects exact deterministic series")
}

# Source the pre-change code into a separate environment. In particular, do
# not let the original function call the currently registered S3 bootstrap.
fixture_dir <- if (dir.exists("fixtures")) "fixtures" else file.path("tests", "fixtures")
legacy <- new.env(parent = bt)
for (file in c("original-CHLT.R", "original-bootstrap-CHLT.R",
                "original-print-CHLT.R", "original-MDF-mlt.R")) {
  sys.source(file.path(fixture_dir, file), envir = legacy)
}
legacy$bootstrap <- legacy$bootstrap.bt_mdfCHLT
outcome <- function(fun, arguments, seed) {
  set.seed(seed)
  warnings <- character()
  result <- withCallingHandlers(
    tryCatch(list(value = do.call(fun, arguments)),
             error = function(e) list(error = conditionMessage(e), class = class(e))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings, rng = .Random.seed)
}

set.seed(103)
single_y <- cumsum(rnorm(80L))
single_break_y <- single_y + 8 * pmax(seq_len(80L) - 40L, 0)
selected_break <- FALSE
single_cases <- list(
  list(y = single_y, max.lag = 0L),
  list(y = matrix(single_y), max.lag = 2L),
  list(y = single_break_y, max.lag = 0L),
  list(y = single_y, max.lag = 0L, boot.cv = TRUE, boot.iter = 21L),
  list(y = single_break_y, max.lag = 0L, boot.cv = TRUE, boot.iter = 21L),
  list(y = single_y, max.lag = 0L, boot.cv = TRUE, boot.iter = 3L),
  list(y = single_y, max.lag = 0L, trim = 0.1),
  list(y = replace(single_y, 3L, NA_real_), max.lag = 0L),
  list(y = single_y, max.lag = 0L, boot.cv = NA)
)
for (i in seq_along(single_cases)) {
  args <- single_cases[[i]]
  expected <- outcome(legacy$uroot.CHLT, args, 900L + i)
  current <- outcome(uroot.CHLT, args, 900L + i)
  explicit <- outcome(uroot.CHLT, c(args, list(breaks = 1L)), 900L + i)
  check(identical(current, expected), paste("single-break legacy output/error and RNG", i))
  check(identical(explicit, expected), paste("explicit single-break legacy behavior", i))
  old_obj <- expected$result$value
  if (!is.null(old_obj)) {
    old_trim <- if (is.null(args$trim)) 0.15 else args$trim
    selected_break <- selected_break || old_obj$params$ADF$tau >= old_trim
    check(identical(capture.output(print(old_obj)),
                    capture.output(legacy$print.bt_mdfCHLT(old_obj))),
          paste("single-break print unchanged", i))
  }
}
check(selected_break, "legacy fixtures exercise an endogenously selected trend break")

# Also exercise both fixed legacy bootstrap branches, including its original
# vector-input error. The baseline is the original method, not a repaired one.
single_obj <- legacy$uroot.CHLT(single_y, max.lag = 0L)
forced <- single_obj
forced$params$MZ <- list(tau = 0.5, cbar = -17.6)
forced$params$ADF <- list(tau = 0.5, cbar = -17.6)
forced$params$tb <- 40L
for (obj in list(single_obj, forced)) {
  for (series in list(single_y, matrix(single_y))) {
    args <- list(obj = obj, iter = 21L, y = series)
    expected <- outcome(legacy$bootstrap.bt_mdfCHLT, args, 407L)
    current <- outcome(bootstrap, args, 407L)
    check(identical(current, expected), "standalone single-break bootstrap and RNG unchanged")
  }
}

# Extending MDF.mlt must leave its full legacy mode unchanged by default.
for (m in 2:3) {
  args <- list(y = y, breaks = m, trim = 0.24)
  expected <- outcome(legacy$MDF.mlt, args, 301L)
  current <- outcome(mdf, args, 301L)
  check(!is.null(expected$result$value), paste("ordinary MDF.mlt fixture succeeds", m))
  check(identical(current, expected), paste("ordinary MDF.mlt legacy behavior", m))
}

# The new explicit lag override is validated in the full MDF.mlt mode too.
# Fractional lags must not silently produce fractional/inconsistent indexing.
for (invalid_lag in list(1.5, -1, NA_real_, Inf)) {
  fails(mdf(y, breaks = 2L, trim = 0.24, max.lag = invalid_lag,
            gls.only = FALSE),
        paste("ordinary MDF.mlt rejects explicit max.lag", invalid_lag))
}
full_with_lag <- mdf(y, breaks = 2L, trim = 0.24, max.lag = 1L,
                     gls.only = FALSE)
equal(full_with_lag$MDF.GLS.2$statistic,
      minimum_reference(y, 2L, 0.24, max.lag = 1L)$statistic,
      "ordinary MDF.mlt honors valid explicit maximum lag")

cat("CHLT multiple-break, independent-bootstrap and legacy checks passed:",
    length(checks), "\n")
