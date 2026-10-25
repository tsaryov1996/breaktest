# Internal implementation of CHLT (January 2014), Algorithm 1 (double-difference wild
# bootstrap). The original sample uses the OLS-MAIC convention of MDF.mlt;
# bootstrap samples use unaugmented DF statistics, as required by Algorithm 1.
# This is the infimum over exactly m candidate trend hinges, allowing their
# coefficients to be zero. It is not a minimum over the tests for 0,...,m.

.chlt.validate <- function(y, max.lag, trim, breaks) {
  if (!is.numeric(y) || !is.null(dim(y)) &&
      (length(dim(y)) != 2L || ncol(y) != 1L) ||
      length(y) < 12L || any(!is.finite(y))) {
    stop("y must be a finite numeric vector or one-column matrix of length >= 12")
  }
  if (length(breaks) != 1L || !is.numeric(breaks) ||
      !is.finite(breaks) || !breaks %in% 1:3) {
    stop("breaks must be one of 1, 2, or 3")
  }
  N <- length(y)
  if (!any(diff(as.numeric(y), differences = 2L) != 0)) {
    stop("Degenerate series: y is an exact constant or linear trend")
  }
  if (length(trim) != 1L || !is.numeric(trim) || !is.finite(trim) ||
      trim <= 0 || trim >= 0.5 || floor(trim * N) < 2L) {
    stop("trim must be in (0, 0.5), with floor(trim * length(y)) >= 2")
  }
  max.lag <- .chlt.max.lag(N, max.lag)
  if (floor((1 - trim) * N) - floor(trim * N) <
      (breaks - 1L) * floor(trim * N)) {
    stop("No admissible break configuration: reduce breaks or trim")
  }
  list(y = matrix(as.numeric(y), ncol = 1L), max.lag = as.integer(max.lag))
}

.chlt.max.lag <- function(N, max.lag) {
  if (is.null(max.lag)) max.lag <- trunc(12 * (N / 100)^(1 / 4))
  if (length(max.lag) != 1L || !is.numeric(max.lag) ||
      !is.finite(max.lag) || max.lag < 0 || max.lag != floor(max.lag) ||
      N - 2 * max.lag - 2 <= 0) {
    stop("max.lag must be a nonnegative integer leaving positive ADF residual degrees of freedom")
  }
  as.integer(max.lag)
}

.chlt.validate.iter <- function(iter) {
  if (length(iter) != 1L || !is.numeric(iter) || !is.finite(iter) ||
      iter < 1 || iter != floor(iter)) {
    stop("The number of bootstrap iterations must be a positive integer")
  }
}

# A deterministic, lexicographic exhaustive grid. The grid and coefficients
# are reused across all bootstrap replications; no break is held fixed at its
# estimate from the observed sample.
.chlt.break.grid <- function(N, breaks, trim) {
  lower <- floor(trim * N)
  upper <- floor((1 - trim) * N)
  separation <- max(1L, lower)
  visit <- function(prefix, first, remaining) {
    last <- upper - (remaining - 1L) * separation
    if (first > last) return(NULL)
    if (remaining == 1L) {
      tail <- seq.int(first, last)
      return(cbind(matrix(rep(prefix, each = length(tail)),
                          nrow = length(tail)), tail))
    }
    do.call(rbind, lapply(seq.int(first, last), function(b) {
      visit(c(prefix, b), b + separation, remaining - 1L)
    }))
  }
  result <- visit(integer(), lower, breaks)
  if (is.null(result)) stop("No admissible break configuration")
  storage.mode(result) <- "integer"
  colnames(result) <- paste0("break", seq_len(breaks))
  result
}

# Evaluate an unaugmented DF t-ratio columnwise, with N-2 residual df.
.chlt.df0 <- function(residuals) {
  N <- nrow(residuals)
  lagged <- residuals[-N, , drop = FALSE]
  delta <- residuals[-1L, , drop = FALSE] - lagged
  xx <- colSums(lagged^2)
  xy <- colSums(lagged * delta)
  coefficient <- xy / xx
  errors <- delta - sweep(lagged, 2L, coefficient, "*")
  rss <- colSums(errors^2)
  value <- xy / sqrt(xx * rss / (N - 2L))
  if (any(!is.finite(value)) || any(xx <= 0) || any(rss <= 0)) {
    stop("Degenerate series: the DF statistic is undefined")
  }
  value
}

.chlt.mdf.scan <- function(y, grid, cbar, max.lag = 0L) {
  y <- as.matrix(y)
  N <- nrow(y)
  B <- ncol(y)
  if (B > 1L && max.lag != 0L) stop("Augmented scans require a single series")
  rho <- 1 - cbar / N
  quasi.y <- rbind(y[1L, , drop = FALSE],
                   y[-1L, , drop = FALSE] - rho * y[-N, , drop = FALSE])
  tr <- seq_len(N)
  candidate.dates <- sort(unique(as.vector(grid)))
  hinges <- outer(tr, candidate.dates, function(t, b) pmax(t - b, 0))
  best <- rep(Inf, B)
  best.row <- rep(NA_integer_, B)
  best.lag <- integer(B)

  for (i in seq_len(nrow(grid))) {
    z <- cbind(1, tr, hinges[, match(grid[i, ], candidate.dates), drop = FALSE])
    quasi.z <- rbind(z[1L, , drop = FALSE],
                     z[-1L, , drop = FALSE] - rho * z[-N, , drop = FALSE])
    fit <- qr(quasi.z)
    if (fit$rank != ncol(z)) stop("Singular GLS deterministic design")
    residuals <- y - z %*% qr.coef(fit, quasi.y)
    if (any(sqrt(colSums(residuals^2)) <=
            100 * .Machine$double.eps * sqrt(colSums(y^2)))) {
      stop("Degenerate series: a candidate broken trend fits y to numerical precision")
    }
    lag <- 0L
    if (max.lag == 0L) {
      statistic <- .chlt.df0(residuals)
    } else {
      # Preserve the existing MDF.mlt lag-selection convention. The final ADF
      # regression is refit over lag+2:N after selection on a common sample.
      ols.residuals <- qr.resid(qr(z), y)
      lag <- uroot.ADF(ols.residuals, const = FALSE, trend = FALSE,
                       max.lag = max.lag, criterion = "aic",
                       modified.criterion = TRUE)$lag
      statistic <- uroot.ADF(residuals, const = FALSE, trend = FALSE,
                             max.lag = lag, criterion = NULL)$t.alpha
      if (!is.finite(statistic)) stop("Non-finite augmented DF statistic")
    }
    improve <- statistic < best
    best[improve] <- statistic[improve]
    best.row[improve] <- i
    best.lag[improve] <- lag
  }
  list(statistic = unname(best), break.points = grid[best.row, , drop = FALSE],
       lag = best.lag)
}

.chlt.multiple <- function(y, max.lag, trim, breaks, boot.cv, boot.iter) {
  if (!is.logical(boot.cv) || length(boot.cv) != 1L || is.na(boot.cv)) {
    stop("boot.cv must be TRUE or FALSE")
  }
  if (boot.cv) .chlt.validate.iter(boot.iter)
  # Use the requested GLS component of the existing multiple-break procedure,
  # without computing its unrelated OLS, pretest and union statistics.
  fit <- MDF.mlt(y, const = FALSE, breaks = breaks, trim = trim,
                 max.lag = max.lag, gls.only = TRUE)
  grid <- fit$grid
  N <- fit$n.obs
  result <- list(
    ADF = list(statistic = fit$statistic[1L], cv = NA_real_),
    params = list(method = "CHLT2014-MDF", breaks = breaks,
                  n.obs = N, trim = trim, separation = floor(trim * N),
                  max.lag = fit$max.lag, lag = fit$lag[1L],
                  break.points = unname(fit$break.points[1L, ]),
                  cbar = fit$cbar, configurations = nrow(grid), level = 0.05)
  )
  class(result) <- "bt_mdfCHLT"
  if (boot.cv) {
    boot <- .chlt.bootstrap.multiple(result, boot.iter, y, grid)
    result$ADF$boot.cv <- boot$ADF
    result$ADF$p.value <- boot$p.value
    result$params$bootstrap <- boot[c("iterations", "difference.order", "lag",
                                      "innovations", "statistics")]
  }
  result
}

.chlt.bootstrap.multiple <- function(obj, iter, y, grid = NULL) {
  .chlt.validate.iter(iter)
  checked <- .chlt.validate(y, obj$params$max.lag,
                            obj$params$trim, obj$params$breaks)
  y <- as.numeric(checked$y)
  N <- length(y)
  if (N != obj$params$n.obs) stop("y must have the length used for the original test")
  if (is.null(grid)) grid <- .chlt.break.grid(N, obj$params$breaks, obj$params$trim)

  # Algorithm 1: no residual resampling, centering, sqrt(2) adjustment, second
  # integration, fitted trend, or AR sieve. Scale cancels in the DF t-ratio.
  innovations <- c(0, 0, diff(y, differences = 2L))
  if (!any(innovations != 0)) stop("Second differences are all zero")
  pseudo <- matrix(rnorm(N * iter), nrow = N, ncol = iter) * innovations
  for (i in seq_len(iter)) pseudo[, i] <- cumsum(pseudo[, i])
  simulated <- .chlt.mdf.scan(pseudo, grid, obj$params$cbar, max.lag = 0L)$statistic
  list(ADF = unname(stats::quantile(simulated, 0.05, type = 1)),
       p.value = mean(simulated <= obj$ADF$statistic),
       iterations = iter, difference.order = 2L, lag = 0L,
       innovations = innovations, statistics = simulated)
}
