library(breaktest)

expect.equal <- function(x, y, tolerance = 1e-9) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
}
expect.error <- function(expr, message) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(message, conditionMessage(error), fixed = TRUE))
}

# Independent level regressions, fitted by QR rather than the package helper.
reference.dates <- function(y, trim, weights = rep(1, length(y) - 1L)) {
  N <- length(y) - 1L
  h <- floor(N * trim)
  fit.window <- function(a, b) {
    lower <- a + max(h + 1L, 2L)
    upper <- b - max(h, 2L)
    if (lower > upper) return(list(date = NA_integer_, profile = data.frame(
      date = integer(), objective = numeric())))
    dates <- seq.int(lower, upper)
    ii <- seq.int(a + 1L, b)
    rss <- vapply(dates, function(k) {
      X <- cbind(y[ii] * (ii <= k), y[ii] * (ii > k))
      fit <- lm.wfit(X, y[ii + 1L], w = weights[ii])
      if (fit$rank < 2L) return(Inf)
      sum(weights[ii] * fit$residuals^2)
    }, numeric(1))
    list(date = dates[which.min(rss)], profile = data.frame(date = dates, objective = rss))
  }
  collapse <- fit.window(0L, N)
  emerge <- fit.window(0L, collapse$date)
  recovery <- fit.window(collapse$date, N)
  list(dates = c(emerge$date, collapse$date, recovery$date),
       search = list(emerge$profile, collapse$profile, recovery$profile))
}
kernel.variance <- function(residuals, bandwidth) {
  distance <- outer(seq_along(residuals), seq_along(residuals), "-")
  K <- exp(-0.5 * (distance / (length(residuals) * bandwidth))^2)
  diag(K) <- 0
  drop(K %*% residuals^2) / rowSums(K)
}

set.seed(42)
N <- 200L
y <- numeric(N + 1L)
y[1L] <- 100
phi <- rep(c(1, 1.04, 0.94, 1), c(70, 50, 40, 40))
sigma <- rep(c(3, 0.7), c(50, 150))
for (t in seq_len(N)) y[t + 1L] <- phi[t] * y[t] + sigma[t] * rnorm(1L)

ols <- bubble.dates.NBCN(y, method = "ols")
reference.ols <- reference.dates(y / max(abs(y)), .05)
expect.equal(ols$break.points, reference.ols$dates + 1L)
for (j in seq_len(3L)) expect.equal(ols$search[[j]], reference.ols$search[[j]])
stopifnot(ols$status == "complete", is.null(ols$variance), is.null(ols$cv))

# Equation (4) and Gaussian leave-one-out smoothing are checked independently.
grid <- c(.2, .04, .1)
wls <- bubble.dates.NBCN(y, bandwidth.grid = grid)
t <- seq_len(N)
explosive <- as.numeric(t > ols$Te_est & t <= ols$Tc_est)
collapsing <- as.numeric(t > ols$Tc_est & t <= ols$Tr_est)
pilot <- lm.fit(cbind(explosive, collapsing,
                     explosive * y[-length(y)], collapsing * y[-length(y)]), diff(y))
expect.equal(wls$pilot$residuals, pilot$residuals)
expect.equal(wls$pilot$coefficients, pilot$coefficients)
outside <- t <= ols$Te_est | t > ols$Tr_est
expect.equal(wls$pilot$residuals[outside], diff(y)[outside])
scores <- vapply(sort(grid), function(b) {
  residuals <- pilot$residuals / max(abs(y))
  mean((residuals^2 - kernel.variance(residuals, b))^2)
}, numeric(1))
expect.equal(wls$cv, scores)
expect.equal(wls$bandwidth, sort(grid)[which.min(scores)])
variance <- kernel.variance(pilot$residuals, wls$bandwidth)
expect.equal(wls$variance, variance)
nearest <- bubble.dates.NBCN(y, bandwidth = 1e-300)
r <- pilot$residuals
adjacent.variance <- c(r[2L]^2, (r[3:N]^2 + r[1:(N - 2L)]^2) / 2, r[N - 1L]^2)
expect.equal(nearest$variance, adjacent.variance)
weights <- min(variance) / variance
reference.wls <- reference.dates(y / max(abs(y)), .05, weights)
expect.equal(wls$break.points, reference.wls$dates + 1L)
for (j in seq_len(3L)) expect.equal(wls$search[[j]], reference.wls$search[[j]])

# Recovery must include the equation immediately following collapse.
ii <- seq.int(wls$Tc_est + 1L, N)
k <- wls$search$recovery$date[1L]
X <- cbind(y[ii] * (ii <= k), y[ii] * (ii > k)) / wls$scale
fit <- lm.wfit(X, y[ii + 1L] / wls$scale, weights[ii])
expect.equal(wls$search$recovery$objective[1L], sum(weights[ii] * fit$residuals^2))
omitted <- lm.wfit(X[-1L, ], y[ii[-1L] + 1L] / wls$scale, weights[ii[-1L]])
stopifnot(abs(wls$search$recovery$objective[1L] -
              sum(weights[ii[-1L]] * omitted$residuals^2)) > 1e-12)

regime <- cut(t, c(0L, reference.wls$dates, N), labels = FALSE)
X <- sapply(seq_len(4L), function(j) y[-length(y)] * (regime == j))
expect.equal(wls$coefficients, lm.wfit(X, y[-1L], weights)$coefficients)
expect.equal(wls$break.points, wls$pilot$break.points)
for (factor in c(-1, 1e-100, 1e100)) {
  scaled <- bubble.dates.NBCN(y * factor, bandwidth.grid = grid)
  expect.equal(scaled$break.points, wls$break.points)
  expect.equal(scaled$variance / factor / factor, wls$variance)
}
expect.equal(bubble.dates.NBCN(matrix(y, ncol = 1), method = "ols")$break.points,
             ols$break.points)
rng <- .Random.seed
invisible(bubble.dates.NBCN(y, bandwidth = .1))
stopifnot(identical(.Random.seed, rng))

# Reproducible low-signal cases exercise both sides of an incomplete OLS fit
# and a WLS fit whose pilot is complete but whose weighted recovery is not.
set.seed(20261109)
innovations <- matrix(rnorm(100 * 400), nrow = 100)
bubble.sample <- function(innovations, sigma = rep(1, 400)) {
  y <- numeric(401)
  for (t in seq_len(400)) {
    phi <- if (t > 160 && t <= 240) 1 + 5 / 400 else
      if (t > 240 && t <= 280) 1 - 6 / 400 else 1
    drift <- if (t <= 160 || t > 280) 1 / 800 else 0
    y[t + 1L] <- drift + phi * y[t] + sigma[t] * innovations[t]
  }
  y
}
for (i in 1:2) {
  short <- bubble.sample(innovations[i, ])
  partial <- bubble.dates.NBCN(short, method = "ols")
  reference <- reference.dates(short / max(abs(short)), .05)
  expect.equal(partial$break.points, reference$dates + 1L)
  stopifnot(partial$status == "partial", !is.na(partial$Tc_est),
            sum(is.na(partial$break.points)) == 1L)
  unavailable <- bubble.dates.NBCN(short)
  stopifnot(unavailable$method == "wls", unavailable$status == "incomplete.pilot",
            all(is.na(unavailable$break.points)), all(is.na(unavailable$coefficients)),
            is.null(unavailable$variance), is.null(unavailable$bandwidth),
            is.null(unavailable$cv), is.null(unavailable$search))
  expect.equal(unavailable$pilot$break.points, partial$break.points)
}
short <- bubble.sample(innovations[10, ], rep(c(5, 1), c(80, 320)))
partial <- bubble.dates.NBCN(short)
stopifnot(partial$status == "partial", partial$pilot$status == "complete",
          partial$pilot$Tc_est == 80L, partial$Tc_est == 366L,
          is.na(partial$Tr_est), !is.na(partial$Te_est),
          all(is.finite(partial$coefficients[1:2])),
          all(is.na(partial$coefficients[3:4])))
reference <- reference.dates(short / max(abs(short)), .05,
                             min(partial$variance) / partial$variance)
expect.equal(partial$break.points, reference$dates + 1L)
for (j in seq_len(3L)) expect.equal(partial$search[[j]], reference$search[[j]])

# Exact ties are deterministic; invalid input remains an error.
ties <- bubble.dates.NBCN(c(2^(0:10), 2^(9:0)), method = "ols")
expect.equal(ties$break.points, c(3L, 11L, 13L))
for (bad in list(c(y, NA), c(y, Inf), y[1:8], matrix(y[1:200], 100, 2), y + 1i)) {
  expect.error(bubble.dates.NBCN(bad), "finite numeric values")
}
for (bad in list(0, .25, NA_real_, c(.05, .1), .05 + 1i)) {
  expect.error(bubble.dates.NBCN(y, trim = bad), "trim must")
}
for (bad in list(0, -1, NA_real_, c(.05, .1), .05 + 1i)) {
  expect.error(bubble.dates.NBCN(y, bandwidth = bad), "Bandwidths must")
}
expect.error(bubble.dates.NBCN(y, bandwidth = .1, bandwidth.grid = .1), "not both")
expect.error(bubble.dates.NBCN(y, method = "ols", bandwidth = .1), "require method")
expect.error(bubble.dates.NBCN(rep(1, 100)), "constant series")
expect.error(bubble.dates.NBCN(c(rep(1, 4), rep(2, 4), rep(1, 4))), "rank deficient")
expect.error(bubble.dates.NBCN(c(1, 1, 1, 2, 4, 2, 1, 1, 1), bandwidth = 1e-300),
             "strictly positive")
expect.error(bubble.dates.NBCN(y * 1e200, bandwidth = .1), "rescale y")
expect.error(bubble.dates.NBCN(y * 1e-200, bandwidth = .1), "rescale y")
