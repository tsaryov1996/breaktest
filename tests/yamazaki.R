library(breaktest)

# Independent deterministic checks for the implementation of Yamazaki (2019).
# Formula checks below deliberately do not call package numerical helpers.
n_checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
  n_checks <<- n_checks + 1L
}
equal <- function(a, b, label, tol = 2e-8) {
  check(isTRUE(all.equal(unname(a), unname(b), tolerance = tol,
                         check.attributes = FALSE)), label)
}
errors <- function(expr, label) {
  err <- tryCatch({force(expr); NULL}, error = identity)
  check(inherits(err, "error"), label)
}
inverse_root <- function(a) {
  eig <- eigen(a, symmetric = TRUE)
  eig$vectors %*% diag(1 / sqrt(eig$values), nrow(a)) %*% t(eig$vectors)
}
rotation <- function(f, root) {
  v <- drop(root %*% f)
  if (length(v) == 1L) return(matrix(-sign(v), 1L))
  stopifnot(length(v) == 2L)
  v <- v / sqrt(sum(v * v))
  target <- rep(-1 / sqrt(2), 2L)
  # A direct planar orthogonal map, independent of the Gram-Schmidt code.
  cbind(target, c(-target[2], target[1])) %*%
    t(cbind(v, c(-v[2], v[1])))
}
logmean <- function(v) {
  m <- max(v)
  m + log(mean(exp(v - m)))
}
independent_at <- function(y, x, z, dates, b, estimated, f, trim = .05) {
  n <- length(y)
  k <- ncol(x)
  design <- cbind(x, x * as.numeric(seq_len(n) > b), z)
  e <- lm.fit(design, y)$residuals
  augmented <- design
  if (abs(b - estimated) >= k) {
    augmented <- cbind(design, x * as.numeric(seq_len(n) > estimated))
  }
  ea <- lm.fit(augmented, y)$residuals
  cum_orig <- cum_mod <- cum_orig_raw <- cum_mod_raw <- matrix(0, n, k)
  for (seg in seq_len(2L)) {
    ii <- if (seg == 1L) seq_len(b) else seq.int(b + 1L, n)
    raw <- x[ii, , drop = FALSE] * e[ii]
    aug <- x[ii, , drop = FALSE] * ea[ii]
    roots <- list(inverse_root(crossprod(raw) / length(ii)),
                  inverse_root(crossprod(aug) / length(ii)))
    for (kind in seq_len(2L)) {
      rr <- roots[[kind]]
      cc <- rotation(f[, seg], rr)
      total <- numeric(k)
      for (i in ii) {
        total <- total + x[i, ] * e[i]
        standardized <- drop(rr %*% total)
        directed <- drop(cc %*% standardized)
        if (kind == 1L) {
          cum_orig[i, ] <- directed
          cum_orig_raw[i, ] <- standardized
        } else {
          cum_mod[i, ] <- directed
          cum_mod_raw[i, ] <- standardized
        }
      }
    }
  }
  aggregate <- function(v, original = FALSE) {
    power <- if (original) rowSums(v^2) else rowSums(abs(v) * v)
    sum(power[seq_len(b)]) / b^2 + sum(power[seq.int(b + 1L, n)]) / (n - b)^2
  }
  ff <- rowSums(abs(cum_mod[dates, , drop = FALSE]) * cum_mod[dates, , drop = FALSE]) / n
  fo <- rowSums(cum_mod_raw[dates, , drop = FALSE]^2) / n
  ff[dates == b] <- fo[dates == b] <- 0
  outside <- abs(dates - b) > floor(trim * n)
  d <- dates[outside] / n
  lambda <- b / n
  omega <- ifelse(d < lambda, d * (lambda - d) / lambda,
                  (d - lambda) * (1 - d) / (1 - lambda))
  list(
    proposed = c(EM = aggregate(cum_orig), MEM = aggregate(cum_mod),
                 SUP = max(ff), AVG = mean(ff),
                 EXP = logmean(-k / 2 * log(omega) + ff[outside] / (2 * omega))),
    original = c(EM = aggregate(cum_orig_raw, TRUE),
                 MEM = aggregate(cum_mod_raw, TRUE), SUP = max(fo), AVG = mean(fo),
                 EXP = logmean(-k / 2 * log(omega) + fo[outside] / (2 * omega)))
  )
}

set.seed(3380913)
n <- 160L
tt <- seq_len(n)
x1 <- matrix(1, n, 1L)
y1 <- .5 + .8 * (tt > 80L) + rnorm(n)
one <- conf.sets.Yamazaki(y1, direction = 1)
check(identical(colnames(one$stats), c("EM", "MEM", "SUP", "AVG", "EXP")),
      "Statistic labels")
check(identical(as.integer(one$candidates), 8L:152L), "Original-row candidate indices")
check(identical(dim(one$masks), c(n, 5L)), "Original-row mask dimensions")
check(all(one$masks[-one$candidates, ] == 0), "No acceptance outside candidate set")
check(all(one$masks[one$candidates, ] == (one$stats <= one$cv)), "Set inversion")
check(all(is.finite(one$stats)) && all(is.finite(one$cv)), "Finite statistics and critical values")
for (b in c(40L, one$estimate, 120L)) {
  expected <- independent_at(y1, x1, NULL, one$candidates, b, one$estimate, matrix(1, 1, 2))
  at <- match(b, one$candidates)
  equal(one$stats[at, ], expected$proposed, paste("Scalar independent formulas at", b))
  equal(one$original.stats[at, ], expected$original, paste("Scalar original formulas at", b))
}

x2 <- cbind(1, rnorm(n))
z2 <- matrix(rnorm(n), ncol = 1)
y2 <- drop(x2 %*% c(.2, -.5) + z2 * .4 +
           (tt > 80L) * (x2 %*% c(1, .5))) + rnorm(n)
two <- conf.sets.Yamazaki(y2, x2, z2, direction = c(1, .5))
for (b in c(35L, two$estimate, 120L)) {
  expected <- independent_at(y2, x2, z2, two$candidates, b, two$estimate,
                             matrix(c(1, .5), 2L, 2L))
  at <- match(b, two$candidates)
  equal(two$stats[at, ], expected$proposed, paste("k=2 independent formulas at", b))
  equal(two$original.stats[at, ], expected$original, paste("k=2 original formulas at", b))
}

scaled <- conf.sets.Yamazaki(y2 * 7, x2, z2, direction = 3 * c(1, .5))
equal(scaled$stats, two$stats, "Positive outcome/direction scaling invariance")
negative <- conf.sets.Yamazaki(-y2, x2, z2, direction = -c(1, .5))
equal(negative$stats, two$stats, "Joint response/direction sign reversal")
opposite <- conf.sets.Yamazaki(y2, x2, z2, direction = -c(1, .5))
equal(opposite$original.stats, two$original.stats, "Original methods ignore direction")
equal(opposite$stats[, c("EM", "MEM", "AVG")],
      -two$stats[, c("EM", "MEM", "AVG")], "Signed directional quadratic scores")
zfirst <- conf.sets.Yamazaki(y2, x2, z2, direction = c(0, 1))
check(all(is.finite(zfirst$stats)), "Zero first direction component remains valid")
byregime <- conf.sets.Yamazaki(y2, x2, z2,
                              direction = cbind(c(1, .5), c(-.3, 1)))
check(all(is.finite(byregime$stats)), "Regime-specific supplied directions")

# Procedure B: unbroken regression, signed integrated squares, signed square roots.
e0 <- lm.fit(cbind(x2, z2), y2)$residuals
partial <- apply(x2 * e0, 2L, cumsum)
btilde <- -colSums(abs(partial) * partial)
fhat <- sign(btilde) * sqrt(abs(btilde))
automatic <- conf.sets.Yamazaki(y2, x2, z2)
explicit <- conf.sets.Yamazaki(y2, x2, z2, direction = fhat)
equal(automatic$stats, explicit$stats, "Procedure B independently reconstructed")
automatic_neg <- conf.sets.Yamazaki(-y2, x2, z2)
equal(automatic_neg$stats, automatic$stats, "Estimated-direction sign reversal")

low <- conf.sets.Yamazaki(y1, conf.level = .90, direction = 1)
equal(low$stats, one$stats, "Confidence level does not change statistics")
check(all(low$masks <= one$masks), "90 percent sets nested in 95 percent sets")
cv_helper <- getFromNamespace(".yamazaki.cv", "breaktest")
legacy_cv <- getFromNamespace("get.cv.KY", "breaktest")
equal(cv_helper(.5, .95, 1), c(.493, .493, .881, .135, 2.670),
      "95 percent k=1 midpoint table coefficients")
equal(cv_helper(.5, .90, 1), c(.350, .350, .709, .095, 2.220),
      "90 percent k=1 midpoint table coefficients")
for (level in c(.90, .95)) for (k in 1L:9L) for (lambda in c(.05, .1, .5, .9, .95)) {
  old <- legacy_cv(lambda, level, k)
  equal(cv_helper(lambda, level, k, modified = FALSE),
        c(old$cval_em, old$cval_em, old$cval_sup, old$cval_avg, old$cval_exp),
        "Original critical values preserve KY table normalization")
  check(all(is.finite(cv_helper(lambda, level, k))), "Directional critical values finite on supported grid")
}
hac <- conf.sets.Yamazaki(y1, direction = 1, covariance = "HAC")
check(all(is.finite(hac$stats)), "HAC path returns finite statistics")

# Independent scalar QS-HAC reference: acf() provides all autocovariances with
# denominator m. This also checks prewhitening/recolouring and bandwidth units.
score <- as.numeric(stats::filter(rnorm(180), .55, method = "recursive"))
ar1 <- sum(head(score, -1) * tail(score, -1)) / sum(head(score, -1)^2)
innovation <- tail(score, -1) - ar1 * head(score, -1)
m <- length(innovation)
aux <- sum(head(innovation, -1) * tail(innovation, -1)) / sum(head(innovation, -1)^2)
aux <- min(.97, max(-.97, aux))
bandwidth <- 1.3221 * (m * 4 * aux^2 / (1 - aux)^4)^.2
auto <- as.numeric(stats::acf(innovation, type = "covariance", demean = FALSE,
                              plot = FALSE, lag.max = m - 1L)$acf)
frequency <- 6 * pi * seq_len(m - 1L) / (5 * bandwidth)
weights <- 3 / frequency^2 * (sin(frequency) / frequency - cos(frequency))
reference_hac <- (auto[1L] + 2 * sum(weights * auto[-1L])) / (1 - ar1)^2
actual_hac <- getFromNamespace(".yamazaki.cov", "breaktest")(matrix(score), "HAC")
equal(actual_hac, matrix(reference_hac), "Independent scalar prewhitened QS-HAC")

# A non-diagonal VAR coefficient makes transpose errors in recolouring visible.
score2 <- cbind(score[1:160], .4 * score[1:160] + rnorm(160))
past <- score2[-nrow(score2), , drop = FALSE]
future <- score2[-1L, , drop = FALSE]
coefficient <- solve(crossprod(past), crossprod(past, future))
innov2 <- future - past %*% coefficient
m2 <- nrow(innov2)
rho <- sigma2 <- numeric(2L)
for (j in 1L:2L) {
  v <- innov2[, j]
  rho[j] <- min(.97, max(-.97, sum(head(v, -1) * tail(v, -1)) / sum(head(v, -1)^2)))
  sigma2[j] <- mean((tail(v, -1) - rho[j] * head(v, -1))^2)
}
alpha <- sum(4 * rho^2 * sigma2^2 / (1 - rho)^8) /
  sum(sigma2^2 / (1 - rho)^4)
bw2 <- 1.3221 * (m2 * alpha)^.2
# A single symmetric kernel matrix computes the double covariance sum directly.
distance <- abs(outer(seq_len(m2), seq_len(m2), "-"))
theta <- 6 * pi * distance / (5 * bw2)
kernel <- matrix(1, m2, m2)
nonzero <- distance > 0
kernel[nonzero] <- 3 / theta[nonzero]^2 *
  (sin(theta[nonzero]) / theta[nonzero] - cos(theta[nonzero]))
omega_innov <- crossprod(innov2, kernel %*% innov2) / m2
left <- solve(diag(2) - t(coefficient))
expected2 <- left %*% omega_innov %*% t(left)
actual2 <- getFromNamespace(".yamazaki.cov", "breaktest")(score2, "HAC")
equal(actual2, expected2, "Independent multivariate HAC covariance orientation")

errors(conf.sets.Yamazaki(rep(0, n)), "Zero variance rejected")
errors(conf.sets.Yamazaki(rep(1, n)), "Nonzero constant response rejected")
errors(conf.sets.Yamazaki(c(rep(0, n / 2), rep(1, n / 2))),
       "Exact noiseless break rejected")
errors(conf.sets.Yamazaki(y1, cbind(x1, x1)), "Rank-deficient regression rejected")
errors(conf.sets.Yamazaki(y1, direction = 0), "Zero direction rejected")
errors(conf.sets.Yamazaki(y1, direction = NA_real_), "Nonfinite direction rejected")
errors(conf.sets.Yamazaki(y1, direction = c(1, 1)), "Direction dimension checked")
errors(conf.sets.Yamazaki(replace(y1, 5L, NA_real_)), "Missing values rejected without reindexing")
errors(conf.sets.Yamazaki(replace(y1, 5L, Inf)), "Infinite response rejected")
errors(conf.sets.Yamazaki(y1, x1[-1, , drop = FALSE]), "Regressor length checked")
errors(conf.sets.Yamazaki(y1, conf.level = .80), "Unsupported level rejected")
errors(conf.sets.Yamazaki(y1, conf.level = c(.90, .95)), "Vector level rejected")
errors(conf.sets.Yamazaki(y1, trim = .10), "Unsupported critical-value trim rejected")
errors(conf.sets.Yamazaki(y1, covariance = "wrong"), "Covariance choice checked")
tiny <- conf.sets.Yamazaki(y1 * 1e-12, direction = 1)
equal(tiny$stats, one$stats, "Small nonzero outcome scale handled relatively", tol = 1e-6)

cat(sprintf("Yamazaki independent checks passed: %d\n", n_checks))
