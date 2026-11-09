#' @title
#' Direction-aware confidence sets for one structural break
#'
#' @description
#' Inverts the signed EM, modified EM, and sup/avg/exp KY statistics proposed
#' by Yamazaki (2019), and returns the corresponding unsigned statistics and
#' confidence sets on the same sample for comparison.
#'
#' @param y A finite numeric response vector or one-column matrix. Missing
#' observations are rejected so that dates retain their original row indices.
#' @param x Regressors whose coefficients can change at the break. For the
#' estimation function, \code{NULL} means a single intercept; otherwise supply
#' all columns explicitly, including an intercept if required. For \code{print},
#' the fitted object.
#' @param z Optional regressors whose coefficients remain fixed, with the same
#' number of rows as \code{y}. No intercept is added automatically.
#' @param conf.level One of \code{0.90} and \code{0.95}.
#' @param direction \code{NULL} estimates the direction by Procedure B in the
#' paper. Alternatively, supply a nonzero numeric vector proportional to
#' \eqn{\Sigma_{xx}\delta}, up to a positive scalar, or a matrix with two columns
#' for the two regimes. With several regressors this is generally not the raw
#' coefficient-change vector \eqn{\delta}. With one regressor its sign specifies
#' an upward or downward coefficient change.
#' @param covariance \code{"iid"} uses empirical score second moments, suitable
#' for serially uncorrelated scores, allowing heteroskedasticity.
#' \code{"HAC"} uses a quadratic spectral estimator with Andrews bandwidth,
#' VAR(1) prewhitening and recolouring.
#' @param trim Only \code{0.05} is supported by the supplied response surfaces.
#'
#' @details
#' The model is \eqn{y_t=x_t'\beta+1\{t>b\}x_t'\delta+z_t'\gamma+u_t}.
#' The returned date is the last observation before the coefficient change.
#' An identifiable, single, nonzero break is assumed. This is inference on the
#' date conditional on a break, not a test that a break exists.
#'
#' The method requires stationary moment conditions. Unit-root regressors,
#' deterministic trends, cointegrating models, and multiple breaks are outside
#' the theory implemented here. The default direction estimator additionally
#' requires equal regressor second moments across the regimes. These statistical
#' assumptions cannot be verified by input validation.
#'
#' For multiple changing coefficients the asymptotic result also assumes a
#' convergent choice of orthogonal rotation (Assumption 3). Numerically completing
#' the basis at a zero first whitened-direction coordinate does not itself prove
#' this assumption for exceptional directions with several zero components when
#' \eqn{k\geq3}.
#'
#' The admissible grid is \code{ceiling(0.05*N):floor(0.95*N)}. Using a ceiling
#' at the lower edge prevents extrapolation of critical-value surfaces below
#' 0.05 for sample sizes not divisible by 20. Rank or covariance problems cause
#' an error. The exclusion gap for the EXP statistic is
#' \code{abs(T2-T1) > floor(0.05*N)}.
#'
#' Signed quadratic forms are \eqn{\sum_j |w_j|w_j}; negative statistics are
#' meaningful. EM and MEM use all observations in each regime. AVG includes the
#' zero contribution at \eqn{T_2=T_1} and divides by the full number of grid
#' points. EXP is calculated by a stable log-mean-exp formula.
#'
#' The test scores always come from the maintained-break regression. EM
#' covariances use those residuals; MEM and KY covariances additionally allow
#' the SSR-minimising break in the residual regression. The implementation
#' uses QR regression and a symmetric covariance inverse square root. Numerically
#' zero residual variances and singular covariance matrices are rejected.
#'
#' The HAC bandwidth's auxiliary AR coefficients are capped in absolute value
#' at 0.97 as a finite-sample safeguard. This convention may differ from the author's
#' GAUSS implementation. Fitted VAR coefficients are not clipped; invertible
#' unstable fits are counted in \code{hac.diagnostics}, as are auxiliary
#' coefficient caps. Singular recolouring fails with an error. Counts cover
#' all four covariance fits per candidate.
#'
#' Modified critical values reproduce Tables 1 and 2 of the May 2019 working
#' paper, for one to nine changing coefficients. Original critical values use
#' the existing KY tables in the package. Rounded response-surface coefficients
#' and finite-sample covariance estimates can affect coverage; nominal coverage
#' is not a finite-sample guarantee. EXP can under-cover in small samples, as
#' reported in the paper.
#'
#' @return An object of class \code{bt_confSetYamazaki}, containing:
#' \item{estimate}{SSR-minimising break date.}
#' \item{candidates}{All tested last-pre-break row indices.}
#' \item{stats, original.stats}{Matrices of signed and unsigned statistics;
#' columns are EM, MEM, SUP, AVG and EXP. Rows match \code{candidates}.}
#' \item{cv, original.cv}{Corresponding critical-value matrices.}
#' \item{csets, original.csets}{Named lists of all accepted row indices,
#' including disconnected components.}
#' \item{masks, original.masks}{Full-sample \eqn{N\times5} logical membership
#' matrices.}
#' \item{direction}{Direction used for the two regimes.}
#' \item{hac.diagnostics}{Counts of auxiliary AR caps and unstable fitted VARs;
#' zero for \code{covariance = "iid"}.}
#' Other fields record the model and covariance settings and diagnostics.
#'
#' @references
#' Yamazaki, D. (2019). Improved confidence sets for the date of a structural
#' break. May 2019 working paper. \doi{10.2139/ssrn.3380913}.
#'
#' Kurozumi, E. and Yamamoto, Y. (2015). Confidence sets for the break date based
#' on optimal tests. The Econometrics Journal, 18, 412--435.
#' \doi{10.1111/ectj.12055}.
#'
#' @seealso [plot_confsets], [fcst.conf.set.EoMorley]
#'
#' @examples
#' set.seed(3380913)
#' y <- rnorm(160) + 1.2 * (seq_len(160) > 80)
#' ans <- conf.sets.Yamazaki(y)
#' ans
#' ans$csets$SUP
#' ans$original.csets$SUP
#' plot_confsets(y, as.list(as.data.frame(ans$masks)),
#'               break_points = ans$estimate)
#'
#' @export
conf.sets.Yamazaki <- function(
  y,
  x = NULL,
  z = NULL,
  conf.level = 0.95,
  direction = NULL,
  covariance = c("iid", "HAC"),
  trim = 0.05
) {
  call <- match.call()
  covariance <- match.arg(covariance)
  if (!is.numeric(y) || is.complex(y) ||
      (!is.null(dim(y)) && (length(dim(y)) != 2L || ncol(y) != 1L)) ||
      !length(y) || any(!is.finite(y))) {
    stop("y must be a finite numeric vector or one-column matrix.")
  }
  time <- if (stats::is.ts(y)) as.numeric(stats::time(y)) else seq_along(y)
  y <- as.numeric(y)
  n <- length(y)
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) || trim != 0.05) {
    stop("The published critical values require trim = 0.05.")
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || !conf.level %in% c(0.90, 0.95)) {
    stop("conf.level must be one of 0.90 or 0.95.")
  }
  x <- if (is.null(x)) matrix(1, n, 1L) else .yamazaki.matrix(x, n, "x")
  z <- if (is.null(z)) matrix(numeric(), n, 0L) else .yamazaki.matrix(z, n, "z")
  k <- ncol(x)
  if (k < 1L || k > 9L) stop("The critical-value tables support 1 to 9 breaking coefficients.")

  # Keep actual candidate fractions inside the calibration domain. Rounding up
  # the lower bound differs from floor(.05*n) only when n is not a multiple of 20.
  lower <- ceiling(trim * n)
  upper <- floor((1 - trim) * n)
  if (lower < k + 2L || n - upper < k + 2L || lower >= upper ||
      n <= 2L * k + ncol(z) + 1L) {
    stop("Too few observations in the trimmed regimes for this regressor dimension.")
  }
  candidates <- seq.int(lower, upper)

  base <- cbind(x, z)
  no.break <- .yamazaki.fit(y, base, "no-break regression")
  estimated.direction <- is.null(direction)
  if (estimated.direction) {
    f <- .yamazaki.direction(x, no.break$residuals)
    f <- cbind(f, f)
  } else {
    if (!is.numeric(direction) || is.complex(direction) || any(!is.finite(direction))) {
      stop("direction must be a finite numeric vector or a k by 2 matrix.")
    }
    if (is.null(dim(direction))) {
      if (length(direction) != k) stop("direction must have length ncol(x).")
      f <- cbind(direction, direction)
    } else {
      if (!identical(dim(direction), c(k, 2L))) stop("direction must be a k by 2 matrix.")
      f <- direction
    }
    for (j in 1:2) f[, j] <- .yamazaki.unit(f[, j], "direction")
  }

  # Fit each maintained break once for the SSR estimate and subsequent scores.
  null.fits <- lapply(candidates, function(tau) {
    .yamazaki.fit(y, cbind(x, x * (seq_len(n) > tau), z),
                 paste("regression at candidate", tau))
  })
  ssr <- vapply(null.fits, function(fit) sum(fit$residuals^2), numeric(1L))
  estimate <- candidates[which.min(ssr)]
  estimated.break.x <- x * (seq_len(n) > estimate)

  methods <- c("EM", "MEM", "SUP", "AVG", "EXP")
  statistics <- original <- critical <- original.critical <- matrix(
    NA_real_, length(candidates), 5L, dimnames = list(candidates, methods)
  )
  augmented <- logical(length(candidates))
  hac.diagnostics <- c(auxiliary.rho.caps = 0L, unstable.VAR.fits = 0L)

  for (i in seq_along(candidates)) {
    tau <- candidates[i]
    before <- seq_len(tau)
    after <- seq.int(tau + 1L, n)
    score <- x * null.fits[[i]]$residuals
    modified.residuals <- null.fits[[i]]$residuals
    # The additional break is used only for covariance estimation.
    # A separating regime needs at least k rows to identify its k coefficients.
    if (abs(tau - estimate) >= k) {
      modified.residuals <- .yamazaki.fit(
        y, cbind(x, x * (seq_len(n) > tau), z, estimated.break.x),
        paste("augmented covariance regression at candidate", tau)
      )$residuals
      augmented[i] <- TRUE
    }
    modified.score <- x * modified.residuals
    em.root <- lapply(list(score[before, , drop = FALSE], score[after, , drop = FALSE]),
                      function(s) .yamazaki.invroot(.yamazaki.cov(s, covariance)))
    mem.root <- lapply(list(modified.score[before, , drop = FALSE],
                            modified.score[after, , drop = FALSE]),
                       function(s) .yamazaki.invroot(.yamazaki.cov(s, covariance)))
    if (covariance == "HAC") {
      for (root in c(em.root, mem.root)) {
        hac.diagnostics <- hac.diagnostics + attr(root, "hac.diagnostics")
      }
    }
    prefix <- list(.yamazaki.cumsum(score[before, , drop = FALSE]),
                   .yamazaki.cumsum(score[after, , drop = FALSE]))
    w.em <- w.mem <- vector("list", 2L)
    w.em.original <- w.mem.original <- vector("list", 2L)
    for (j in 1:2) {
      w.em.original[[j]] <- prefix[[j]] %*% em.root[[j]]
      w.mem.original[[j]] <- prefix[[j]] %*% mem.root[[j]]
      w.em[[j]] <- w.em.original[[j]] %*% t(.yamazaki.rotation(f[, j], em.root[[j]]))
      w.mem[[j]] <- w.mem.original[[j]] %*% t(.yamazaki.rotation(f[, j], mem.root[[j]]))
    }
    em <- sum(.yamazaki.q(w.em[[1L]])) / tau^2 +
      sum(.yamazaki.q(w.em[[2L]])) / (n - tau)^2
    mem <- sum(.yamazaki.q(w.mem[[1L]])) / tau^2 +
      sum(.yamazaki.q(w.mem[[2L]])) / (n - tau)^2
    em.original <- sum(w.em.original[[1L]]^2) / tau^2 +
      sum(w.em.original[[2L]]^2) / (n - tau)^2
    mem.original <- sum(w.mem.original[[1L]]^2) / tau^2 +
      sum(w.mem.original[[2L]]^2) / (n - tau)^2
    # Equation (14): the negative tail on the left is the positive prefix,
    # because the maintained regression makes the segment's total score zero.
    # Set T2=T1 exactly to zero to avoid dependence on floating-point OLS noise.
    path <- path.original <- numeric(length(candidates))
    left <- candidates < tau
    right <- candidates > tau
    path[left] <- .yamazaki.q(w.mem[[1L]][candidates[left], , drop = FALSE]) / n
    path[right] <- .yamazaki.q(w.mem[[2L]][candidates[right] - tau, , drop = FALSE]) / n
    path.original[left] <- rowSums(w.mem.original[[1L]][candidates[left], , drop = FALSE]^2) / n
    path.original[right] <- rowSums(w.mem.original[[2L]][candidates[right] - tau, , drop = FALSE]^2) / n
    lambda <- tau / n
    lambda2 <- candidates / n
    omega <- ifelse(left, lambda2 * (lambda - lambda2) / lambda,
                    (lambda2 - lambda) * (1 - lambda2) / (1 - lambda))
    exp.keep <- abs(candidates - tau) > floor(trim * n)
    if (!any(exp.keep)) stop("No alternative dates remain for the EXP statistic.")
    expo <- .yamazaki.logmeanexp(-k / 2 * log(omega[exp.keep]) +
                                  path[exp.keep] / (2 * omega[exp.keep]))
    expo.original <- .yamazaki.logmeanexp(-k / 2 * log(omega[exp.keep]) +
                                           path.original[exp.keep] / (2 * omega[exp.keep]))
    statistics[i, ] <- c(em, mem, max(path), mean(path), expo)
    original[i, ] <- c(em.original, mem.original, max(path.original),
                       mean(path.original), expo.original)
    critical[i, ] <- .yamazaki.cv(lambda, conf.level, k, modified = TRUE)
    original.critical[i, ] <- .yamazaki.cv(lambda, conf.level, k, modified = FALSE)
  }
  accepted <- statistics <= critical
  original.accepted <- original <= original.critical
  masks <- original.masks <- matrix(FALSE, n, 5L, dimnames = list(seq_len(n), methods))
  masks[candidates, ] <- accepted
  original.masks[candidates, ] <- original.accepted
  csets <- setNames(lapply(seq_along(methods), function(j) candidates[accepted[, j]]), methods)
  original.csets <- setNames(lapply(seq_along(methods), function(j) candidates[original.accepted[, j]]), methods)
  result <- list(call = call, nobs = n, k = k, conf.level = conf.level,
                 estimate = estimate, candidates = candidates, time = time,
                 critical.lambda = candidates / n, trim = trim,
                 direction = f, direction.estimated = estimated.direction,
                 covariance = covariance, augmented = augmented, ssr = ssr,
                 hac.diagnostics = hac.diagnostics,
                 stats = statistics, original.stats = original,
                 cv = critical, original.cv = original.critical,
                 csets = csets, original.csets = original.csets,
                 masks = masks, original.masks = original.masks)
  class(result) <- "bt_confSetYamazaki"
  result
}

#' @rdname conf.sets.Yamazaki
#' @param ... Unused arguments to the print method.
#' @export
print.bt_confSetYamazaki <- function(x, ...) {
  cat("Yamazaki (2019) confidence sets for one stationary-regression break\n")
  cat("Confidence level:", x$conf.level, " | SSR break estimate:", x$estimate, "\n")
  cat("Direction:", if (x$direction.estimated) "estimated (Procedure B)" else "supplied",
      " | Score covariance:", x$covariance, "\n")
  print(data.frame(method = colnames(x$stats),
                   directional.dates = lengths(x$csets),
                   two.sided.dates = lengths(x$original.csets)), row.names = FALSE)
  cat("Each set is in $csets; indices denote the last observation before the break.\n")
  invisible(x)
}

.yamazaki.matrix <- function(a, n, label) {
  if (!is.numeric(a) || is.complex(a) || any(!is.finite(a)) ||
      (!is.null(dim(a)) && length(dim(a)) != 2L)) {
    stop(label, " must be a finite numeric vector or matrix.")
  }
  a <- as.matrix(a)
  if (nrow(a) != n || ncol(a) < 1L) stop(label, " must have length(y) rows and at least one column.")
  a
}

.yamazaki.fit <- function(y, design, context) {
  q <- qr(design, tol = 1e-10)
  if (q$rank < ncol(design)) stop("Rank-deficient ", context, ".")
  if (nrow(design) <= q$rank) stop("No residual degrees of freedom in ", context, ".")
  residuals <- as.numeric(qr.resid(q, y))
  # Reject exact fits using a tolerance relative to the response scale.
  scale <- max(abs(y))
  if (scale == 0 || sqrt(mean((residuals / scale)^2)) <=
      100 * .Machine$double.eps * sqrt(length(y))) {
    stop("Zero or numerically zero residual variance in ", context, ".")
  }
  list(residuals = residuals)
}

.yamazaki.cumsum <- function(a) {
  out <- a
  for (j in seq_len(ncol(a))) out[, j] <- cumsum(a[, j])
  out
}

.yamazaki.q <- function(a) rowSums(abs(a) * a)

.yamazaki.unit <- function(v, label) {
  scale <- max(abs(v))
  if (!is.finite(scale) || scale == 0) stop("Cannot determine a nonzero ", label, ".")
  v <- v / scale
  v / sqrt(sum(v^2))
}

.yamazaki.direction <- function(x, residuals) {
  s <- .yamazaki.cumsum(x * residuals)
  scale <- max(abs(s))
  if (!is.finite(scale) || scale == 0) stop("Cannot estimate break direction from zero scores.")
  s <- s / scale
  b <- -colSums(abs(s) * s)
  .yamazaki.unit(sign(b) * sqrt(abs(b)), "estimated break direction")
}

.yamazaki.invroot <- function(omega) {
  if (any(!is.finite(omega))) stop("Non-finite score covariance estimate.")
  e <- eigen((omega + t(omega)) / 2, symmetric = TRUE)
  if (max(e$values) <= 0 || min(e$values) <= 1e-10 * max(e$values)) {
    stop("Score covariance is singular or not positive definite; inspect the data/regressors.")
  }
  root <- e$vectors %*% (t(e$vectors) / sqrt(e$values))
  attr(root, "hac.diagnostics") <- attr(omega, "hac.diagnostics")
  root
}

.yamazaki.basis <- function(v) {
  k <- length(v)
  basis <- matrix(.yamazaki.unit(v, "rotation direction"), k, 1L)
  if (k == 1L) return(basis)
  # This is the paper's Gram--Schmidt order wherever it is nondegenerate.
  # Trying e1 last also covers a zero first component, omitted by Procedure A.
  axes <- diag(k)[, c(seq.int(2L, k), 1L), drop = FALSE]
  for (j in seq_len(k)) {
    a <- axes[, j]
    a <- a - as.vector(basis %*% crossprod(basis, a))
    a <- a - as.vector(basis %*% crossprod(basis, a))
    norm <- sqrt(sum(a^2))
    if (norm > 1e-10) basis <- cbind(basis, a / norm)
    if (ncol(basis) == k) break
  }
  if (ncol(basis) != k) stop("Unable to construct an orthogonal direction rotation.")
  basis
}

.yamazaki.rotation <- function(f, invroot) {
  v <- as.vector(invroot %*% f)
  k <- length(v)
  -.yamazaki.basis(rep(1 / sqrt(k), k)) %*% t(.yamazaki.basis(v))
}

.yamazaki.logmeanexp <- function(v) {
  peak <- max(v)
  if (!is.finite(peak)) stop("Non-finite EXP log weights; inspect scale or covariance.")
  peak + log(mean(exp(v - peak)))
}

.yamazaki.cov <- function(score, type) {
  n <- nrow(score)
  k <- ncol(score)
  if (type == "iid") return(crossprod(score) / n)
  # VAR(1) prewhitening without intercept, appropriate for the zero-mean scores.
  prev <- score[-n, , drop = FALSE]
  next.score <- score[-1L, , drop = FALSE]
  q <- qr(prev, tol = 1e-10)
  if (q$rank < k) stop("Rank-deficient VAR(1) in HAC prewhitening.")
  coefficient <- qr.coef(q, next.score)
  unstable.var <- as.integer(max(Mod(eigen(coefficient, only.values = TRUE)$values)) >= 1)
  innovation <- next.score - prev %*% coefficient
  recolour <- diag(k) - t(coefficient)
  if (rcond(recolour) <= 1e-10) stop("Singular HAC recolouring matrix.")
  inverse <- solve(recolour)
  m <- nrow(innovation)
  alpha.numerator <- alpha.denominator <- 0
  auxiliary.caps <- 0L
  for (j in seq_len(k)) {
    v <- innovation[, j]
    vprev <- v[-m]
    vnext <- v[-1L]
    denom <- sum(vprev^2)
    if (!is.finite(denom) || denom <= 0) stop("Degenerate HAC bandwidth regression.")
    rho <- sum(vprev * vnext) / denom
    if (!is.finite(rho)) stop("Non-finite auxiliary AR(1) in HAC bandwidth selection.")
    # Small edge regimes can yield |rho| >= 1. Record the use of the same
    # 0.97 bound used elsewhere in breaktest.
    auxiliary.caps <- auxiliary.caps + as.integer(abs(rho) > 0.97)
    rho <- max(-0.97, min(0.97, rho))
    sigma2 <- mean((vnext - rho * vprev)^2)
    alpha.numerator <- alpha.numerator + 4 * rho^2 * sigma2^2 / (1 - rho)^8
    alpha.denominator <- alpha.denominator + sigma2^2 / (1 - rho)^4
  }
  if (alpha.denominator <= 0) stop("Degenerate HAC bandwidth estimate.")
  bandwidth <- 1.3221 * (m * alpha.numerator / alpha.denominator)^0.2
  omega <- crossprod(innovation) / m
  if (bandwidth > sqrt(.Machine$double.eps)) {
    for (lag in seq_len(m - 1L)) {
      s <- 6 * pi * lag / (5 * bandwidth)
      weight <- 3 / s^2 * (sin(s) / s - cos(s))
      gamma <- crossprod(innovation[seq_len(m - lag), , drop = FALSE],
                         innovation[seq.int(lag + 1L, m), , drop = FALSE]) / m
      omega <- omega + weight * (gamma + t(gamma))
    }
  }
  result <- inverse %*% omega %*% t(inverse)
  attr(result, "hac.diagnostics") <- c(auxiliary.rho.caps = auxiliary.caps,
                                      unstable.VAR.fits = unstable.var)
  result
}
