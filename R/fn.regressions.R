#' @title
#' Custom OLS with extra information
#'
#' @description
#' Getting OLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#'
#' @return A list of:
#' * beta: estimates of coefficients,
#' * resid: estimated residuals,
#' * predict: forecasted values,
#' * t.beta: \eqn{t}-statistics for `beta`.
#'
#' @importFrom Rfast rowAll
#' @importFrom Rfast lmfit
#' @importFrom Rfast spdinv
#'
#' @keywords internal
OLS.reg <- function(y, x) {
  y <- as.matrix(y)
  x <- as.matrix(x)

  N <- length(y)

  rows <- rowAll(!is.na(y)) & rowAll(!is.na(x))

  y <- .msub(y, rows)
  x <- .msub(x, rows)

  cf <- solve(crossprod(x)) %*% crossprod(x, y)
  r <- y - x %*% cf

  cf <- drop(cf)
  s.sq <- sum(r^2) / (nrow(x) - ncol(x))
  se.cf <- sqrt(diag(s.sq * spdinv(t(x) %*% x)))
  t.beta <- cf / se.cf

  resid <- rep(NA, N)
  resid[rows] <- r

  fitted <- rep(NA, N)
  fitted[rows] <- y - r

  result <- list(
    coefficients = cf,
    s.sq = s.sq,
    se.coefs = se.cf,
    t.stats = t.beta,
    residuals = resid,
    fitted.values = fitted,
    endog = y,
    exog = x
  )

  class(result) <- "bt_ols"
  result
}


#' @title
#' Custom GLS with extra information
#'
#' @description
#' Getting GLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param c A coefficient for \eqn{\rho} calculation.
#'
#' @return A list of:
#' * beta: estimates of coefficients,
#' * resid: estimated residuals,
#' * predict: forecasted values,
#' * t.beta: \eqn{t}-statistics for `beta`.
#'
#' @keywords internal
GLS.reg <- function(y, z, c) {
  y <- as.matrix(y)
  z <- as.matrix(z)

  N <- nrow(y)
  Nc <- ncol(y)

  rho <- 1 + c / N

  y_hat <- y - rho * .lagn(y, 1)
  y_hat[1, ] <- y[1, ]

  z_hat <- z - rho * .lagn(z, 1)
  z_hat[1, ] <- z[1, ]

  betas <- NULL
  resids <- NULL
  fitted <- NULL
  t.betas <- NULL

  for (i in 1:Nc) {
    .model <- OLS.reg(y_hat[, i, drop = FALSE], z_hat)
    betas <- cbind(betas, .model$coefficients)
    t.betas <- cbind(t.betas, .model$t.stats)
    fitted <- cbind(fitted, z %*% .model$coefficients)
    resids <- cbind(resids, y[, i, drop = FALSE] - fitted)
  }

  result <- list(
    coefficients = betas,
    t.stats = drop(t.betas),
    residuals = resids,
    fitted.values = fitted
  )
  class(result) <- c("bt_gls", "bt_ols")
  result
}


#' @title
#' Custom AR with extra information
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param max.lag The maximum number of lags.
#' @param criterion A criterion for lag number estimation.
#'
#' @return A list of:
#' * beta: estimates of coefficients,
#' * residuals: estimated residuals,
#' * predict: forecasted values,
#' * t.beta: \eqn{t}-statistics for `beta`,
#' * lag: estimated number of lags.
#'
#' @importFrom Rfast rowAll
#' @keywords internal
AR.reg <- function(
  y,
  x,
  max.lag,
  criterion = "aic"
) {
  if (!is.null(criterion)) {
    if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
      stop("ERROR! .AR: Unknown criterion")
    }
  }

  y <- as.matrix(y)
  x <- if (!is.null(x)) as.matrix(x)

  Nx <- if (!is.null(x)) ncol(x) else 0

  for (l in seq_len(max.lag)) {
    x <- cbind(x, .lagn(y, l))
  }

  rows <- rowAll(!is.na(y)) & rowAll(!is.na(x))

  if (is.null(criterion)) {
    resLag <- max.lag
  } else {
    resLag <- 0
    minIC <- Inf

    for (l in 0:max.lag) {
      loopIC <- info.criterions(
        OLS.reg(y[rows], x[rows, 1:(Nx + l)]),
        criterion
      )

      if (loopIC < minIC) {
        minIC <- loopIC
        resLag <- l
      }
    }
  }

  rows <- rowAll(!is.na(y)) &
    rowAll(!is.na(x[, 1:(Nx + resLag), drop = FALSE]))
  result <- OLS.reg(y[rows], x[rows, 1:(Nx + resLag)])

  result$lag <- resLag
  result$criterion <- minIC
  result$criterion.name <- criterion

  class(result) <- c("bt_ar", "bt_ols")
  result
}


#' @title
#' Estimating DOLS regression for multiple known break points
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param z A matrix of extra explanatory stochastic regressors.
#' These variables are used as is without taking lagged or lead values.
#' @param model A scalar or vector of break types:
#' * 1: for the break in const.
#' * 2: for the break in trend.
#' * 3: for the break in const and trend.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#' @param criterion A criterion for lag/lead number selection.
#'
#' @return A list of:
#' * Estimates of coefficients,
#' * Estimates of residuals,
#' * A set of informational criterions values,
#' * \eqn{t}-statistics for the estimates of coefficients.
#'
#' @importFrom Rfast rowAll
#' @keywords internal
DOLS.reg <- function(
  y,
  x,
  z,
  n.lags,
  n.leads,
  criterion = "aic"
) {
  y <- as.matrix(y)

  if (is.null(x)) {
    stop("ERROR! DOLS.multiple: explanatory variables needed for DOLS")
  }
  x <- as.matrix(x)

  z <- if (!is.null(z)) z <- as.matrix(z)

  dX <- .diffn(x)

  xL <- NULL
  for (l in seq_len(n.lags)) {
    xL <- cbind(xL, .lagn(dX, l))
  }

  xF <- NULL
  for (l in seq_len(n.leads)) {
    xF <- cbind(xF, .lagn(dX, -l))
  }

  mX <- cbind(z, x, xL, xF)

  rows <- rowAll(!is.na(y)) & rowAll(!is.na(mX))

  if (is.null(criterion)) {
    resLag <- n.lags
    resLead <- n.leads
  } else {
    minIC <- Inf
    resLag <- 0
    resLead <- 0

    for (l in c(0, seq_len(n.lags))) {
      for (f in c(0, seq_len(n.leads))) {
        mX <- cbind(z, x, xL[, seq_len(l)], xF[, seq_len(f)])
        loopIC <- info.criterions(
          OLS.reg(y[rows], mX[rows, ]),
          criterion
        )

        if (loopIC < minIC) {
          minIC <- loopIC
          resLag <- l
          resLead <- f
        }
      }
    }
  }

  mX <- cbind(z, x, xL[, seq_len(resLag)], xF[, seq_len(resLead)])
  rows <- rowAll(!is.na(y)) & rowAll(!is.na(mX))

  result <- OLS.reg(y[rows], mX[rows, ])
  result$lags <- resLag
  result$leads <- resLead

  class(result) <- c("bt_dols", "bt_ols")
  result
}


#' @title
#' Nadaraya–Watson kernel regression.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * unif: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * gauss: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @return A list of arguments as well as the estimated coefficient vector and
#' residuals.
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @keywords internal
NW.reg <- function(
  y,
  x,
  h,
  kernel = "unif"
) {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.reg: unknown kernel")
  }

  N <- length(y)

  rho <- numeric(N)

  for (k in 1:N) {
    .w <- NW.kernel(k, (1:N) / N, h, kernel)
    rho[k] <- sum(x * .w * y) / sum(x * .w * x)
  }

  list(
    my = y,
    mx = x,
    h = h,
    kernel = kernel,
    rr1.est = rho,
    u.hat = y - rho * x
  )
}


#' @title
#' Nadaraya–Watson kernel volatility estimation
#'
#' @param e A series of interest.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * unif: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * gauss: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @return A list of arguments as well as the estimated omega and s.e.
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor.
#' “Lag Length Selection for Unit Root Tests in the Presence
#' of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @keywords internal
NW.variance <- function(
  e,
  h,
  kernel = "unif"
) {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.variance: unknown kernel")
  }

  .e <- na.omit(e)
  N <- length(.e)
  omega2 <- numeric(N)

  for (k in 1:N) {
    .w <- NW.kernel(k, (1:N) / N, h, kernel)
    omega2[k] <- sum(.w * e^2) / sum(.w)
  }

  dN <- length(e) - length(.e)
  omega2 <- c(rep(NA, dN), omega2)

  list(
    me = e,
    h = h,
    kernel = kernel,
    omega.sq = omega2,
    se = sqrt(omega2)
  )
}


#' @title
#' LOO-CV for h in Nadaraya–Watson kernel regression.
#'
#' @param y A dependent variable.
#' @param x An explanatory variable.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * unif: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * gauss: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @return A list of arguments as well as the estimated bandwidth `h`.
#'
#' @keywords internal
NW.bandwidth <- function(y, x, kernel = "unif") {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.bandwidth: unknown kernel")
  }

  .y <- na.omit(y)
  .x <- x[!is.na(y)]
  N <- length(.y)

  h_candidates <- seq(N^(-0.5), N^(-0.3), by = 0.01)
  rss <- Inf
  h <- NULL

  for (.h in h_candidates) {
    rho <- numeric(N)
    for (k in 1:N) {
      .w <- NW.kernel(k, (1:N) / N, .h, kernel)
      .w[k] <- 0
      rho[k] <- sum(.x * .w * .y) / sum(.x * .w * .x)
    }

    .rss <- sum((y - rho * x)^2)
    if (.rss < rss) {
      rss <- .rss
      h <- .h
    }
  }

  list(
    my = y,
    mx = x,
    kernel = kernel,
    h = h
  )
}


#' @title
#' Returning a kernel value for a specific point
#'
#' @param i An index of the base point.
#' @param x A series for kernel calculations.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * unif: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * gauss: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @importFrom stats pnorm
#'
#' @keywords internal
NW.kernel <- function(i, x, h, kernel = "unif") {
  switch(kernel,
    unif = ifelse((abs((x - x[i]) / h) <= 1), 1, 0),
    gauss = pnorm((x - x[i]) / h)
  )
}
