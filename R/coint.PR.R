#' @title
#' A set of residual based tests for cointegration
#'
#' @details
#' Perron and Rodríguez (2016) provide generalized least-squares (GLS) detrended versions
#' of single-equation static regression or residuals-based tests
#' for testing whether or not non-stationary time series are cointegrated.
#' Their approach is to consider nearly optimal tests for unit roots
#' and to apply them in the cointegration context.
#' Authors derive the local asymptotic power functions of all tests
#' considered for a triangular data-generating process,
#' imposing a directional restriction such that the regressors are pure integrated processes.
#' Their GLS versions of the tests do indeed provide substantial power improvements
#' over their ordinary least-squares counterparts.
#'
#' The code provided is the original GAUSS code by Perron and Rodríguez
#' ported to R.
#'
#' @param y,x Variables of interest. `x` can be a matrix of several variables.
#' @param deter A value equal to
#' * 1: quasi-demeaned y and x,
#' * 2: quasi-detrended y and x,
#' * 3: quasi-demeaned y and quasi-detrended x.
#'
#' @return A list of:
#' * 7x1-matrix of test statistics values,
#' * estimated number of lags.
#'
#' @references
#' Perron, Pierre, and Gabriel Rodríguez.
#' “Residuals‐based Tests for Cointegration with Generalized Least‐squares
#' Detrended Data.”
#' The Econometrics Journal 19, no. 1 (February 1, 2016): 84–111.
#' https://doi.org/10.1111/ectj.12056.
#'
#' @export
coint.PR <- function(y, x, deter, kmax = round(4 * (N / 100)^(1 / 4)), signif = 0.05) {
  if (!signif %in% c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2)) {
    stop(
      "ERROR! `signif` should be one of (0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2)!"
    )
  }
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  N <- nrow(y)
  Nc <- min(ncol(x), 5)
  i_sign <- which(c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2) == signif)

  opt_cbar <- .cbar_PR[[Nc]][deter]



  zy <- cbind(
    .const(N),
    if (deter == 2) .trend(N) else NULL
  )

  zx <- cbind(
    .const(N),
    if (deter %in% c(2, 3)) .trend(N) else NULL
  )

  y.d <- GLS.reg(y, zy, opt_cbar)$residuals
  x.d <- GLS.reg(x, zx, opt_cbar)$residuals

  model <- OLS.reg(y.d, x.d)
  uhat <- cbind(model$residuals)

  c_values <- .cval_PR[[deter]]
  result <- statistics.PR(uhat, kmax, opt_cbar, deter)

  result[["MZ(rho)"]]$c.value <- c_values[["Z"]][i_sign, Nc]
  result[["MSB"]]$c.value <- c_values[["MSB"]][i_sign, Nc]
  result[["MZ(t.rho)"]]$c.value <- c_values[["Zt"]][i_sign, Nc]
  result[["ADF"]]$c.value <- c_values[["Zt"]][i_sign, Nc]
  result[["Z(rho)"]]$c.value <- c_values[["Z"]][i_sign, Nc]
  result[["Z(t.rho)"]]$c.value <- c_values[["Zt"]][i_sign, Nc]
  if (deter == 1 || deter == 3) {
    result[["MP(T, demeaned)"]]$c.value <- c_values[["MPt"]][i_sign, Nc]
  } else if (deter == 2) {
    result[["MP(T, detrended)"]]$c.value <- c_values[["MPt"]][i_sign, Nc]
  }

  result
}


#' @title
#' Internal procedure for calculating test statistics from
#' Perron-Rodriguez (2016)
#'
#' @param ud A vector of residuals for testing.
#' @param kmax Maximum lag number.
#' @param c.bar A `c` parameter used for GLS detrending purposes.
#' @param deter A value equal to
#' * 1: quasi-demeaned y and x,
#' * 2: quasi-detrended y and x,
#' * 3: quasi-demeaned y and quasi-detrended x.
#'
#' @return A list of:
#' * 7x1-matrix of test statistics values,
#' * estimated number of lags.
#'
#' @keywords internal
statistics.PR <- function(u, kmax, c.bar, deter) {
  N <- length(u) - 1

  lm_model <- OLS.reg(u, .lagn(u, 1))
  rho <- lm_model$coefficients[1]
  ee_1 <- na.omit(lm_model$residuals)

  su2 <- sum(ee_1^2) / N
  t_rho <- lm_model$t.stats[1]
  sum_ud <- sum(u[1:N]^2)
  uT <- u[N + 1]

  klag <- uroot.ADF(u, FALSE, FALSE, kmax, "aic", TRUE)$lag
  model_2 <- uroot.ADF(u, FALSE, FALSE, klag, NULL)

  e2 <- na.omit(model_2$model$residuals)
  sk2 <- sum(e2^2) / N
  sumb <- if (klag == 0) 0 else sum(model_2$model$coefficients[-1])
  s2 <- sk2 / ((1 - sumb)^2)

  result <- list()

  result[["MZ(rho)"]] <- list(statistic = (uT^2 / N - s2) / (2 * sum_ud / N^2))
  result[["MSB"]] <- list(statistic = sqrt(sum_ud / s2 / N^2))
  result[["MZ(t.rho)"]] <- list(
    statistic = (uT^2 / N - s2) / sqrt(4 * s2 * sum_ud / N^2)
  )
  result[["ADF"]] <- list(statistic = model_2$t.alpha)
  result[["Z(rho)"]] <- list(
    statistic = N * (rho - 1) - (s2 - su2) / 2 / sum_ud / N^2
  )
  result[["Z(t.rho)"]] <- list(
    statistic = sqrt(su2 / s2) *
      t_rho -
      (s2 - su2) / sqrt(4 * s2 * sum_ud / N^2)
  )
  if (deter == 1 || deter == 3) {
    result[["MP(T, demeaned)"]] <- list(
      statistic = (c.bar^2 * sum_ud / N^2 - c.bar * uT^2 / N) /
        s2
    )
  } else if (deter == 2) {
    result[["MP(T, detrended)"]] <- list(
      statistic = (c.bar^2 * sum_ud / N^2 + (1 - c.bar) * uT^2 / N) /
        s2
    )
  }

  result$lag <- klag
  class(result) <- "bt_cointPR"

  result
}
