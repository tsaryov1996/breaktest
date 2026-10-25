#' @title
#' MDF test for a single break and possible heteroscedasticity
#'
#' @details
#' The code provided is the original GAUSS code by Cavaliere et al.
#' ported to R.
#'
#' @param y A time series of interest.
#' @param max.lag The maximum possible lag.
#' @param trim Trimming parameter for lag selection
#' @param boot.cv Whether we need bootstrapped critical values.
#' @param boot.iter Number of bootstrap iterations.
#'
#' @return An object of type `mdfCHLT`. It's a list of four sublists
#' each containing:
#' * The value of \eqn{MZ_\alpha}, \eqn{MSB}, \eqn{MZ_t}, or \eqn{ADF},
#' * The asymptotic c.v.,
#' * The bootstrapped c.v.
#'
#' @references
#' Cavaliere, Giuseppe, David I. Harvey, Stephen J. Leybourne,
#' and A.M. Robert Taylor.
#' “Testing for Unit Roots in the Presence of a Possible Break in Trend and
#' Nonstationary Volatility.”
#' Econometric Theory 27, no. 5 (October 2011): 957–91.
#' https://doi.org/10.1017/S0266466610000605.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
uroot.CHLT <- function(
  y,
  max.lag = NULL,
  trim = 0.15,
  boot.cv = FALSE,
  boot.iter = 999
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  ## CV ##
  g_brk_MZ <- .cval_MDF_CHLT$g.brk.MZ
  g_brk_ADF <- .cval_MDF_CHLT$g.brk.ADF
  cv_MZa_lim <- .cval_MDF_CHLT$cv.MZa.lim
  cv_MSB_lim <- .cval_MDF_CHLT$cv.MSB.lim
  cv_MZt_lim <- .cval_MDF_CHLT$cv.MZt.lim
  cv_ers_lim <- .cval_MDF_CHLT$cv.ers.lim
  # MZa_cv_vals_lim <- .cval_MDF_CHLT$MZa.cv.vals.lim
  # MSB_cv_vals_lim <- .cval_MDF_CHLT$MSB.cv.vals.lim
  # MZt_cv_vals_lim <- .cval_MDF_CHLT$MZt.cv.vals.lim
  # ADF_cv_vals_lim <- .cval_MDF_CHLT$ADF.cv.vals.lim
  # cbar_vals <- .cval_MDF_CHLT$cbar.vals
  tau_cbar_MZa_cv_lim <- .cval_MDF_CHLT$tau.cbar.MZa.cv.lim
  tau_cbar_MSB_cv_lim <- .cval_MDF_CHLT$tau.cbar.MSB.cv.lim
  tau_cbar_MZt_cv_lim <- .cval_MDF_CHLT$tau.cbar.MZt.cv.lim
  tau_cbar_ADF_cv_lim <- .cval_MDF_CHLT$tau.cbar.ADF.cv.lim

  ## Start
  N <- nrow(y)
  tr <- .trend(N)

  if (is.null(max.lag)) {
    max.lag <- trunc(12 * ((N / 100)^(1 / 4)))
  }

  tb_L <- trunc(trim * N)
  tb_U <- trunc((1 - trim) * N)

  ## Break date estimation
  ssr_dy <- Inf
  tb_dy <- 0

  for (tb in tb_L:tb_U) {
    b1 <- (y[tb] - y[1]) / (tb - 1)
    b2 <- (y[N] - y[tb]) / (N - tb) - b1
    loop_ssr <- (N - 1) *
      b1^2 +
      (N - tb) * b2^2 -
      2 * b1 * (y[N] - y[1]) -
      2 * b2 * (y[N] - y[tb]) +
      2 * b1 * b2 * (N - tb)
    if (loop_ssr < ssr_dy) {
      ssr_dy <- loop_ssr
      tb_dy <- tb
    }
  }

  tau_dy <- tb_dy / N
  # du_tb_dy <- .dt(tb_dy, N)
  dt_tb_dy <- .dt(tb_dy, N)

  z <- cumsum(y)
  x <- cbind(tr, cumsum(tr), cumsum(dt_tb_dy))
  ru <- OLS.reg(z, x)$residuals
  x <- cbind(tr, cumsum(tr))
  rr <- OLS.reg(z, x)$residuals
  W_stat_dy <- sum(rr^2) / sum(ru^2) - 1

  lam_MZ_brk_tau_dy <- exp(-g_brk_MZ * W_stat_dy / sqrt(N))
  tau_lam_MZ <- (1 - lam_MZ_brk_tau_dy) * tau_dy

  lam_ADF_brk_tau_dy <- exp(-g_brk_ADF * W_stat_dy / sqrt(N))
  tau_lam_ADF <- (1 - lam_ADF_brk_tau_dy) * tau_dy

  ## Unit root test
  x_OLS_t <- cbind(.const(N), .trend(N))

  r_GLS_t <- GLS.reg(y, x_OLS_t, -13.5)$residuals
  r_OLS_t <- OLS.reg(y, x_OLS_t)$residuals
  k_t <- uroot.ADF(
    r_OLS_t,
    const = FALSE,
    trend = FALSE,
    max.lag = max.lag,
    criterion = "aic",
    modified.criterion = TRUE
  )$lag
  ers_t_k <- uroot.ADF(
    r_GLS_t,
    const = FALSE,
    trend = FALSE,
    max.lag = k_t,
    criterion = NULL
  )$t.alpha

  tmp.MZ <- .mz.statistics(r_GLS_t, k_t)
  MZa <- tmp.MZ$mza
  MSB <- tmp.MZ$msb
  MZt <- tmp.MZ$mzt
  rm(tmp.MZ)

  ## MZ
  if (tau_lam_MZ < trim) {
    HHLT_MZa_kMAIC <- MZa
    HHLT_MZa_kMAIC_cv <- cv_MZa_lim
    HHLT_MSB_kMAIC <- MSB
    HHLT_MSB_kMAIC_cv <- cv_MSB_lim
    HHLT_MZt_kMAIC <- MZt
    HHLT_MZt_kMAIC_cv <- cv_MZt_lim
  } else {
    if (tau_lam_MZ == trim) {
      cbar_tau_lam_MZ <- tau_cbar_MZa_cv_lim[1, 2]
      cv_MZa_tau_lam_MZ_lim <- tau_cbar_MZa_cv_lim[1, 3]
      cv_MSB_tau_lam_MZ_lim <- tau_cbar_MSB_cv_lim[1, 3]
      cv_MZt_tau_lam_MZ_lim <- tau_cbar_MZt_cv_lim[1, 3]
    } else if (tau_lam_MZ == 1 - trim) {
      cbar_tau_lam_MZ <- tau_cbar_MZa_cv_lim[15, 2]
      cv_MZa_tau_lam_MZ_lim <- tau_cbar_MZa_cv_lim[15, 3]
      cv_MSB_tau_lam_MZ_lim <- tau_cbar_MSB_cv_lim[15, 3]
      cv_MZt_tau_lam_MZ_lim <- tau_cbar_MZt_cv_lim[15, 3]
    } else {
      tau_near_index <- which.min(
        abs(tau_lam_MZ - tau_cbar_MZa_cv_lim[, 1])
      )
      tau_near <- tau_cbar_MZa_cv_lim[tau_near_index, 1]

      if (tau_lam_MZ > tau_near) {
        tau_l <- tau_near
        tau_l_index <- tau_near_index
        tau_u <- tau_near + 0.05
        tau_u_index <- tau_near_index + 1
      } else {
        tau_l <- tau_near - 0.05
        tau_l_index <- tau_near_index - 1
        tau_u <- tau_near
        tau_u_index <- tau_near_index
      }

      weight_l <- 1 - (tau_lam_MZ - tau_l) / 0.05
      weight_u <- 1 - (tau_u - tau_lam_MZ) / 0.05

      cbar_tau_lam_MZ <-
        weight_l *
        tau_cbar_MZa_cv_lim[tau_l_index, 2] +
        weight_u * tau_cbar_MZa_cv_lim[tau_u_index, 2]
      cv_MZa_tau_lam_MZ_lim <-
        weight_l *
        tau_cbar_MZa_cv_lim[tau_l_index, 3] +
        weight_u * tau_cbar_MZa_cv_lim[tau_u_index, 3]
      cv_MSB_tau_lam_MZ_lim <-
        weight_l *
        tau_cbar_MSB_cv_lim[tau_l_index, 3] +
        weight_u * tau_cbar_MSB_cv_lim[tau_u_index, 3]
      cv_MZt_tau_lam_MZ_lim <-
        weight_l *
        tau_cbar_MZt_cv_lim[tau_l_index, 3] +
        weight_u * tau_cbar_MZt_cv_lim[tau_u_index, 3]
    }

    r_GLS_bt <- GLS.bt(y, tau_lam_MZ, cbar_tau_lam_MZ)$residuals

    tb_lam_MZ <- trunc(tau_lam_MZ * N)
    r_OLS_bt <- OLS.reg(
      y,
      cbind(.const(N), .trend(N), .dt(tb_lam_MZ, N))
    )$residuals

    k_bt <- uroot.ADF(
      r_OLS_bt,
      const = FALSE,
      trend = FALSE,
      max.lag = max.lag,
      criterion = NULL
    )$lag

    tmp.MZ <- .mz.statistics(r_GLS_bt, k_bt)
    MZa_tau_lam_MZ_k <- tmp.MZ$mza
    MSB_tau_lam_MZ_k <- tmp.MZ$msb
    MZt_tau_lam_MZ_k <- tmp.MZ$mzt
    rm(tmp.MZ)

    HHLT_MZa_kMAIC <- MZa_tau_lam_MZ_k
    HHLT_MZa_kMAIC_cv <- cv_MZa_tau_lam_MZ_lim
    HHLT_MSB_kMAIC <- MSB_tau_lam_MZ_k
    HHLT_MSB_kMAIC_cv <- cv_MSB_tau_lam_MZ_lim
    HHLT_MZt_kMAIC <- MZt_tau_lam_MZ_k
    HHLT_MZt_kMAIC_cv <- cv_MZt_tau_lam_MZ_lim
  }

  ## ADF
  if (tau_lam_ADF < trim) {
    HHLT_ADF_kMAIC <- ers_t_k
    HHLT_ADF_kMAIC_cv <- cv_ers_lim
  } else {
    if (tau_lam_ADF == trim) {
      cbar_tau_lam_ADF <- tau_cbar_ADF_cv_lim[1, 2]
      cv_ADF_tau_lam_ADF_lim <- tau_cbar_ADF_cv_lim[1, 3]
    } else if (tau_lam_ADF == 1 - trim) {
      cbar_tau_lam_ADF <- tau_cbar_ADF_cv_lim[15, 2]
      cv_ADF_tau_lam_ADF_lim <- tau_cbar_ADF_cv_lim[15, 3]
    } else {
      tau_near_index <- which.min(
        abs(tau_lam_ADF - tau_cbar_ADF_cv_lim[, 1])
      )
      tau_near <- tau_cbar_ADF_cv_lim[tau_near_index, 1]

      if (tau_lam_ADF > tau_near) {
        tau_l <- tau_near
        tau_l_index <- tau_near_index
        tau_u <- tau_near + 0.05
        tau_u_index <- tau_near_index + 1
      } else {
        tau_l <- tau_near - 0.05
        tau_l_index <- tau_near_index - 1
        tau_u <- tau_near
        tau_u_index <- tau_near_index
      }

      weight_l <- 1 - (tau_lam_ADF - tau_l) / 0.05
      weight.u <- 1 - (tau_u - tau_lam_ADF) / 0.05

      cbar_tau_lam_ADF <-
        weight_l *
        tau_cbar_ADF_cv_lim[tau_l_index, 2] +
        weight.u * tau_cbar_ADF_cv_lim[tau_u_index, 2]
      cv_ADF_tau_lam_ADF_lim <-
        weight_l *
        tau_cbar_ADF_cv_lim[tau_l_index, 3] +
        weight.u * tau_cbar_ADF_cv_lim[tau_u_index, 3]
    }

    r_GLS_bt <- GLS.bt(y, tau_lam_ADF, cbar_tau_lam_ADF)$residuals

    tb_lam_ADF <- trunc(tau_lam_ADF * N)
    r_OLS_bt <- OLS.reg(
      y,
      cbind(.const(N), .trend(N), .dt(tb_lam_ADF, N))
    )$residuals

    k_bt <- uroot.ADF(
      r_OLS_bt,
      const = FALSE,
      trend = FALSE,
      max.lag = max.lag,
      criterion = NULL
    )$lag

    # ers_tau_lam_ADF_0 <- ADF.test(
    #   r_GLS_bt,
    #   const = FALSE,
    #   trend = FALSE,
    #   max.lag = 0,
    #   criterion = NULL
    # )$t.alpha

    ers_tau_lam_ADF_k <- uroot.ADF(
      r_GLS_bt,
      const = FALSE,
      trend = FALSE,
      max.lag = k_bt,
      criterion = NULL
    )$t.alpha

    HHLT_ADF_kMAIC <- ers_tau_lam_ADF_k
    HHLT_ADF_kMAIC_cv <- cv_ADF_tau_lam_ADF_lim
  }

  result <- list(
    MZa = list(
      statistic = HHLT_MZa_kMAIC,
      cv = HHLT_MZa_kMAIC_cv
    ),
    MSB = list(
      statistic = HHLT_MSB_kMAIC,
      cv = HHLT_MSB_kMAIC_cv
    ),
    MZt = list(
      statistic = HHLT_MZt_kMAIC,
      cv = HHLT_MZt_kMAIC_cv
    ),
    ADF = list(
      statistic = HHLT_ADF_kMAIC,
      cv = HHLT_ADF_kMAIC_cv
    ),
    params = list(
      MZ = list(
        tau = tau_lam_MZ,
        cbar = if (tau_lam_MZ >= trim) cbar_tau_lam_MZ else NULL
      ),
      ADF = list(
        tau = tau_lam_ADF,
        cbar = if (tau_lam_MZ >= trim) cbar_tau_lam_ADF else NULL
      ),
      trim = trim,
      tb = tb_dy
    )
  )

  class(result) <- "bt_mdfCHLT"

  if (boot.cv) {
    boot_cv <- bootstrap(result, boot.iter, y)
    result$MZa$boot.cv <- boot_cv$MZa
    result$MSB$boot.cv <- boot_cv$MSB
    result$MZt$boot.cv <- boot_cv$MZt
    result$ADF$boot.cv <- boot_cv$ADF
  }

  result
}

#' @title
#' GLS fitering with a break
#'
#' @param y A time series of interest.
#' @param lambda A break relative position.
#' @param c A coefficient for \eqn{\rho} calculation.
#'
#' @keywords internal
GLS.bt <- function(y, lambda, c) {
  N <- nrow(y)
  tb <- trunc(lambda * N)
  x <- cbind(.const(N), .dt(tb, N))
  GLS.reg(y, x, c)
}
