bootstrap.bt_mdfCHLT <- function(obj, iter = 999, y) {
  N <- length(y)
  dy <- diff(y)

  trim <- obj$params$trim
  tb_dy <- obj$params$tb
  tau_lam_MZ <- obj$params$MZ$tau
  cbar_tau_lam_MZ <- obj$params$MZ$cbar
  tau_lam_ADF <- obj$params$ADF$tau
  cbar_tau_lam_ADF <- obj$params$ADF$cbar

  r <- c(
    0,
    OLS.reg(
      dy,
      cbind(.const(N), .du(tb_dy, N))[-1, ]
    )$residuals
  )
  result <- NULL
  for (i in 1:iter) {
    z <- rnorm(N)
    y_wb <- cumsum(r * z)
    if (tau_lam_MZ < trim) {
      r_GLS_t_wb <- GLS.reg(
        y_wb,
        cbind(.const(N), .trend(N)),
        -13.5
      )$residuals

      MZ_wb <- .mz.statistics(r_GLS_t_wb, 0)
      MZa_wb <- MZ_wb$mza
      MSB_wb <- MZ_wb$msb
      MZt_wb <- MZ_wb$mzt

      rm(MZ_wb)
    } else {
      resid.wb <- GLS.bt(y, tau_lam_MZ, cbar_tau_lam_MZ)$residuals

      MZ_wb <- .mz.statistics(resid.wb, 0)
      MZa_wb <- MZ_wb$mza
      MSB_wb <- MZ_wb$msb
      MZt_wb <- MZ_wb$mzt

      rm(MZ_wb)
    }
    if (tau_lam_ADF < trim) {
      r_GLS_t_wb <- GLS.reg(
        y_wb,
        cbind(.const(N), .trend(N)),
        -13.5
      )$residuals

      ers_ADF_wb <- uroot.ADF(
        r_GLS_t_wb,
        const = FALSE,
        trend = FALSE,
        max.lag = 0,
        criterion = NULL
      )$t.alpha
    } else {
      r_GLS_wb <- GLS.bt(y, tau_lam_ADF, cbar_tau_lam_ADF)$residuals

      ers_ADF_wb <- uroot.ADF(
        r_GLS_wb,
        const = FALSE,
        trend = FALSE,
        max.lag = 0,
        criterion = NULL
      )$t.alpha
    }
    result <- rbind(result, c(MZa_wb, MSB_wb, MZt_wb, ers_ADF_wb))
  }
  list(
    MZa = sort(result[, 1])[trunc(0.05 * iter)],
    MSB = sort(result[, 2])[trunc(0.05 * iter)],
    MZt = sort(result[, 3])[trunc(0.05 * iter)],
    ADF = sort(result[, 4])[trunc(0.05 * iter)]
  )
}
