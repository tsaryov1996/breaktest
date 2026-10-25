#' Functions to receive bootstrapped p-values.
#'
#' @description
#' `bootstrap` is a generic function aimed to handle objects of classes provided in this package.
#'
#' @param obj An object of one of the following classes:
#' * `bt_adf`,
#' * `bt_kpss`,
#' * `bt_SADF` and `bt_GSADF`,
#' * `bt_mdfCHLT`.
#'
#' @details
#' `bootstrap` is a generic function aimed to handle objects of classes provided in this package.
#' Normally one should use special function parameters to obtain bootstrapped critical values or
#' p-values (depends on test).
#' If later you decide to calculate one then you call `bootstrap` function on the object
#' with corresponding results.
#'
#' @return
#' Bootstrapped critical value or p-value.
#'
#' @export
bootstrap <- function(obj, ...) UseMethod("bootstrap")


#' @rdname bootstrap
#' @param iter Number of bootstrapping iterations.
#' @exportS3Method
bootstrap.bt_adf <- function(obj, iter = 999) {
  vCoefs <- obj$model$coefficients[-1]
  vEps <- obj$model$residuals

  cLag <- obj$lag

  cN <- length(vEps)
  mDeter <- cbind(
    if (obj$const) .const(cN) else NULL,
    if (obj$trend) .trend(cN) else NULL
  )

  result <- NULL
  for (i in 1:iter) {
    u <- rep(0, cLag + cN)
    eps <- sample(vEps, cN, replace = TRUE)
    if (cLag > 0) {
      for (s in 1:cN) {
        u[cLag + s] <- u[(cLag + s - 1):s] %*% vCoefs + eps[s]
      }
      u <- u[-(1:cLag)]
    } else {
      for (s in 1:cN) {
        u[s] <- eps[s]
      }
    }
    tmp.y <- as.matrix(cumsum(u))
    if (obj$recursive) {
      tmp.y <- detrend.recursively(
        tmp.y,
        mDeter,
        obj$recursive.params$cc,
        obj$recursive.params$gamma,
        obj$recursive.params$trim
      )
    }
    tmp.res <- OLS.reg(
      tmp.y,
      obj$model$exog
    )
    result <- c(result, tmp.res$t.stats[1])
  }

  sum(result < obj$t.alpha) / iter
}


#' @rdname bootstrap
#' @param boot.type Type of bootstrapping:
#' \itemize{
#' \item{"sample"}{sampling from residuals with replacement,}
#' \item{"Cavaliere-Taylor"}{multiplying residuals by \eqn{N(0, 1)}-distributed}
#' variable,
#' \item{"Rademacher"}{multiplying residuals by Rademacher-distributed variable.}
#' }
#' @importFrom stats rnorm
#' @exportS3Method
bootstrap.bt_kpss <- function(obj, iter = 999, boot.type = "sample") {
  xreg <- obj$exog
  u <- obj$residuals
  max.lag <- obj$lr.var.max.lag
  kernel <- obj$lr.var.kernel

  result <- NULL
  for (i in 1:iter) {
    y.loop <- switch(
      boot.type,
      "sample" = sample(u, length(u), replace = TRUE),
      "Cavaliere-Taylor" = rnorm(length(u)) * u,
      "Rademacher" = sample(c(-1, 1), length(u), replace = TRUE) * u,
      stop(
        "ERROR! bootstsrap.bt_kpss: Unknown bootstrap type '",
        boot.type,
        "'"
      )
    )
    resids <- OLS.reg(y.loop, xreg)$residuals

    result <- c(
      result,
      if (kernel == "Kurozumi") {
        .kpss.statistic(resids, .lr.var.kurozumi(resids))
      } else {
        .kpss.statistic(resids, .lr.var.spc(resids, max.lag, kernel))
      }
    )
  }

  (1 / iter) * sum(obj$statistic <= result)
}


#' @rdname bootstrap
#' @exportS3Method
bootstrap.bt_SADF <- function(obj, iter = 999) {
  y <- obj$y
  N <- length(y)

  trim <- obj$trim
  const <- obj$const

  model <- uroot.SADF(y, trim, const)
  SADF.value <- model$SADF.value

  result <- NULL
  for (i in 1:iter) {
    y.star <- cumsum(rnorm(N - 1) * .diffn(y, na = 0))
    result <- c(result, uroot.SADF(y.star, trim, const)$SADF.value)
  }

  sum(result > SADF.value) / iter
}


#' @rdname bootstrap
#' @exportS3Method
bootstrap.bt_GSADF <- function(obj, iter = 999) {
  y <- obj$y

  trim <- obj$trim

  const <- obj$const

  N <- length(y)
  ## Find GSADF.value.

  model <- uroot.GSADF(y, trim, const)
  GSADF.value <- model$GSADF.value

  result <- NULL
  for (i in 1:iter) {
    y.star <- cumsum(rnorm(N - 1) * .diffn(y, na = 0))
    result <- c(result, uroot.GSADF(y.star, trim, const)$GSADF.value)
  }
  sum(result > GSADF.value) / iter
}


#' @rdname bootstrap
#' @param y Original series, required for `bt_mdfCHLT` objects.
#' @details For `uroot.CHLT(..., breaks = 2)` or `3`, uses the second-difference
#' wild bootstrap of Algorithm 1 in Cavaliere et al. (January 2014), repeating
#' the complete date search with zero ADF lags. Returns `ADF` (the empirical
#' 5% critical value), `p.value`, `iterations`, `difference.order`, `lag`,
#' `innovations` and `statistics`. The one-break bootstrap is unchanged.
#' @exportS3Method
bootstrap.bt_mdfCHLT <- function(obj, iter = 999, y) {
  if (identical(obj$params$method, "CHLT2014-MDF")) {
    return(.chlt.bootstrap.multiple(obj, iter, y))
  }

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
