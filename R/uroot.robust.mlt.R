#' @title
#' A wrapping function around [breaktest.KP] and [MDF.mlt].
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y a time series of interest.
#' @param const Whether the constant term is allowed to change (Structural changes in intercept, otherwise Structural changes only in slope).
#' @param season Whether the seasonal adjustment is needed.
#' @param breaks Number of breaks.
#' @param trim Trimming value for a possible break date bounds.
#'
#' @references
#' Skrobotov, Anton.
#' “On Trend Breaks and Initial Condition in Unit Root Testing”
#' Journal of Time Series Econometrics, 2018, Issue 1.
#' https://doi.org/10.1515/jtse-2016-0014
#'
#' @export
uroot.robust.mlt <- function(
  y,
  const = TRUE,
  season = FALSE,
  breaks = 2,
  trim = 0.15
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  ## Start ##
  N <- nrow(y)

  if (season) {
    SEAS <- cbind(
      .const(N),
      seasonal.dummies(N)
    )
    y <- OLS.reg(y, SEAS)$residuals
  }

  m.star <- breaktest.KP(
    y = y,
    const = const,
    breaks = breaks,
    criterion = "bic",
    trim = trim
  )

  result <- MDF.mlt(
    y = y,
    const = const,
    breaks = breaks,
    breaks.star = m.star$breaks,
    trim = trim,
    ZA = FALSE
  )

  result$season <- season
  result$KP.sequential <- m.star

  class(result) <- "bt_robustURN"

  result
}


#' @title
#' MDF procedure for multiple unknown breaks.
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y a time series of interest.
#' @param const Whether level shifts accompany the trend breaks. The
#' deterministic regression always includes a constant and a linear trend.
#' @param breaks Number of breaks.
#' @param breaks.star Number of breaks got from the Kejrival-Perron procedure.
#' @param trim Trimming value for a possible break date bounds.
#' @param ZA Whether ZA variant should be used.
#' @param max.lag Maximum ADF lag. NULL retains the original sample-size rule.
#' @param gls.only If TRUE, compute only the GLS minimum for exactly `breaks`
#' candidate trend breaks, without OLS/union statistics or asymptotic tables.
#' Requires `const = FALSE` and `ZA = FALSE`. This component is used by
#' [uroot.CHLT] for the multiple-break wild bootstrap procedure.
#'
#' @return A list of sublists each containing
#' * The value of statistic: \eqn{MDF-GLS}, \eqn{MDF-OLS},
#' * The asymptotic critical values.
#' \eqn{UR} values are included as well.
#' With `gls.only = TRUE`, returns `statistic`, `break.points`, `lag`, the
#' candidate `grid`, positive `cbar`, `max.lag` and `n.obs` instead.
MDF.mlt <- function(
  y,
  const = FALSE,
  breaks = 1,
  breaks.star = 1,
  trim = 0.15,
  ZA = FALSE,
  max.lag = NULL,
  gls.only = FALSE
) {
  if (!is.logical(gls.only) || length(gls.only) != 1L || is.na(gls.only)) {
    stop("gls.only must be TRUE or FALSE")
  }
  if (gls.only) {
    if (!identical(const, FALSE) || !identical(ZA, FALSE)) {
      stop("gls.only requires const = FALSE and ZA = FALSE")
    }
    checked <- .chlt.validate(y, max.lag, trim, breaks)
    grid <- .chlt.break.grid(nrow(checked$y), breaks, trim)
    cbar <- c(17.6, 21.5, 25.5)[breaks]
    fit <- .chlt.mdf.scan(checked$y, grid, cbar, checked$max.lag)
    return(c(fit, list(grid = grid, cbar = cbar, max.lag = checked$max.lag,
                       n.obs = nrow(checked$y))))
  }

  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- nrow(y)

  ## Critical values ##
  model <- paste0(if (const) "c" else "n", if (ZA) "z" else "")

  cv.DF.OLS.t <- .cval_MDF_multiple[[model]]$cv.DF.OLS.t
  cv.MDF.OLS1 <- .cval_MDF_multiple[[model]]$cv.MDF.OLS1
  cv.MDF.OLS2 <- .cval_MDF_multiple[[model]]$cv.MDF.OLS2
  cv.MDF.OLS3 <- .cval_MDF_multiple[[model]]$cv.MDF.OLS3
  cv.DF.GLS.t <- .cval_MDF_multiple[[model]]$cv.DF.GLS.t
  cv.MDF.GLS1 <- .cval_MDF_multiple[[model]]$cv.MDF.GLS1
  cv.MDF.GLS2 <- .cval_MDF_multiple[[model]]$cv.MDF.GLS2
  cv.MDF.GLS3 <- .cval_MDF_multiple[[model]]$cv.MDF.GLS3
  sap.ur3.ols <- .cval_MDF_multiple[[model]]$sap.ur3.ols
  sap.ur3.olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.olsgls
  sap.ur3.1ols <- .cval_MDF_multiple[[model]]$sap.ur3.1ols
  sap.ur3.1olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.1olsgls
  sap.ur3.2ols <- .cval_MDF_multiple[[model]]$sap.ur3.2ols
  sap.ur3.2olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.2olsgls
  sap.ur3.3ols <- .cval_MDF_multiple[[model]]$sap.ur3.3ols
  sap.ur3.3olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.3olsgls
  sap.ur2.ols <- .cval_MDF_multiple[[model]]$sap.ur2.ols
  sap.ur2.olsgls <- .cval_MDF_multiple[[model]]$sap.ur2.olsgls
  sap.ur2.1ols <- .cval_MDF_multiple[[model]]$sap.ur2.1ols
  sap.ur2.1olsgls <- .cval_MDF_multiple[[model]]$sap.ur2.1olsgls
  sap.ur2.2ols <- .cval_MDF_multiple[[model]]$sap.ur2.2ols
  sap.ur2.2olsgls <- .cval_MDF_multiple[[model]]$sap.ur2.2olsgls
  sap.cv.ur.2 <- .cval_MDF_multiple[[model]]$sap.cv.ur.2
  sap.cv.ur.3 <- .cval_MDF_multiple[[model]]$sap.cv.ur.3

  ## Start ##
  kmax <- if (is.null(max.lag)) trunc(12 * (N / 100)^(1 / 4)) else
    .chlt.max.lag(N, max.lag)

  tb_L <- trunc(trim * N)
  sep <- tb_L
  tb_U <- trunc((1 - trim) * N)

  x <- cbind(.const(N), .trend(N))

  ## GLS case
  r_GLS_t <- GLS.reg(y, x, -13.5)$residuals
  r_OLS_t <- OLS.reg(y, x)$residuals
  k_t <- uroot.ADF(
    r_OLS_t,
    const = FALSE,
    trend = FALSE,
    max.lag = kmax,
    criterion = "aic",
    modified.criterion = TRUE
  )$lag

  DF_GLS_t <- uroot.ADF(
    r_GLS_t,
    const = FALSE,
    trend = FALSE,
    max.lag = k_t,
    criterion = NULL
  )$t.alpha

  ## OLS case
  DF_OLS_t <- uroot.ADF(
    r_OLS_t,
    const = FALSE,
    trend = FALSE,
    max.lag = k_t,
    criterion = NULL
  )$t.alpha

  ## OLS-GLS ##
  ## One break
  MDF_OLS1 <- Inf
  MDF_GLS1 <- Inf

  for (tb1 in tb_L:tb_U) {
    x <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tb1, N) else NULL,
      .dt(tb1, N)
    )

    r_OLS <- OLS.reg(y, x)$residuals
    k_t <- uroot.ADF(
      r_OLS,
      const = FALSE,
      trend = FALSE,
      max.lag = kmax,
      criterion = "aic",
      modified.criterion = TRUE
    )$lag

    DF1 <- uroot.ADF(
      r_OLS,
      const = FALSE,
      trend = FALSE,
      max.lag = k_t,
      criterion = NULL
    )
    MDF_OLS1 <- min(
      MDF_OLS1,
      if (!ZA) {
        denom <- 1 - sum(DF1$model$coefficients) + DF1$alpha
        N * DF1$alpha / denom
      } else {
        DF1$t.alpha
      }
    )

    r_GLS <- GLS.reg(y, x, -17.6)$residuals
    DF1_tb <- uroot.ADF(
      r_GLS,
      const = FALSE,
      trend = FALSE,
      max.lag = k_t,
      criterion = NULL
    )$t.alpha
    MDF_GLS1 <- min(MDF_GLS1, DF1_tb)
  }

  ## Two breaks
  MDF_OLS2 <- Inf
  MDF_GLS2 <- Inf

  for (tb1 in tb_L:(tb_U - sep)) {
    for (tb2 in (tb1 + sep):tb_U) {
      x <- cbind(
        .const(N),
        .trend(N),
        if (const) .du(tb1, N) else NULL,
        .dt(tb1, N),
        if (const) .du(tb2, N) else NULL,
        .dt(tb2, N)
      )

      r_OLS <- OLS.reg(y, x)$residuals
      k_t <- uroot.ADF(
        r_OLS,
        const = FALSE,
        trend = FALSE,
        max.lag = kmax,
        criterion = "aic",
        modified.criterion = TRUE
      )$lag

      DF2 <- uroot.ADF(
        r_OLS,
        const = FALSE,
        trend = FALSE,
        max.lag = k_t,
        criterion = NULL
      )
      MDF_OLS2 <- min(
        MDF_OLS2,
        if (!ZA) {
          denom <- 1 - sum(DF2$model$coefficients) + DF2$alpha
          N * DF2$alpha / denom
        } else {
          DF2$t.alpha
        }
      )

      r_GLS <- GLS.reg(y, x, -21.5)$residuals
      DF2_tb <- uroot.ADF(
        r_GLS,
        const = FALSE,
        trend = FALSE,
        max.lag = k_t,
        criterion = NULL
      )$t.alpha
      MDF_GLS2 <- min(MDF_GLS2, DF2_tb)
    }
  }

  ## Three breaks
  MDF_OLS3 <- Inf
  MDF_GLS3 <- Inf

  for (tb1 in tb_L:(tb_U - 2 * sep)) {
    for (tb2 in (tb1 + sep):(tb_U - sep)) {
      for (tb3 in (tb2 + sep):tb_U) {
        x <- cbind(
          .const(N),
          .trend(N),
          if (const) .du(tb1, N) else NULL,
          .dt(tb1, N),
          if (const) .du(tb2, N) else NULL,
          .dt(tb2, N),
          if (const) .du(tb3, N) else NULL,
          .dt(tb3, N)
        )

        r_OLS <- OLS.reg(y, x)$residuals
        k_t <- uroot.ADF(
          r_OLS,
          const = FALSE,
          trend = FALSE,
          max.lag = kmax,
          criterion = "aic",
          modified.criterion = TRUE
        )$lag

        DF3 <- uroot.ADF(
          r_OLS,
          const = FALSE,
          trend = FALSE,
          max.lag = k_t,
          criterion = NULL
        )
        MDF_OLS3 <- min(
          MDF_OLS3,
          if (!ZA) {
            denom <- 1 - sum(DF3$model$coefficients) + DF3$alpha
            N * DF3$alpha / denom
          } else {
            DF3$t.alpha
          }
        )

        r_GLS <- GLS.reg(y, x, -25.5)$residuals
        DF3_tb <- uroot.ADF(
          r_GLS,
          const = FALSE,
          trend = FALSE,
          max.lag = k_t,
          criterion = NULL
        )$t.alpha
        MDF_GLS3 <- min(MDF_GLS3, DF3_tb)
      }
    }
  }

  ## Alternative break selection
  if (breaks == 2) {
    tbs <- segments.GLS(y, const, TRUE, 2)

    x <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tbs[1], N) else NULL,
      .dt(tbs[1], N),
      if (const) .du(tbs[2], N) else NULL,
      .dt(tbs[2], N)
    )

    tmp.OLS <- OLS.reg(y, x)
    bb <- tmp.OLS$coefficients
    res <- tmp.OLS$residuals
    t_alpha_2_id <- as.numeric(bb[1] / sqrt(drop(t(res) %*% res) / N) > 1)
  }
  if (breaks == 3) {
    tbs <- segments.GLS(y, const, TRUE, 3)

    x <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tbs[1], N) else NULL,
      .dt(tbs[1], N),
      if (const) .du(tbs[2], N) else NULL,
      .dt(tbs[2], N),
      if (const) .du(tbs[3], N) else NULL,
      .dt(tbs[3], N)
    )

    tmp.OLS <- OLS.reg(y, x)
    bb <- tmp.OLS$coefficients
    res <- tmp.OLS$residuals
    t_alpha_3_id <- as.numeric(bb[1] / sqrt(drop(t(res) %*% res) / N) > 1)
  }

  ## breaks.star
  Tbb <- if (breaks.star == 0) {
    0
  } else {
    segments.GLS(y, const, TRUE, breaks.star)
  }

  if (breaks == 2) {
    ur2ols_sa <- as.numeric(
      (DF_OLS_t < (sap.cv.ur.2 * sap.ur2.ols * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS1)) ||
        (MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS2))
    )

    ur2olsgls_sa <- as.numeric(
      (DF_GLS_t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.GLS.t)) ||
        (MDF_GLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS1)) ||
        (DF_OLS_t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS1)) ||
        (MDF_GLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS2))
    )

    ur2_1ols_sa <- as.numeric(
      (MDF_OLS1 < (sap.cv.ur.2 * sap.ur2.1ols * cv.MDF.OLS1)) ||
        (MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.1ols * cv.MDF.OLS2))
    )

    ur2_1olsgls_sa <- as.numeric(
      (MDF_GLS1 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.GLS1)) ||
        (MDF_OLS1 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.OLS1)) ||
        (MDF_GLS2 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.OLS2))
    )

    ur2_2ols_sa <- as.numeric(
      MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.2ols * cv.MDF.OLS2)
    )

    ur2_2olsgls_sa <- as.numeric(
      (MDF_GLS2 < (sap.cv.ur.2 * sap.ur2.2olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.2olsgls * cv.MDF.OLS2))
    )

    UR <- if (breaks.star == 2) {
      ((1 - t_alpha_2_id) * ur2olsgls_sa + t_alpha_2_id * ur2ols_sa)
    } else if (breaks.star == 1) {
      ((1 - t_alpha_2_id) * ur2_1olsgls_sa + t_alpha_2_id * ur2_1ols_sa)
    } else {
      ((1 - t_alpha_2_id) * ur2_2olsgls_sa + t_alpha_2_id * ur2_2ols_sa)
    }

    ## without pre-test    for breaks
    ur2ols_sa <- as.numeric(
      (DF_OLS_t < (sap.cv.ur.2 * sap.ur2.ols * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS1)) ||
        (MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS2))
    )

    ur2olsgls_sa <- as.numeric(
      (DF_GLS_t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.GLS.t)) ||
        (MDF_GLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS1)) ||
        (DF_OLS_t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS1)) ||
        (MDF_GLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS2))
    )

    UR1 <- (1 - t_alpha_2_id) * ur2olsgls_sa + t_alpha_2_id * ur2ols_sa
  } else if (breaks == 3) {
    ur3ols_sa <- as.numeric(
      (DF_OLS_t < (sap.cv.ur.3 * sap.ur3.ols * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS1)) ||
        (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS2)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS3))
    )

    ur3olsgls_sa <- as.numeric(
      (DF_GLS_t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.GLS.t)) ||
        (MDF_GLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS1)) ||
        (DF_OLS_t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS1)) ||
        (MDF_GLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS2)) ||
        (MDF_GLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS3)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS3))
    )

    ur3_1ols_sa <- as.numeric(
      (MDF_OLS1 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS1)) ||
        (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS2)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS3))
    )

    ur3_1olsgls_sa <- as.numeric(
      (MDF_GLS1 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS1)) ||
        (MDF_OLS1 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS1)) ||
        (MDF_GLS2 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS2)) ||
        (MDF_GLS3 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS3)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS3))
    )

    ur3_2ols_sa <- as.numeric(
      (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.2ols * cv.MDF.OLS2)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.2ols * cv.MDF.OLS3))
    )

    ur3_2olsgls_sa <- as.numeric(
      (MDF_GLS2 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.OLS2)) ||
        (MDF_GLS3 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.GLS3)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.OLS3))
    )

    ur3_3ols_sa <- as.numeric(
      MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.3ols * cv.MDF.OLS3)
    )

    ur3_3olsgls_sa <- as.numeric(
      (MDF_GLS3 < (sap.cv.ur.3 * sap.ur3.3olsgls * cv.MDF.GLS3)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.3olsgls * cv.MDF.OLS3))
    )

    UR <- if (breaks.star == 3) {
      ((1 - t_alpha_3_id) * ur3olsgls_sa + t_alpha_3_id * ur3ols_sa)
    } else if (breaks.star == 2) {
      ((1 - t_alpha_3_id) * ur3_1olsgls_sa + t_alpha_3_id * ur3_1ols_sa)
    } else if (breaks.star == 1) {
      ((1 - t_alpha_3_id) * ur3_2olsgls_sa + t_alpha_3_id * ur3_2ols_sa)
    } else {
      ((1 - t_alpha_3_id) * ur3_3olsgls_sa + t_alpha_3_id * ur3_3ols_sa)
    }

    ## without pre-test for breaks

    ur3ols_sa <- as.numeric(
      (DF_OLS_t < (sap.cv.ur.3 * sap.ur3.ols * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS1)) ||
        (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS2)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS3))
    )

    ur3olsgls_sa <- as.numeric(
      (DF_GLS_t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.GLS.t)) ||
        (MDF_GLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS1)) ||
        (DF_OLS_t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.OLS.t)) ||
        (MDF_OLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS1)) ||
        (MDF_GLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS2)) ||
        (MDF_OLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS2)) ||
        (MDF_GLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS3)) ||
        (MDF_OLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS3))
    )

    UR1 <- (1 - t_alpha_3_id) * ur3olsgls_sa + t_alpha_3_id * ur3ols_sa
  }

  result <- list(
    const = const,
    breaks.star = breaks.star,
    breaks.tbb = Tbb,
    breaks = breaks
  )

  result$MDF.GLS.1 <- list(
    statistic = MDF_GLS1,
    cv = cv.MDF.GLS1
  )
  result$MDF.GLS.2 <- list(
    statistic = MDF_GLS2,
    cv = cv.MDF.GLS2
  )
  if (breaks == 3) {
    result$MDF.GLS.3 <- list(
      statistic = MDF_GLS3,
      cv = cv.MDF.GLS3
    )
  }

  result$MDF.OLS.1 <- list(
    statistic = MDF_OLS1,
    cv = cv.MDF.OLS1
  )
  result$MDF.OLS.2 <- list(
    statistic = MDF_OLS2,
    cv = cv.MDF.OLS2
  )
  if (breaks == 3) {
    result$MDF.OLS.3 <- list(
      statistic = MDF_OLS3,
      cv = cv.MDF.OLS3
    )
  }

  result$UR1 <- UR1
  result$UR <- UR

  class(result) <- "bt_mdfHLTN"

  result
}


#' @title
#' Kejrival-Perron procedure of breaks number detection
#'
#' @param y An input series of interest.
#' @param const Whether the break in constant is allowed.
#' @param breaks Number of breaks.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming value for a possible break date bounds.
#'
#' @return The estimated optimal break point.
#'
#' @references
#' Kejriwal, Mohitosh, and Pierre Perron.
#' “A Sequential Procedure to Determine the Number of Breaks in Trend
#' with an Integrated or Stationary Noise Component:
#' Determination of Number of Breaks in Trend.”
#' Journal of Time Series Analysis 31, no. 5 (September 2010): 305–28.
#' https://doi.org/10.1111/j.1467-9892.2010.00666.x.
#'
#' Skrobotov, Anton.
#' “Survey on structural breaks and unit root tests.”
#' Applied Econometrics 58 (2020): 96–141.
#' https://doi.org/10.22394/1993-7601-2020-58-96-141.
breaktest.KP <- function(
  y,
  const = FALSE,
  breaks = 1,
  criterion = "bic",
  trim = 0.15
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- nrow(y)
  kmax <- trunc(12 * (N / 100)^(1 / 4))

  model <- as.numeric(const) + 1
  trim.pos <- which(c(0.01, 0.05, 0.1, 0.15, 0.25) == trim)

  res <- 0

  for (l in 0:(breaks - 1)) {
    test.stat <- KP.seq.statistic(y, const, l, criterion, trim, kmax)
    c.v <- .cval_KP[[model]][[trim.pos]]

    if (test.stat < c.v[2, l + 1]) {
      res <- l
      break
    }
  }

  result <- list(
    breaks = res,
    statistic = test.stat,
    cr.val = c.v[2, res + 1],
    const = const,
    trim = trim
  )
  class(result) <- "bt_KP"

  result
}


#' @title
#' Sequential statistic for breaks at unknown date.
#'
#' @details
#' This procedure is based on ideas of Perron & Yabu (2009).
#'
#' @param y A time series of interest.
#' @param const Allowing the break in constant.
#' @param breaks A number of breaks.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming value for a possible break date bounds.
#' @param max.lag The maximum possible lag in the model.
#'
#' @return An estimated Wald statistic.
#'
#' @references
#' Kejriwal, Mohitosh, and Pierre Perron.
#' “A Sequential Procedure to Determine the Number of Breaks in Trend
#' with an Integrated or Stationary Noise Component:
#' Determination of Number of Breaks in Trend.”
#' Journal of Time Series Analysis 31, no. 5 (September 2010): 305–28.
#' https://doi.org/10.1111/j.1467-9892.2010.00666.x.
#'
#' Skrobotov, Anton.
#' “Survey on structural breaks and unit root tests.”
#' Applied Econometrics 58 (2020): 96–141.
#' https://doi.org/10.22394/1993-7601-2020-58-96-141.
#'
#' @export
KP.seq.statistic <- function(
  y,
  const = FALSE,
  breaks = 1,
  criterion = "aic",
  trim = 0.15,
  max.lag = trunc(12 * (length(y) / 100)^(1 / 4))
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  model <- ifelse(const, 2, 1)
  R <- .cval_PY_sequential[[model]]$R
  v_t <- .cval_PY_sequential[[model]]$v.t

  N <- nrow(y)
  h <- trunc(trim * N)

  if (breaks == 0) {
    datevec <- c(0, N)
  } else {
    SSR.data <- SSR.matrix(y, cbind(.const(N), .trend(N)), h)
    dates <- segments.BP(
      y,
      cbind(.const(N), .trend(N)),
      breaks,
      h,
      SSR.data
    )
    datevec <- c(0, sort(drop(dates$break.point)), N)
  }
  wald <- NULL

  for (i in seq_len(breaks + 1)) {
    T_i <- datevec[i + 1] - datevec[i]
    vect1 <- NULL

    tbL <- max(trunc(datevec[i] + T_i * trim), max.lag + 2)
    tbH <- trunc(datevec[i + 1] - T_i * trim)

    if (tbL < tbH - 1) {
      for (tb in tbL:tbH) {
        lam1 <- tb / datevec[i + 1]

        reg <- cbind(
          .const(N),
          if (const) .du(tb, N) else NULL,
          .trend(N) - datevec[i],
          .dt(tb, N)
        )

        y_i <- .msub(y, (datevec[i] + 1):datevec[i + 1])
        reg_i <- .msub(reg, (datevec[i] + 1):datevec[i + 1])

        khat <- max(1, AR.reg(y_i, reg_i, max.lag, criterion)$lag)

        #u <- OLS.reg(y_i, reg_i)$residuals
        u <- y_i - reg_i %*% solve(t(reg_i) %*% reg_i, t(reg_i) %*% y_i)
        du <- .diffn(u, na = 0)

        regu <- .lagn(u, 1, na = 0)
        for (l in seq_len(khat - 1)) {
          regu <- cbind(regu, .lagn(du, l, na = 0))
        }

        depu <- u[khat:length(u)]
        regu <- .msub(regu, khat:length(u))

        tmp.OLS <- OLS.reg(depu, regu)
        b <- tmp.OLS$coefficients
        ehat <- tmp.OLS$residuals

        VCV <- spdinv(t(regu) %*% regu) * sum(ehat^2) / length(ehat)

        ahat <- b[1]
        vahat <- VCV[1, 1]
        tau1 <- (ahat - 1) / sqrt(vahat)

        # Upper Biased Estimator
        t05 <- v_t[ceiling(lam1 * 10)]

        IP <- trunc((khat + 1) / 2)
        r <- ncol(reg)
        k <- 10

        c1 <- sqrt((1 + r) * T_i)
        c2 <- ((1 + r) * T_i - t05^2 * (IP + T_i)) /
          (t05 * (t05 + k) * (IP + T_i))

        if (tau1 > t05) {
          c.tau <- -tau1
        }
        if (tau1 <= t05 && tau1 > -k) {
          c.tau <- IP * tau1 / N - (r + 1) / (tau1 + c2 * (tau1 + k))
        }
        if (tau1 <= -k && tau1 > -c1) {
          c.tau <- IP * tau1 / N - (r + 1) / tau1
        }
        if (tau1 <= -c1) {
          c.tau <- 0
        }

        amus <- ahat + c.tau * sqrt(vahat)
        if (amus >= 1) {
          amus <- 1
        } else if (amus <= -1) {
          amus <- -0.99
        }

        CR <- sqrt(N) * abs(amus - 1)
        if (CR <= 1) {
          amus <- 1
        }

        gdep <- rbind(
          .msub(y, datevec[i] + 1),
          .msub(y, (datevec[i] + 2):datevec[i + 1]) -
            amus * .msub(y, (datevec[i] + 1):(datevec[i + 1] - 1)) # nolint
        )
        greg <- rbind(
          .msub(reg, datevec[i] + 1),
          .msub(reg, (datevec[i] + 2):datevec[i + 1]) -
            amus * .msub(reg, (datevec[i] + 1):(datevec[i + 1] - 1)) # nolint
        )

        tmp.OLS <- OLS.reg(gdep, greg)
        b <- tmp.OLS$coefficients
        v <- tmp.OLS$residuals

        if (khat == 1) {
          h0 <- sum(v^2) / length(v)
        } else {
          if (amus == 1) {
            regv <- NULL
            for (l in seq_len(khat - 1)) {
              regv <- cbind(regv, .lagn(v, l, na = 0))
            }

            depv <- v[(khat - 1):length(v)]
            regv <- .msub(regv, (khat - 1):length(v))

            tmp.OLS <- OLS.reg(depv, regv)
            beta <- tmp.OLS$coefficients
            e <- tmp.OLS$residuals

            if (!const) {
              h0 <- (sum(e^2) / (T_i - khat)) / ((1 - sum(beta))^2)
            }
            if (const) {
              vbeta <- matrix(0, nrow = khat - 1, ncol = 4)
              for (ki in seq_len(khat - 1)) {
                regki <- cbind(
                  .const(N),
                  .du(tb - ki, N),
                  .trend(N),
                  .du(tb - ki, N) * (.trend(N) - tb)
                )
                gdepki <- rbind(
                  .msub(y, datevec[i] + 1),
                  .msub(y, (datevec[i] + 2):datevec[i + 1]) - # nolint
                    amus *
                      .msub(y, (datevec[i] + 1):(datevec[i + 1] - 1)) # nolint
                )
                gregki <- rbind(
                  reg[datevec[i] + 1, ],
                  regki[(datevec[i] + 2):datevec[i + 1], ] - # nolint
                    amus * regki[(datevec[i] + 1):(datevec[i + 1] - 1), ] # nolint
                )
                vbeta[ki, ] <- drop(OLS.reg(gdepki, gregki)$coefficients)
              }
              sige <- sum(e^2) / (T_i - khat)
              h0 <- sige / ((1 - sum(beta))^2)
              b[2] <- (sqrt(h0) / sqrt(sige)) *
                (b[2] - drop(t(vbeta[, 2]) %*% beta))
            }
          }

          if (abs(amus) < 1) {
            h0 <- .lr.var.quad(v)
          }
        }
        VCV <- h0 * qr.solve(t(greg) %*% greg)
        vect1 <- c(
          vect1,
          t(R %*% b) %*%
            qr.solve(R %*% VCV %*% t(R)) %*%
            (R %*% b)
        )
      }

      wald <- max(wald, log(sum(exp(vect1 / 2)) / T_i))
    }
  }

  wald
}
