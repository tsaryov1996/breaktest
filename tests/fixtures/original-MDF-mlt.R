MDF.mlt <- function(
  y,
  const = FALSE,
  breaks = 1,
  breaks.star = 1,
  trim = 0.15,
  ZA = FALSE
) {
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
  kmax <- trunc(12 * (N / 100)^(1 / 4))

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
