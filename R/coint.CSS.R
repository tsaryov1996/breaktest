#' @title
#' KPSS-based stationarity and cointegration test with multiple structural breaks
#'
#' @description
#' Test for stationarity and cointegration in the presence of multiple structural breaks using the KPSS test.
#'
#' This procedure extends the single- and two-break tests of Carrion-i-Silvestre
#' and Sansó (2006, 2007) to the case of an arbitrary (known) number of breaks.
#'
#' @details
#' This function is a generalization of two methods proposed in Carrion-i-Silvestre and Sansó (2006, 2007).
#'
#' Carrion-i-Silvestre and Sansó (2006) proposed an LM-Type statistic
#' to test the null hypothesis of cointegration allowing for the possibility of a structural break,
#' in both the deterministic component and the cointegration vector.
#' The test has been designed to be used as a complement to the usual non-cointegration tests
#' in order to obtain stronger evidence of cointegration.
#' The cases of known and unknown break dates are considered.
#'
#' Carrion-i-Silvestre and Sansó (2007) generalized the KPSS-type test for stationarity
#' to allow for two structural breaks.
#' Seven models have been defined depending on the way in which the structural breaks
#' affect the time series behaviour.
#'
#' At the first step, the residuals are obtained from a regression of the series \eqn{y}
#' on the set of deterministic components and (if any) exogenous regressors \eqn{x}.
#' Deterministic components may include
#' * a constant,
#' * a trend,
#' * constant with breaks \eqn{{DU}_{t,\tau} = \mathrm{1}[t > \tau]},
#' * trends with breaks \eqn{{DT}_{t,\tau} = \mathrm{1}[t > \tau](t - \tau)}.
#'
#' If \eqn{x} is weakly exogenous, the regression is estimated by OLS; otherwise, DOLS regression is applied.
#' The procedure also allows for the break in the cointegrating equation;
#' in that case, products of \eqn{{DU}_{t,\tau}} and \eqn{x} are added to the regressors list.
#'
#' After the residuals are obtained, the long-run variance of the KPSS test statistic is calculated using a kernel
#' specified via `lr.kernel`. By default (`lr.kernel = "Kurozumi"`), the QS kernel is used
#' with bandwidth selected as in Andrews (1991) and trancated following Kurozumi's (2002) proposal.
#'
#' The p-value (if needed) is calculated by a bootstrapping procedure following Cavaliere and Taylor (2006).
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' If it's NULL then the result of this function is just a KPSS test for y.
#' @param const,trend Whether a constant or trend should be included.
#' @param break.type A single value or vector of
#' * `c`: for the break in const,
#' * `t`: for the break in trend,
#' * `ct`: for the break in const and trend.
#' * , 'NULL' if there are no breaks
#' @param break.point Array of structural break moments, 'NULL' if there are no breaks.
#' @param break.coint Whether breaks in cointegrating relation are needed to be included.
#' See Carrion-i-Silvestre & Sansó (2006) for details.
#' @param weakly.exog Boolean where we specify whether the stochastic regressors are exogenous or not
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous (DOLS is used in this case).
#' @param max.lags,max.leads Scalars defininig the initial number of lags and leads for DOLS.
#' @param lr.kernel Kernel for calculating long-run variance
#' * `Kurozumi` for the Kurozumi's proposal, using Quadratic kernel (default).
#' * `Bartlett`: for Bartlett kernel,
#' * `Quadratic`: for Quadratic Spectral kernel,
#' @param lr.lag scalar showing the bandwidth of long run variance estimator.
#' If negative or `NULL` then the bandwidth is selected as in Andrews (1991).
#' @param criterion Information criterion for DOLS lags and leads selection: aic, bic, hq, or lwz,
#' @param boot.p Whether bootstrapped p-values should be returned.
#' @param boot.iter The number of bootstrap iterations,
#' @param boot.type Method of building synthetic \eqn{y}:
#' * `sample` ---  permutations bootstrap,
#' * `Cavaliere-Taylor` --- wild bootstrap with \eqn{N(0, 1)} multipliers,
#' * `Rademacher` --- wild bootstrap with \eqn{\pm 1} multipliers with equal probabilities.
#' @param ... A dummy parameter for technical purposes. Just ignore it.
#'
#' If you want to get results like in Carrion-i-Silvestre & Sansó (2006) then you should
#' * provide a `break.point` of length 1;
#' * for the specification of `An` set `const=TRUE`, `trend=FALSE`, `break.type=c("c")`;
#' * for the specification of `A` set `const=TRUE`, `trend=TRUE`, `break.type=c("c")`;
#' * for the specification of `B` set `const=TRUE`, `trend=TRUE`, `break.type=c("t")`;
#' * for the specification of `C` set `const=TRUE`, `trend=TRUE`, `break.type=c("ct")`;
#' * for the specification of `D` set `const=TRUE`, `trend=FALSE`, `break.type=c("c")`;
#' * for the specification of `E` set `const=TRUE`, `trend=TRUE`, `break.type=c("ct")`, and `break.coint=TRUE`.
#'
#' If you want to get results like in Carrion-i-Silvestre & Sansó (2007) then you should
#' * provide a `break.point` of length 2;
#' * for the specification of `AAn` set `const=TRUE`, `trend=FALSE`, `break.type=c("c", "c")`;
#' * for the specification of `AA` set `const=TRUE`, `trend=TRUE`, `break.type=c("c", "c")`;
#' * for the specification of `BB` set `const=TRUE`, `trend=TRUE`, `break.type=c("t", "t")`;
#' * for the specification of `CC` set `const=TRUE`, `trend=TRUE`, `break.type=c("ct", "ct")`;
#' * for the specification of `AB-BA` set `const=TRUE`, `trend=FALSE`, `break.type=c("c", "t")`;
#' * for the specification of `AC-CA` set `const=TRUE`, `trend=TRUE`, `break.type=c("c", "ct")`;
#' * for the specification of `BC-CB` set `const=TRUE`, `trend=TRUE`, `break.type=c("t", "ct")`.
#'
#' @return An object of class `bt_kpss` containing
#' * `statistic`: the value of test statistic,
#' * `coefficients`: OLS/DOLS estimates of the coefficients,
#' * `se.coefs`: standard errors of the coefficients above,
#' * `t.stats`: \eqn{t}-statistics for the coefficients above,
#' * `residuals`: Residuals of the model,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `fitted.values`: fitted values of the estimated model,
#' * `endog`: final \eqn{y} vector used in the estimated model,
#' * `exog`: final \eqn{x} matrix used in the estimated model,
#' * `DOLS.lags`: The estimated number of lags and leads in DOLS,
#' * `break.type`, `break.point`, `break.coint`: breaks specification,
#' * `criterions`: values of the information criterions for the estimated model,
#' * `lags`, `leads`: number of lags and leads in the estimated model.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' Cavaliere, G., & Robert Taylor, A. M. (2006).
#' “Testing the Null of Co‐integration in the Presence of Variance Breaks.”
#' Journal of Time Series Analysis, 27(4), 613-636.
#'
#' Andrews, Donald W. K.
#' “Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation.”
#' Econometrica 59, no. 3 (1991): 817–58.
#' https://doi.org/10.2307/2938229.
#'
#' Kurozumi, Eiji.
#' “Testing for Stationarity with a Break.”
#' Journal of Econometrics 108, no. 1 (May 1, 2002): 63–99.
#' https://doi.org/10.1016/S0304-4076(01)00106-3.
#'
#' @export
coint.CSS <- function(
  y,
  x = NULL,
  const = TRUE,
  trend = FALSE,
  break.type,
  break.point,
  break.coint = FALSE,
  weakly.exog = TRUE,
  max.lags,
  max.leads,
  criterion = "bic",
  lr.kernel = "Kurozumi",
  lr.lag = NULL,
  boot.p = FALSE,
  boot.type = "sample",
  boot.iter = 999
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.null(x) && !is.matrix(x)) {
    x <- as.matrix(x)
  }

  N <- nrow(y)

  deter <- trend.variables(break.type, N, break.point, const, trend)
  xdu <- NULL
  if (!is.null(x) && break.coint) {
    for (bp in break.point) {
      xdu <- cbind(xdu, sweep(x, 1, .du(bp, N), `*`))
    }
  }

  result <- if (weakly.exog || is.null(x)) {
    mX <- cbind(deter, xdu, x)
    result <- OLS.reg(y, mX)
    c(
      result,
      list(
        break.type = break.type,
        break.point = break.point,
        break.coint = break.coint,
        lags = 0,
        leads = 0
      )
    )
  } else {
    z <- cbind(deter, xdu)
    result <- DOLS.reg(
      y,
      x,
      z,
      max.lags,
      max.leads,
      criterion
    )
    c(
      result,
      list(
        break.type = break.type,
        break.point = break.point,
        break.coint = break.coint
      )
    )
  }

  test <- if (lr.kernel == "Kurozumi") {
    .kpss.statistic(result$residuals, .lr.var.kurozumi(result$residuals))
  } else {
    .kpss.statistic(
      result$residuals,
      .lr.var.spc(result$residuals, lr.lag, lr.kernel)
    )
  }

  result$statistic <- test
  result$lr.var.max.lag <- lr.lag
  result$lr.var.kernel <- lr.kernel
  class(result) <- "bt_kpss"

  if (boot.p) {
    result$p.value <- bootstrap(result, boot.iter, boot.type)
  }

  result
}
