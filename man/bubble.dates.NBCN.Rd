% Please edit documentation in R/bubble.dates.NBCN.R
\name{bubble.dates.NBCN}
\alias{bubble.dates.NBCN}
\title{Bubble dates with optional volatility correction}
\usage{
bubble.dates.NBCN(
  y,
  trim = 0.05,
  method = c("wls", "ols"),
  bandwidth = NULL,
  bandwidth.grid = NULL
)
}
\arguments{
\item{y}{A finite numeric vector, univariate time series or one-column
matrix containing the initial level and at least eight observations.}
\item{trim}{Fraction of the full sample excluded at each end of a search
window, strictly between zero and one quarter.}
\item{method}{Use adaptive weighted least squares (\code{"wls"}) or ordinary
least squares (\code{"ols"}).}
\item{bandwidth}{Positive Gaussian bandwidth as a fraction of the full
sample size. NULL selects it by leave-one-out cross-validation.}
\item{bandwidth.grid}{Positive bandwidth candidates, used only when
\code{bandwidth} is NULL. The default is 25 logarithmically spaced values
from \code{1 / (length(y) - 1)} to \code{0.5}. Bandwidth arguments require WLS.}
}
\description{
Estimate the emergence, collapse and recovery of one bubble using
sequential OLS or adaptive WLS sample splitting.
}
\details{
With input \eqn{y_0,\ldots,y_T}, \code{y[1]} is the initial level. A date
\eqn{k} is the last equation before a regime change. \code{Te_est}, \code{Tc_est}
and \code{Tr_est} use this convention; \code{break.points} adds one to give row
indices in \code{y}, suitable for \code{\link[=plot_bubble]{plot_bubble()}}.

Each split minimises the weighted residual sum of squares of two AR(1)
regressions without intercepts. Collapse is estimated on the full sample;
emergence on equations \code{1:Tc_est}; recovery on \code{(Tc_est + 1):T}, keeping
the lag at \code{Tc_est}. For a window \eqn{(a,b]}, candidates run from
\eqn{a+\max(h+1,2)} to \eqn{b-\max(h,2)}, where
\eqn{h=\lfloor T\,trim\rfloor}. This follows the paper's recovery
search convention and requires at least two equations in each regime.
Exact ties select the earliest date.

WLS first obtains OLS dates. In the pilot residual regression,
\code{diff(y)} has an intercept and a lag coefficient in each of the explosive
and collapsing regimes, and fitted value zero outside them (equation 4
of Kurozumi and Skrobotov). Gaussian leave-one-out smoothing of squared
residuals supplies inverse-variance weights for all three searches.
CV minimises the squared prediction error of the squared residuals;
the explicit criterion and default grid are implementation choices since
the paper does not specify them. Exact CV ties select the smallest
bandwidth. Extremely small bandwidths use the nearest-neighbour kernel
limit if every Gaussian weight underflows. No variance floor is added.

If a split leaves too little data or cannot identify emergence or recovery,
that date is NA and the other dates are retained (\code{status = "partial"}).
If the OLS pilot is incomplete, WLS cannot estimate its variance:
all final dates are NA (\code{status = "incomplete.pilot"}) and the available
OLS dates remain in \code{pilot}. An unidentified collapse, a rank-deficient
pilot regression or undefined variances cause an error.

These are conditional point estimates for one four-regime bubble, not a
bubble-existence test or confidence sets. The model assumes independent
innovations with a positive deterministic variance function. WLS need not
outperform OLS for every volatility path, and its dates do not validate
weighted confidence sets from \code{\link[=cs.bubble.emerge]{cs.bubble.emerge()}}
and related functions.
}
\value{
An object of class \code{bt_bubbleDates} with:
\itemize{
\item \code{method}, \code{status}: requested method and \code{"complete"},
\code{"partial"} or \code{"incomplete.pilot"}; \code{reason} explains missing dates.
\item \code{Te_est}, \code{Tc_est}, \code{Tr_est}: equation indices;
\code{break.points}: corresponding row indices named \code{emerge},
\code{collapse}, \code{recovery}. Missing dates are NA.
\item \code{coefficients}: four AR coefficients at the selected dates, or NA
where a regime boundary is missing; \code{phi_a}, \code{phi_b}: the two middle
coefficients. \code{regime.check} checks \code{phi_a > 1} and
\code{0 < phi_b < 1}, not significance.
\item \code{pilot}: preliminary OLS dates, row indices, status and reasons;
for WLS, also residuals and coefficients if its residual regression was fitted.
\item \code{variance}, \code{bandwidth}, \code{bandwidth.grid}, \code{cv}:
WLS smoothing results, NULL when unused. CV scores use residuals divided by
\code{scale}.
\item \code{search}: final candidate-date/objective profiles, NULL if WLS was
not run. Objectives use \code{y / scale} and inverse variances divided by
their maximum.
\item \code{y}, \code{trim}, \code{n.obs}, \code{scale}: data, trimming fraction,
number of equations and numerical scaling factor.
}
}
\examples{
set.seed(42)
y <- numeric(201)
y[1] <- 100
phi <- rep(c(1, 1.04, 0.94, 1), c(70, 50, 40, 40))
for (t in seq_along(phi)) y[t + 1] <- phi[t] * y[t] + rnorm(1)
bubble.dates.NBCN(y, bandwidth = 0.1)$break.points
}
\references{
Kurozumi, E. and A. Skrobotov (2026). Improving the accuracy of bubble date
estimators under time-varying volatility. arXiv:2306.02977v2, Sections 3--6.
\url{https://arxiv.org/abs/2306.02977v2}.
}
