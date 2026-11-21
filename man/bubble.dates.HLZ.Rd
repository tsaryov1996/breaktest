% Please edit documentation in R/bubble.dates.HLZ.R
\name{bubble.dates.HLZ}
\alias{bubble.dates.HLZ}
\title{Sign-based bubble date estimation}
\usage{
bubble.dates.HLZ(y, min.window = 0.1, epsilon = 0.01, trim = 0)
}
\arguments{
\item{y}{A finite numeric series or one-column numeric matrix. Observations
must be consecutive; missing values are not removed.}

\item{min.window}{Minimum window length as a fraction of \code{length(y)}.}

\item{epsilon}{Power in the modified residual variance, in \eqn{(0,1]}.}

\item{trim}{Fraction excluded at each outer end of the search. The default
allows both sample boundaries; \code{trim = 0.05} excludes 5 percent at each end.}
}
\value{
An object of class \code{bt_bubbleDates}: a list with \code{y}, \code{method}, the
tuning parameters, \code{n.obs}, \code{min.length}, \code{break.points} (emergence and
collapse), \code{window}, the maximum \code{statistic}, \code{C.t}, and a \code{profile} of the
best window ending at each candidate endpoint. \code{diagnostics} gives the
numbers of candidate, singular, nonpositive-variance and valid windows,
the number of zero increments, and whether emergence precedes the sample.
}
\description{
Estimate the emergence and collapse of one explosive episode using the
modified sign-based criterion in Section 6 of Harvey, Leybourne and Zu (2020).
}
\details{
For \eqn{N=length(y)}, form the cumulative signs
\eqn{C_1=0}, \eqn{C_t=\sum_{i=2}^t[1-2I(y_i-y_{i-1}\leq0)]}.
Zero increments have sign -1. On each window \eqn{(a,b]}, regress
\eqn{\Delta C_t} on \eqn{C_{t-1}} without an intercept, using the same
cumulative series throughout. The dates maximize
\deqn{\hat\rho(a,b)\sqrt{\sum_{t=a+1}^b C_{t-1}^2}/
       [\tilde s^2(a,b)]^{\epsilon/2},}
where \eqn{\tilde s^2(a,b)=[b\hat s^2(0,b)-a\hat s^2(0,a)]/(b-a-1)}
and \eqn{\hat s^2(0,j)=SSR(0,j)/(j-1)}.

Candidate boundaries satisfy \eqn{a\geq ceiling(N\,trim)},
\eqn{b\leq floor(N(1-trim))}, and
\eqn{b-a\geq max(3,ceiling(N\,min.window))}.
The initial row has \eqn{C_0=\Delta C_1=0}; prefix terms for
\eqn{a=0,1} are zero. A prefix with a zero regressor uses coefficient zero.
Windows with zero regressor sum of squares or nonpositive modified variance
are excluded and counted in \code{diagnostics}. If none remain, the function
errors. Exact ties select the earliest end, then the earliest start.

\code{break.points} gives boundaries in the rows of \code{y}: emergence is the last
row before the explosive window, and collapse is its last row. An emergence
of zero means before the first observed row; it has no observation or calendar
date. No index offset is needed. The \code{window} is open at its start and closed
at its end.

These are conditional point estimates for one bubble, with no recovery date
or confidence sets. The function does not test whether a bubble exists.
The theory assumes independent innovations with zero mean and median,
positive deterministic time-varying volatility, and the other regularity
conditions of HLZ. The minimum window must not exceed the true episode length.
}
\references{
Harvey, D. I., S. J. Leybourne, and Y. Zu (2020).
Sign-Based Unit Root Tests for Explosive Financial Bubbles in the Presence
of Deterministically Time-Varying Volatility. Econometric Theory 36, 122--169.
\url{https://doi.org/10.1017/S0266466619000057}.
}
\examples{
set.seed(42)
y <- cumsum(rnorm(100))
fit <- bubble.dates.HLZ(y)
fit$break.points
}
