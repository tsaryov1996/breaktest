# breaktest

Due to the fact that this package is not in CRAN yet you need to the following.

1. First, install `devtools`
   ```R
   install.packages("devtools")
   library(devtools)
   ```

2. Clone this repository to any palce convenient for you.
   E.g. `C:\git\breaktest`

3. Go to the folder **in** R
   ```R
   setwd("C:/git/breaktest")
   devtools::install()
   ```

4. Load the package
   ```R
   library(breaktest)
   ```


## Dating a single bubble

`bubble.dates.NBCN()` estimates the emergence, collapse and market recovery dates.
Use `method = "ols"` for ordinary least squares or `method = "wls"` to give
less weight to periods of high volatility. `bubble.dates.HLZ()` estimates the
start and end of the explosive episode from the signs of price changes. Both
methods assume a single bubble and return its boundaries in `break.points`
as observation positions.

```r
ols <- bubble.dates.NBCN(y, method = "ols")
wls <- bubble.dates.NBCN(y, method = "wls")
hlz <- bubble.dates.HLZ(y)
```

If a subsample is too short to estimate a date, NBCN reports `NA` for that
date. When the preliminary OLS fit is incomplete, its available dates are in
`wls$pilot$break.points`.

## Optional endpoint trimming for bubble confidence sets

The three `cs.bubble.*()` functions accept `trim_endpoints = TRUE` to apply
5% trimming of estimated regime windows, as described in
[Kurozumi and Skrobotov (2026), Section 4](https://arxiv.org/html/2511.16172v2#S4).
The default is `FALSE`. The existing numeric `trim` controls comparison distances.

```r
cs_e <- cs.bubble.emerge(
  y[seq_len(Tc + 1L)], phi_a = fit$phi_a, s2 = fit$s.sq,
  trim_endpoints = TRUE
)
```

Pass the original window, including its lag observation. Returned vectors keep
their original lengths and indices, with zeros at excluded dates. See
`?bubbles.nbcn` for all three functions, the scaling convention and a runnable
example.

After installing the package, run `Rscript --vanilla tests/bubble-trim.R`.
