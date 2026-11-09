# Run against the installed package:
# LC_ALL=en_US.UTF-8 Rscript --vanilla tests/bubble-trim.R
library(breaktest)

checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
  checks <<- checks + 1L
}
expect_error <- function(expr, pattern, label) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  check(inherits(err, "error") && grepl(pattern, conditionMessage(err)), label)
}
args_for <- function(event, y) {
  switch(event,
         emerge = list(y = y, phi_a = 1.03, s2 = 2),
         collapse = list(y = y, lambda_e = 0.25, phi_a = 1.03,
                         phi_b = 0.97, s2 = 2),
         recovery = list(y = y, lambda_e = 0.25, lambda_c = 0.50,
                         phi_a = 1.03, phi_b = 0.97, s2 = 2))
}
run <- function(event, y, trim_endpoints = FALSE, trim = 0.1) {
  fun <- getExportedValue("breaktest", paste0("cs.bubble.", event))
  do.call(fun, c(args_for(event, y),
                 list(trim = trim, trim_endpoints = trim_endpoints)))
}
events <- c("emerge", "collapse", "recovery")
series <- function(n) {
  # A finite, nonzero series with a rise and a decline, plus irregular noise.
  t <- 0:n
  10 + sin(t / 4) + 0.04 * t + 2 * sin(t / 13)
}

# The optional flag preserves the public default and accepts column matrices.
for (event in events) {
  fun <- getExportedValue("breaktest", paste0("cs.bubble.", event))
  y <- series(61)
  check(identical(do.call(fun, args_for(event, y)), run(event, y, FALSE)),
        paste(event, "default differs from FALSE"))
  check(identical(run(event, y, TRUE), run(event, matrix(y), TRUE)),
        paste(event, "vector and column-matrix results differ"))
  # These larger candidate gaps were already supported before endpoint trimming.
  gaps <- if (event == "recovery") 0.5 else c(0.5, 0.6)
  for (gap in gaps) {
    original_call <- do.call(fun, c(args_for(event, y), list(trim = gap)))
    check(identical(original_call, run(event, y, FALSE, trim = gap)),
          paste(event, "changed the default at trim", gap))
    check(all(lengths(run(event, y, TRUE, trim = gap)) == 61L),
          paste(event, "rejected a valid candidate gap", gap))
  }
}

# Explicit expected sample positions test the article's 5% endpoint rule. The
# helper is observed at the statistical boundary, without calling its window
# construction helper to produce expected values.
for (n in c(100L, 101L)) {
  y <- seq_len(n + 1L) + sin(seq_len(n + 1L)) / 10
  for (event in events) {
    observed <- getExportedValue("breaktest", paste0("cs.bubble.", event))
    observer_env <- new.env(parent = asNamespace("breaktest"))
    environment(observed) <- observer_env
    expected_rows <- switch(event,
                            emerge = seq.int(1L, n - 4L),
                            collapse = seq.int(6L, n - 4L),
                            recovery = seq.int(6L, n + 1L))
    base_stat <- get(".nbcn.t.stat", envir = asNamespace("breaktest"))
    calls <- 0L
    seen_lag <- FALSE
    seen_last <- FALSE
    observer_env$.nbcn.t.stat <- function(x, t1, t2) {
      if (!identical(as.numeric(x), y[expected_rows])) {
        stop("Statistical calculation received observations outside its window")
      }
      if (min(t1, t2) < 2L || max(t1, t2) > length(x)) {
        stop("Statistical calculation indexed outside its window")
      }
      calls <<- calls + 1L
      seen_lag <<- seen_lag || min(t1, t2) == 2L
      seen_last <<- seen_last || max(t1, t2) == length(x)
      base_stat(x, t1, t2)
    }
    trimmed <- do.call(observed, c(args_for(event, y), list(trim_endpoints = TRUE)))
    plain <- run(event, y, FALSE)
    check(calls > 0L && seen_lag && seen_last,
          paste(event, n, "did not exercise both retained sample edges"))
    check(identical(names(trimmed), names(plain)), paste(event, "changed method names"))
    check(all(lengths(trimmed) == n), paste(event, "changed output lengths"))
    check(all(vapply(trimmed, function(z) all(z %in% 0:1), logical(1))),
          paste(event, "returned nonbinary membership"))
    excluded <- switch(event,
                       emerge = seq.int(n - 4L, n),
                       collapse = c(1:5, seq.int(n - 4L, n)),
                       recovery = 1:5)
    check(all(vapply(trimmed, function(z) all(z[excluded] == 0L), logical(1))),
          paste(event, "included trimmed dates"))
  }
}

# floor(.05 * N) == 0 must be an exact no-op, not merely similar estimates.
for (event in events) {
  for (n in c(8L, 19L)) {
    check(identical(run(event, series(n), FALSE),
                    run(event, series(n), TRUE)),
          paste(event, "changed results when no observations are removed"))
  }
  for (bad in list(NA, c(TRUE, FALSE), 1, "TRUE", logical(0), NULL)) {
    args <- args_for(event, series(30))
    args["trim_endpoints"] <- list(bad)
    expect_error(do.call(getExportedValue("breaktest", paste0("cs.bubble.", event)), args),
                 "trim_endpoints", paste(event, "accepted invalid logical flag"))
  }
  expect_error(run(event, c(1, 2), TRUE), "three finite",
               paste(event, "accepted an unusable short sample"))
  expect_error(run(event, series(30), TRUE, trim = 1), "trim",
               paste(event, "accepted an invalid candidate-gap trim"))
}

# Numerical regression at one collapse candidate. Compute the AR numerator
# directly in original transition coordinates, independently of the runtime
# window/helper implementation. Choose s2 so that omitting the left-time offset
# or replacing original_N by retained_N changes the binary test decision.
n <- 100L
t <- 0:n
y <- 100 * 1.03^pmin(t, 50) * 0.97^pmax(t - 50, 0)
candidate <- 70L
alternatives <- 79:95  # Retained interval [6,95], gap floor(.1*90)=9.
raw_t <- vapply(alternatives, function(v) {
  tt <- seq.int(candidate, v)
  lagged <- y[tt]
  sum(lagged * (y[tt + 1L] - lagged)) / sqrt(sum(lagged^2))
}, numeric(1))
cv <- -sqrt(0.25 * qchisq(0.05, df = 1))
scale <- sqrt(100 * 1.03^(2 * candidate) * 0.03 / 2)
target <- 0.98 * cv
s2 <- (mean(raw_t) / (target * scale))^2
collapse <- cs.bubble.collapse(y, 0.25, 1.03, 0.97, s2,
                               trim_endpoints = TRUE)
check(collapse$cset_EMa12[candidate] == 0L,
      "Collapse scaling lost the original phase origin or original sample scale")
check(target * 1.03^5 < cv && target * sqrt(100 / 90) < cv,
      "Numerical fixture does not discriminate missing offset/sample scaling")

# A second numerical regression verifies the recovery phase offset. Its
# left-sided extremum has a denominator dependent on the selected alternative.
candidate <- 70L
alternatives <- 14:61  # [offset + Teps, candidate - Teps] with offset=5, Teps=9.
# Use a declining series throughout this check so the relevant extremum has
# the same sign as the one-sided critical value.
y_recovery <- 100 * 0.97^t
raw_t <- vapply(alternatives, function(v) {
  tt <- seq.int(v, candidate)
  lagged <- y_recovery[tt]
  sum(lagged * (y_recovery[tt + 1L] - lagged)) / sqrt(sum(lagged^2))
}, numeric(1))
selected <- alternatives[which.max(raw_t)]
scale <- 1.03^25 * 0.97^selected * sqrt(100 * 0.03 / 2)
target <- 1.05 * cv
s2 <- (max(raw_t) / (target * scale))^2
recovery <- cs.bubble.recovery(y_recovery, 0.25, 0.50, 1.03, 0.97, s2,
                               trim_endpoints = TRUE)
check(recovery$cset_EMb21s[candidate] == 1L,
      "Recovery scaling lost the original collapse-time origin")
check(target * 0.97^5 > cv,
      "Recovery fixture does not discriminate a missing left-time offset")
target <- 0.98 * cv
s2 <- (max(raw_t) / (target * scale))^2
recovery <- cs.bubble.recovery(y_recovery, 0.25, 0.50, 1.03, 0.97, s2,
                               trim_endpoints = TRUE)
check(recovery$cset_EMb21s[candidate] == 0L,
      "Recovery trimming changed the original explosive duration/sample scale")
check(target * 1.03^2 < cv && target * sqrt(100 / 95) < cv,
      "Recovery fixture does not discriminate changed explosive duration/sample scale")

# Trimming must actually alter the calculation, not just add a silent flag.
changed <- vapply(events, function(event) {
  !identical(run(event, series(100), FALSE),
             run(event, series(100), TRUE))
}, logical(1))
check(all(changed), "Trimming had no observable effect for at least one event")

# Plotting's pre-existing explicit input-to-row mapping must continue working.
# This is optional only because the calculation tests require no plotting deps.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  y <- series(200)
  sets <- list(
    Te = run("emerge", y[1:101], TRUE)$cset_LE,
    Tc = run("collapse", y[41:161], TRUE)$cset_EMa,
    Tr = run("recovery", y[101:201], TRUE)$cset_LE
  )
  plot <- plot_bubble_confsets(y, sets,
    break_points = c(Te = 41, Tc = 101, Tr = 161),
    set_index = list(Te = 2:101, Tc = 42:161, Tr = 102:201))
  check(inherits(plot, "ggplot"), "Trimmed sets no longer map to plotting rows")
  check(length(ggplot2::ggplot_build(plot)$data) > 0,
        "Plot with trimmed sets did not build")
}
cat(sprintf("PASS: %d bubble endpoint trimming checks\n", checks))
