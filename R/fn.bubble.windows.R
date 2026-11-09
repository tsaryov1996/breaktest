# Prepare one of the two-regime windows from Kurozumi and Skrobotov, Section 4.
# The first element is always the lag observation, not a candidate break date.
.bubble.trim.window <- function(y, trim, trim_endpoints, event) {
  if (!is.logical(trim_endpoints) || length(trim_endpoints) != 1L ||
      is.na(trim_endpoints)) {
    stop("'trim_endpoints' must be TRUE or FALSE.")
  }
  is_column <- is.null(dim(y)) || (is.matrix(y) && ncol(y) == 1L)
  if (!is.numeric(y) || !is_column || length(y) < 3L || any(!is.finite(y))) {
    stop("'y' must contain at least three finite numeric observations, ",
         "as a vector or one-column matrix (including the initial lag).")
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0) {
    stop("'trim' must be a nonnegative finite number.")
  }
  original_N <- length(y) - 1L
  k <- if (trim_endpoints) floor(0.05 * original_N) else 0L
  left <- if (event == "emerge") 0L else k
  right <- if (event == "recovery") 0L else k
  N <- original_N - left - right
  Teps <- max(1L, floor(N * trim))
  if (Teps >= N || (event == "recovery" && 2L * Teps > N)) {
    stop("The retained window is too short for the selected 'trim'.")
  }
  list(
    y = matrix(y[seq.int(left + 1L, original_N - right + 1L)], ncol = 1L),
    original_N = original_N, left = left, right = right
  )
}

# Keep membership entries on the input window's time grid. Excluded dates are 0.
.bubble.pad.sets <- function(sets, window) {
  if (window$left == 0L && window$right == 0L) return(sets)
  lapply(sets, function(x) c(integer(window$left), x, integer(window$right)))
}
