#' Confidence sets for the emergence, collapse, and recovery date of a bubble
#' @name bubbles.nbcn
#'
#' @param y      Time series of interest. Should be a column matrix; if not is converted internally.
#' @param phi_a    AR coefficient before collapse (explosive phase).
#' @param phi_b    AR coefficient after collapse (recovery phase).
#' @param lambda_e Scalar, fraction for critical value scaling.
#' @param lambda_c  Scalar lambda_c = Tc / T_yall.
#' @param s2     Residual variance.
#' @param trim   Trimming parameter to exclude observations from breakpoint candidates.
#'
#' @details
#' Both the "estimated" and "true" parameter variants are handled via the `phi_a`, `phi_b`, and `s2` arguments,
#' which are to be taken from the results of [segments.AR1] or [segments.NBCN] functions.
#'
#' @examples
#' fit <- segments.NBCN(y, trim = 0.1)
#' estimated_dates <- c(Te = fit$Te_est, Tc = fit$Tc_est, Tr = fit$Tr_est)
#' Te <- unname(estimated_dates["Te"])
#' Tc <- unname(estimated_dates["Tc"])
#' Tr <- unname(estimated_dates["Tr"])
#' cs_e <- cs.bubble.emerge(
#'   y = y[seq_len(Tc + 1L)], 
#'   phi_a = fit$phi_a, s2 = fit$s.sq, trim = 0.1
#' )
#' cs_c <- cs.bubble.collapse(
#'   y = y[seq.int(Te + 1L, Tr + 1L)], lambda_e = Te / n, 
#'   phi_a = fit$phi_a, phi_b = fit$phi_b, s2 = fit$s.sq, trim = 0.1
#' )
#' cs_r <- cs.bubble.recovery(
#'   y = y[seq.int(Tc + 1L, n + 1L)], 
#'   lambda_e = Te / n, lambda_c = Tc / n,
#'   phi_a = fit$phi_a, phi_b = fit$phi_b, s2 = fit$s.sq, trim = 0.1
#' )
#' conf_sets <- list(Te = cs_e$cset_LE, Tc = cs_c$cset_EMa, Tr = cs_r$cset_LE)
#' set_index <- list(
#'   Te = seq.int(2L, Tc + 1L),
#'   Tc = seq.int(Te + 2L, Tr + 1L),
#'   Tr = seq.int(Tc + 2L, n + 1L)
#' )
#' break_points <- estimated_dates + 1L
#' p_bands <- plot_bubble_confsets(
#'   y, conf_sets = conf_sets, break_points = break_points,
#'   set_index = set_index, time = dates, shade_phases = TRUE,
#'   date_format = "%Y-%m", y_label = "Time series",
#'   title = "Confidence sets for bubble dates"
#' )
#' p_rows <- plot_confsets(
#'   y, conf_sets = conf_sets, break_points = break_points,
#'   set_index = set_index, time = dates, relative_heights = c(2, 1),
#'   date_format = "%Y-%m", x_label = "Date", y_label = "Time series",
#'   title = "Confidence sets for bubble dates"
#' )
#' 
#' Kurozumi, Eiji, and Anton Skrobotov. 2026.
#' "Confidence Sets for the Emergence, Collapse, and Recovery Dates of a Bubble".
#' arXiv:2511.16172.
#' Preprint, arXiv. https://doi.org/10.48550/arXiv.2511.16172.
#' @details


#' @rdname bubbles.nbcn
#' @order 1
#' @return Named list of binary vectors, where 1 means that this observation is in the corresponding confidence set.
#' @export
cs.bubble.emerge <- function(
  y,
  phi_a,
  s2,
  trim = 0.1
) {
  y <- as.matrix(y)

  N <- length(y) - 1
  Teps <- max(1L, floor(N * trim))
  rho_a <- phi_a - 1

  cset_EMa12 <- integer(N)
  cset_EMb12 <- integer(N)
  cset_LRa12 <- integer(N)
  cset_LRb12 <- integer(N)
  cset_EMa21 <- integer(N)
  cset_EMb21 <- integer(N)
  cset_LRa21 <- integer(N)

  # --- T1 < T2 ---
  for (T1 in 2:(N + 1 - Teps)) {
    T2_range <- (T1 + Teps):(N + 1)
    tstat <- numeric(length(T2_range))
    statLRa <- numeric(length(T2_range))
    statLRb <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      T2 <- T2_range[k]
      rv <- .nbcn.t.stat(y, T1, T2)

      tstat[k] <- rv$t.stat / sqrt(s2)
      statLRa[k] <- (2 * rv$S_y1dy - rho_a * rv$S_yL2) / s2
      statLRb[k] <- (y[T2]^2 - rho_a * rv$S_yL2) / s2
    }

    ind_minEM <- which.min(tstat)
    T2_EM <- ind_minEM + T1 + Teps - 1
    ind_minLRa <- which.min(statLRa)
    T2_LR <- ind_minLRa + T1 + Teps - 1
    ind_minLRb <- which.min(statLRb)
    T2_LRb <- ind_minLRb + T1 + Teps - 1

    EMb_stat <- min(tstat) /
      sqrt(abs(rho_a) * phi_a^(2 * (T2_EM - T1)) * N / 2)
    EMa_stat <- sum(tstat) /
      sqrt(abs(phi_a^(2 * (N - T1 + 1)) * N / 2 / rho_a))
    LRa_stat <- min(statLRa) / (N * phi_a^(2 * (T2_LR - T1)) / 2)
    LRb_stat <- min(statLRb) / (N * phi_a^(2 * (T2_LRb - T1)) / 2)

    cval_EM12 <- get.cv.emergence(T1, T2_EM, N, (T1 - 1) / N)
    cval_LRa <- get.cv.emergence(T1, T2_LR, N)
    cval_LRb <- get.cv.emergence(T1, T2_LRb, N)

    cset_EMa12[T1 - 1] <- as.integer(EMa_stat > cval_EM12)
    cset_EMb12[T1 - 1] <- as.integer(EMb_stat > cval_EM12)
    cset_LRa12[T1 - 1] <- as.integer(LRa_stat > cval_LRa)
    cset_LRb12[T1 - 1] <- as.integer(LRb_stat > cval_LRb)
  }

  # --- T2 < T1 ---
  for (T1 in (Teps + 2):(N + 1)) {
    T2_range <- 2:(T1 - Teps)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      T2 <- T2_range[k]
      rv <- .nbcn.t.stat(y, T1, T2)
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy - rho_a * rv$S_yL2) / s2
    }

    EMa_stat <- sum(tstat) / N
    EMb_stat <- max(tstat)
    LR_stat <- max(statLR) / (N^2 * rho_a)

    cval_EM21 <- get.cv.emergence(T1, T1 - Teps, N, (T1 - 1) / N)
    cval_LR21 <- get.cv.emergence(T1, T1 - Teps, N)

    cset_EMa21[T1 - 1] <- as.integer(EMa_stat < cval_EM21["EMa21"])
    cset_EMb21[T1 - 1] <- as.integer(EMb_stat < cval_EM21["EMb21"])
    cset_LRa21[T1 - 1] <- as.integer(LR_stat < cval_LR21)
  }

  list(
    cset_EMa = pmin(cset_EMa12, cset_EMa21),
    cset_EMb = pmin(cset_EMb12, cset_EMb21),
    cset_LRa = pmin(cset_LRa12, cset_LRa21),
    cset_LE = pmin(cset_LRb12, cset_EMa21),
    cset_EMa12 = cset_EMa12,
    cset_EMb12 = cset_EMb12,
    cset_LRa12 = cset_LRa12,
    cset_LRb12 = cset_LRb12,
    cset_LE12 = cset_LRb12,
    cset_EMa21 = cset_EMa21,
    cset_EMb21 = cset_EMb21,
    cset_LRa21 = cset_LRa21,
    cset_LE21 = cset_EMa21
  )
}


#' @rdname bubbles.nbcn
#' @order 2
#' @export
cs.bubble.collapse <- function(
  y,
  lambda_e,
  phi_a,
  phi_b,
  s2,
  trim = 0.1
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- length(y) - 1
  Teps <- max(1L, floor(N * trim))

  rho_a <- phi_a - 1
  rho_b <- 1 - phi_b

  # T1 < T2
  cval_EM_lt <- get.cv.collapse(lambda_e, tail = "right")
  cval_LR_lt <- get.cv.collapse(lambda_e, tail = "right", LR = TRUE)

  # T2 < T1
  cval_EM_gt <- get.cv.collapse(lambda_e, tail = "left")
  cval_LR_gt <- get.cv.collapse(lambda_e, tail = "left", LR = TRUE)

  cset_EMa12 <- integer(N)
  cset_EMb12 <- integer(N)
  cset_LRa12 <- integer(N)

  cset_EMa21 <- integer(N)
  cset_EMb21 <- integer(N)
  cset_LRa21 <- integer(N)

  # --- T1 < T2 ---
  for (T1 in 2:(N + 1 - Teps)) {
    T2_range <- (T1 + Teps):(N + 1)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      rv <- .nbcn.t.stat(y, T1, T2_range[k])
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + (2 - phi_a - phi_b) * rv$S_yL2) / s2
    }

    scale_EM <- sqrt(N * phi_a^(2 * (T1 - 1)) * abs(rho_b) / 2)
    scale_LR <- N * (phi_a - phi_b) * phi_a^(2 * (T1 - 1)) / (2 * rho_b)

    cset_EMa12[T1 - 1] <- as.integer(mean(tstat) / scale_EM < cval_EM_lt)
    cset_EMb12[T1 - 1] <- as.integer(max(tstat) / scale_EM < cval_EM_lt)
    cset_LRa12[T1 - 1] <- as.integer(max(statLR) / scale_LR < cval_LR_lt)
  }

  # --- T2 < T1 ---
  for (T1 in (Teps + 2):(N + 1)) {
    T2_range <- 2:(T1 - Teps)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      rv <- .nbcn.t.stat(y, T1, T2_range[k])
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + (2 - phi_a - phi_b) * rv$S_yL2) / s2
    }

    scale_EM <- sqrt(N * phi_a^(2 * (T1 - 1)) * abs(rho_a) / 2)
    scale_LR <- N * (phi_a - phi_b) * phi_a^(2 * (T1 - 1)) / (2 * rho_a)

    idx <- T1 - 1
    cset_EMa21[idx] <- as.integer(mean(tstat) / scale_EM > cval_EM_gt)
    cset_EMb21[idx] <- as.integer(min(tstat) / scale_EM > cval_EM_gt)
    cset_LRa21[idx] <- as.integer(min(statLR) / scale_LR > cval_LR_gt)
  }

  list(
    cset_EMa = pmin(cset_EMa12, cset_EMa21),
    cset_EMb = pmin(cset_EMb12, cset_EMb21),
    cset_LRa = pmin(cset_LRa12, cset_LRa21),
    cset_EMa12 = cset_EMa12,
    cset_EMb12 = cset_EMb12,
    cset_LRa12 = cset_LRa12,
    cset_EMa21 = cset_EMa21,
    cset_EMb21 = cset_EMb21,
    cset_LRa21 = cset_LRa21
  )
}

#' @rdname bubbles.nbcn
#' @order 3
#' @export
cs.bubble.recovery <- function(
  y,
  lambda_e,
  lambda_c,
  phi_a,
  phi_b,
  s2,
  trim = 0.1
) {
  y <- as.matrix(y)

  N <- length(y) - 1
  Teps <- max(1L, floor(N * trim))
  rho_b <- 1 - phi_b

  exp_ac <- floor(N * (lambda_c - lambda_e))

  cset_LRa12 <- integer(N)
  cset_EMa12s <- integer(N)
  cset_EMb12s <- integer(N)
  cset_LRa21 <- integer(N)
  cset_LRb21 <- integer(N)
  cset_EMa21s <- integer(N)
  cset_EMb21s <- integer(N)

  # --- T1 < T2 ---
  for (T1 in 2:(N + 1 - Teps)) {
    T2_range <- (T1 + Teps):(N + 1)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      rv <- .nbcn.t.stat(y, T1, T2_range[k])
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + rho_b * rv$S_yL2) / s2
    }

    ind_minLR <- which.min(statLR)
    T2_LR <- ind_minLR + T1 + Teps - 1
    scale_LRa <- N *
      (T2_LR - T1) *
      rho_b *
      phi_a^(2 * exp_ac) *
      phi_b^(2 * (T1 - 1))

    LRa_stat <- min(statLR) / scale_LRa
    EMas_stat <- sum(tstat) / N
    EMbs_stat <- min(tstat)

    cv_EM <- get.cv.recovery(lambda_e, (T1 - 1) / N, tail = "right")
    cv_LR <- get.cv.recovery(lambda_e, (T1 - 1) / N, tail = "right", LR = TRUE)

    cset_LRa12[T1 - 1] <- as.integer(LRa_stat > cv_LR["LRas"])
    cset_EMa12s[T1 - 1] <- as.integer(EMas_stat > cv_EM["EMas"])
    cset_EMb12s[T1 - 1] <- as.integer(EMbs_stat > cv_EM["EMbs"])
  }

  # --- T2 < T1 ---
  for (T1 in (2 * Teps + 1):(N + 1)) {
    T2_range <- (Teps + 1):(T1 - Teps)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))
    statLRb <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      T2 <- T2_range[k]
      rv <- .nbcn.t.stat(y, T1, T2)
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + rho_b * rv$S_yL2) / s2
      statLRb[k] <- (-y[T2 - 1]^2 - rv$S_dy2 + rho_b * rv$S_yL2) / s2
    }

    ind_maxLR <- which.max(statLR)
    T2_LR <- ind_maxLR + Teps
    ind_maxLRb <- which.max(statLRb)
    T2_LRb <- ind_maxLRb + Teps
    ind_maxEM <- which.max(tstat)
    T2_EM <- ind_maxEM + Teps

    scale_LRa <- phi_a^(2 * exp_ac) * phi_b^(2 * (T2_LR - 1)) * N / 2
    scale_LRb <- phi_a^(2 * exp_ac) * phi_b^(2 * (T2_LRb - 1)) * N / 2
    scale_EMa <- phi_a^(exp_ac) * phi_b^(Teps) * sqrt(N / abs(rho_b) / 2)
    scale_EMb <- phi_a^(exp_ac) *
      phi_b^(T2_EM - 1) *
      sqrt(N * abs(rho_b) / 2)

    LRa_stat <- max(statLR) / scale_LRa
    LRb_stat <- max(statLRb) / scale_LRb
    EMas_stat <- sum(tstat) / scale_EMa
    EMbs_stat <- max(tstat) / scale_EMb

    cv_EM <- get.cv.recovery(lambda_e, lambda1 = 0, tail = "left")
    cv_LR <- get.cv.recovery(lambda_e, lambda1 = 0, tail = "left", LR = TRUE)

    cset_EMa21s[T1 - 1] <- as.integer(EMas_stat < cv_EM)
    cset_EMb21s[T1 - 1] <- as.integer(EMbs_stat < cv_EM)
    cset_LRa21[T1 - 1] <- as.integer(LRa_stat < cv_LR)
    cset_LRb21[T1 - 1] <- as.integer(LRb_stat < cv_LR)
  }

  list(
    cset_EMas = pmin(cset_EMa12s, cset_EMa21s),
    cset_EMbs = pmin(cset_EMb12s, cset_EMb21s),
    cset_LRa = pmin(cset_LRa12, cset_LRa21),
    cset_LE = pmin(cset_EMb12s, cset_LRb21),
    cset_EMa12s = cset_EMa12s,
    cset_EMb12s = cset_EMb12s,
    cset_LRa12 = cset_LRa12,
    cset_LE12 = cset_EMb12s,
    cset_EMa21s = cset_EMa21s,
    cset_EMb21s = cset_EMb21s,
    cset_LRa21 = cset_LRa21,
    cset_LRb21 = cset_LRb21,
    cset_LE21 = cset_LRb21
  )
}


#' @rdname bubbles.nbcn
#' @order 4
#' @details
#' [segments.AR1] function fits
#' \deqn{y_t = \phi_1 y_{t-1} I(t \leq bp) + \phi_2 y_{t-1} I(t > bp) + e_t}
#' over all candidate breakpoints in \eqn{[N \times trim, N \times (1 - trim)]}.
#' @return [segments.AR1] returns a named list of:
#' * break.point: index of estimates breakpoint,
#' * coefficients: AR(1) coefficients for pre-bubble and exploding regimes,
#' * SSR: minimum value of SSR,
#' * s.sq: estimated variance of internal model resuduals \eqn{\hat{e}_t}.
segments.AR1 <- function(
  y,
  trim = 0.1
) {
  y <- as.matrix(y)
  yL <- na.omit(.lagn(y, 1))
  y <- y[-1]

  N <- length(y)

  SSR <- Inf
  resBreak <- N
  resBeta <- NULL
  resSsq <- NULL

  tr <- .trend(N)
  bp_min <- max(1, floor(N * trim))
  bp_max <- N - bp_min

  for (bp in bp_min:bp_max) {
    mX <- cbind(
      yL * (tr <= bp),
      yL * (tr > bp)
    )

    loopModel <- OLS.reg(y, mX)
    loopSSR <- sum(loopModel$residuals^2)

    if (loopSSR < SSR) {
      resBreak <- bp
      SSR <- loopSSR
      resBeta <- loopModel$coefficients
      resSsq <- loopModel$s.sq
    }
  }

  list(
    SSR = SSR,
    s.sq = resSsq,
    break.point = resBreak,
    coefficients = resBeta
  )
}

#' @rdname bubbles.nbcn
#' @order 5
#' @details
#' [segments.NBCN] function uses a sequential three-step search:
#'   1. Find Tc on the full series.
#'   2. Find Te on y\[1:(Tc+1)\].
#'   3. Find Tr on y\[(Tc+2):end\].
#' Then refits the 4-regime model to get `phi_a`, `phi_b`, and `s.sq`.
#' @return [segments.NBCN] returns a named list of:
#' * Te_est, Tc_est, Tr_est: estimates moments of bubble emerging, break, and post-break restoration,
#' * phi_a, phi_b: AR(1) coefficients for exploding and post-break regimes,
#' * s.sq: estimated variance of internal model resuduals \eqn{\hat{e}_t}.
segments.NBCN <- function(y, trim = 0.1) {
  y <- as.matrix(y)

  N <- length(y)

  Tc_est <- segments.AR1(y, trim)$break.point
  Te_est <- segments.AR1(y[1:(Tc_est + 1)], trim)$break.point
  Tr_est <- (Tc_est + 1) + segments.AR1(y[(Tc_est + 2):N], trim)$break.point

  yL <- na.omit(.lagn(y, 1))
  y <- y[-1]
  N <- N - 1
  tr <- .trend(N)

  mX <- cbind(
    yL * (tr <= Te_est),
    yL * (tr > Te_est) * (tr <= Tc_est),
    yL * (tr > Tc_est) * (tr <= Tr_est),
    yL * (tr > Tr_est)
  )

  tmpModel <- OLS.reg(y, mX)

  list(
    Te_est = Te_est,
    Tc_est = Tc_est,
    Tr_est = Tr_est,
    phi_a = tmpModel$coefficients[2],
    phi_b = tmpModel$coefficients[3],
    s.sq = tmpModel$s.sq
  )
}


#' Plotting break dates for bubbles
#'
#' @param y Time series of interest.
#' @param trim Trimming parameter for breakpoints estimation procedure.
#' @param break_points Breakpoints to be plotted;
#' if NULL (default) then these points are estimated.
#' If `date_first` is not NULL `break_points` should bu numeric or Date,
#' otherwise it should be numeric only.
#' @param plot_type Select the style of the resulting graph.
#' @param date_first If not NULL then the X-axis will be labeled with dates.
#' @param date_by Time delta of your data.
#' @param date_format Format of dates on the graph.
#' @param y_label An optional label for Y-axis.
#' @param title An optional title of the graph.
#'
#' @import ggplot2
#' @import patchwork
#'
#' @references
#' Kurozumi, Eiji, and Anton Skrobotov. 2026.
#' "Confidence Sets for the Emergence, Collapse, and Recovery Dates of a Bubble".
#' arXiv:2511.16172.
#' Preprint, arXiv. https://doi.org/10.48550/arXiv.2511.16172.
#'
#' @export
plot_bubble <- function(
  y,
  trim = 0.1,
  break_points = NULL,
  plot_type = c("paper", "presentation"),
  date_first = NULL,
  date_format = "%Y-%m",
  date_by = "month",
  y_label = NULL,
  title = NULL
) {
  plot_type <- match.arg(plot_type)
  params <- plot_bubbles_style(plot_type)

  N <- length(y)

  # X-axis values
  x_data <- if (!is.null(date_first)) {
    seq(as.Date(date_first), by = date_by, length.out = N)
  } else {
    1:N
  }

  # Break dates
  if (is.null(break_points)) {
    segments <- segments.NBCN(y, trim)
    break_points <- c(
      segments$Te_est,
      segments$Tc_est,
      segments$Tr_est
    )
  }

  # X-axis values
  if (!is.null(date_first)) {
    if (is.numeric(break_points)) {
      break_points <- x_data[break_points]
    }
    stopifnot(inherits(break_points, "Date"))
    vline_lbl <- format(break_points, date_format)
  } else {
    stopifnot(is.numeric(break_points))
    vline_lbl <- break_points
  }

  # Vertical lines label positioning
  y_range <- range(y, na.rm = TRUE)
  vline_y <- y_range[1] + 0.10 * diff(y_range)
  x_mid <- mean(range(x_data))
  label_hjust <- ifelse(break_points < x_mid, 1.1, -0.1)

  # Plot construction
  p <- ggplot(mapping = aes(x = x_data, y = y)) +
    # X-axis
    (if (!is.null(date_first)) {
      scale_x_date(
        date_labels = date_format,
        expand = expansion(mult = 0.01)
      )
    } else {
      scale_x_continuous(
        expand = expansion(mult = 0.01)
      )
    }) +
    # Y-axis
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    # Main line
    geom_line(
      colour = params$line_colour,
      linewidth = params$line_size
    ) +
    # Vertical lines
    geom_vline(
      xintercept = as.numeric(break_points),
      colour = params$vline_colour,
      linetype = params$vline_type,
      linewidth = params$vline_size
    ) +
    # Vertical lines labels
    geom_text(
      mapping = aes(
        x = break_points,
        y = vline_y,
        label = vline_lbl,
        hjust = label_hjust
      ),
      angle = params$label_angle,
      size = params$label_size,
      colour = params$label_colour,
      inherit.aes = FALSE
    ) +

    labs(
      x = NULL,
      y = y_label,
      title = title
    ) +
    params$base_theme() +
    theme(
      panel.background = element_rect(fill = params$panel_bg, colour = NA),
      panel.grid.major = element_line(colour = params$grid_major),
      panel.grid.minor = element_line(colour = params$grid_minor),
      axis.text = element_text(size = params$axis_text_sz),
      axis.title.y = element_text(
        size = params$axis_title_sz,
        margin = margin(r = 6)
      ),
      plot.title = element_text(
        size = params$title_sz,
        face = "bold",
        margin = margin(b = 8)
      ),
      plot.margin = margin(10, 15, 8, 8)
    )

  suppressWarnings(print(p))
}


plot_bubbles_style <- function(plot_type) {
  switch(
    plot_type,

    paper = list(
      # цвета
      line_colour = "black",
      vline_colour = "#CC0000",
      label_colour = "grey30",
      # размеры
      line_size = 0.55,
      vline_size = 0.55,
      vline_type = "dashed",
      label_size = 2.8,
      label_angle = 90,
      # оси
      axis_text_sz = 9,
      axis_title_sz = 10,
      title_sz = 11,
      # тема
      base_theme = theme_bw,
      panel_bg = "white",
      grid_major = "grey88",
      grid_minor = "grey94"
    ),

    presentation = list(
      # цвета
      line_colour = "#1a1a1a",
      vline_colour = "#E8000D",
      label_colour = "#333333",
      # размеры
      line_size = 0.85,
      vline_size = 0.85,
      vline_type = "dashed",
      label_size = 4.0,
      label_angle = 90,
      # оси
      axis_text_sz = 13,
      axis_title_sz = 14,
      title_sz = 15,
      # тема
      base_theme = theme_gray,
      panel_bg = "#F0F0F0",
      grid_major = "white",
      grid_minor = "white"
    )
  )
}
