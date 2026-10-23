#' Plotting confidence sets for bubble dates as shaded bands
#'
#' @description Plot a time series with separate shaded components of confidence
#' sets for the emergence, collapse and recovery dates of a bubble.
#'
#' @param y Numeric time series or a one-column matrix. Missing observations
#' are allowed and retain their original positions.
#' @param conf_sets Named list of binary or logical membership vectors, named
#' `Te`, `Tc` and `Tr`. Select one already computed method per date.
#' One-column membership matrices are also accepted; missing membership is not.
#' @param break_points Optional estimated dates. Numeric values are **one-based
#' row indices in y**, even when `time` is supplied. Alternatively supply Date
#' or POSIXct values matching `time` exactly. Names, if supplied, must be
#' `Te`, `Tc`, `Tr`; otherwise this order is used. Use NA to omit an estimate.
#' @param set_index Optional indices mapping membership entries to rows of `y`:
#' a common integer vector or a named list with the same names as `conf_sets`.
#' Required for sets calculated on shorter windows. Without it each vector must
#' have `length(y)` entries. Indices must be strictly increasing and unique.
#' @param time Optional strictly increasing numeric, Date or POSIXct vector of
#' length `length(y)`. Defaults to `time(y)` for a ts object or observation numbers.
#' @param plot_type Style of the resulting graph: `paper` or `presentation`.
#' @param date_first Optional first date, used with `date_by` to construct a
#' regular Date axis. Cannot be combined with `time`.
#' @param date_format Format of calendar labels.
#' @param date_by Time increment passed to `seq.Date`, such as `"month"`.
#' @param colours One colour or one colour per confidence set; names may be used.
#' @param conf_alpha Opacity of the confidence bands, between zero and one.
#' @param shade_phases Also show pale estimated phase spans between Te and Tc,
#' between Tc and Tr, and after Tr. Requires all three `break_points`.
#' @param phase_alpha Opacity of these optional phase spans.
#' @param y_label Optional Y-axis label.
#' @param title Optional plot title.
#'
#' @details
#' These functions visualize supplied sets; they do not estimate dates or compute
#' confidence levels. A disconnected set is never replaced by its convex hull.
#' For bands, each accepted observation occupies a cell bounded by midpoints
#' between adjacent observed times. The first and last cells extend by half of
#' the adjacent time gap. Thus singleton sets remain visible. Connectivity is
#' defined by consecutive rows of `y`, not by consecutive calendar days.
#'
#' Bands describe uncertainty about a **date**, not the duration of a phase.
#' Optional phase shading is determined by the point estimates and is labelled
#' separately in the caption. The three supplied sets need not constitute a
#' joint confidence region; this function does not infer simultaneous coverage.
#'
#' Index mapping is explicit. When a bubble function receives `y[a:b]` including
#' its initial lag observation, its length `b-a` membership vector is mapped with
#' `set_index = seq.int(a + 1L, b)`. For a full series including `y0`, estimates
#' from `segments.NBCN()` must be mapped to displayed rows by adding one.
#' Do not apply this offset to objects already indexed by rows of `y`.
#'
#' @return A ggplot object, displayed by R when evaluated at the console.
#' Inside loops or functions use `print()`. The `confsets_data` attribute contains
#' the time grid, membership table, connected components, mapped estimates and
#' optional phase spans, for inspecting the exact plotted data.
#'
#' @references
#' Kurozumi, E. and Skrobotov, A. (2026). Confidence Sets for the Emergence,
#' Collapse, and Recovery Dates of a Bubble. Econometric Reviews, 2026, 45, 1529-1563.
#'
#' @seealso [plot_confsets()], [plot_bubble()], [cs.bubble.emerge()]
#' @examples
#' t <- seq_len(100)
#' y <- sin(t / 8) + pmax(t - 30, 0) / 15 - pmax(t - 65, 0) / 8
#' sets <- list(Te = t %in% c(26:33, 38), Tc = t %in% 62:68,
#'              Tr = t %in% 83:91)
#' plot_bubble_confsets(y, sets, break_points = c(30, 65, 87))
#' @export
plot_bubble_confsets <- function(
  y,
  conf_sets,
  break_points = NULL,
  set_index = NULL,
  time = NULL,
  plot_type = c("paper", "presentation"),
  date_first = NULL,
  date_format = "%Y-%m",
  date_by = "month",
  colours = c(Te = "#F4B6C2", Tc = "#6DEB6B", Tr = "#8277F5"),
  conf_alpha = 0.4,
  shade_phases = FALSE,
  phase_alpha = 0.08,
  y_label = NULL,
  title = NULL
) {
  plot_type <- match.arg(plot_type)
  .confsets_alpha(conf_alpha, "conf_alpha")
  .confsets_alpha(phase_alpha, "phase_alpha")
  if (!is.logical(shade_phases) || length(shade_phases) != 1L || is.na(shade_phases)) {
    stop("'shade_phases' must be TRUE or FALSE.")
  }
  d <- .confsets_prepare(y, conf_sets, break_points, set_index, time,
                         date_first, date_by, bubble = TRUE)
  colours <- .confsets_colours(colours, d$set_names)
  if (shade_phases) {
    if (nrow(d$estimates) != 3L) stop("Phase shading requires all three estimated dates.")
    d$phases <- data.frame(set = d$set_names, xmin = d$estimates$x,
                          xmax = c(d$estimates$x[-1L], d$limits[2L]))
  }
  x <- y <- xmin <- xmax <- set <- NULL
  p <- ggplot2::ggplot()
  if (shade_phases) {
    p <- p + ggplot2::geom_rect(
      data = d$phases,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = set),
      alpha = phase_alpha, show.legend = FALSE
    )
  }
  p <- p + ggplot2::geom_rect(
    data = d$components,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = set),
    alpha = conf_alpha
  ) + ggplot2::geom_line(
    data = data.frame(x = d$time, y = d$y), ggplot2::aes(x = x, y = y),
    colour = "black", linewidth = if (plot_type == "paper") 0.55 else 0.85,
    na.rm = TRUE
  )
  yr <- range(d$y, na.rm = TRUE)
  yrange <- diff(yr)
  if (yrange == 0) yrange <- max(1, abs(yr[1L]) * 0.1)
  if (nrow(d$estimates)) {
    p <- p + ggplot2::geom_segment(
      data = d$estimates,
      ggplot2::aes(x = x, xend = x),
      y = yr[1L] - 0.04 * yrange, yend = yr[2L] + 0.08 * yrange,
      arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed"),
      linewidth = 0.55, colour = "black"
    ) + ggplot2::geom_text(
      data = d$estimates, ggplot2::aes(x = x, label = set),
      y = yr[2L] + 0.13 * yrange,
      size = if (plot_type == "paper") 3 else 4, colour = "black"
    ) + ggplot2::expand_limits(y = c(yr[1L] - 0.04 * yrange,
                                    yr[2L] + 0.20 * yrange))
  }
  empty <- setdiff(d$set_names, d$components$set)
  caption <- c(if (length(empty)) paste("Empty confidence set:", paste(empty, collapse = ", ")),
               if (shade_phases) "Pale background: estimated phase spans; darker bands: confidence sets.")
  p <- p + .confsets_xscale(d, date_format) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    ggplot2::scale_fill_manual(
      values = colours, limits = d$set_names, drop = FALSE,
      labels = c("Te: emergence", "Tc: collapse", "Tr: recovery"),
      name = "Confidence sets for dates"
    ) + .confsets_theme(plot_type) +
    ggplot2::labs(x = NULL, y = y_label, title = title,
                  caption = if (length(caption)) paste(caption, collapse = "\n") else NULL)
  d$mode <- "bands"
  attr(p, "confsets_data") <- d
  p
}
