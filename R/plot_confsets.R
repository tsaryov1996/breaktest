#' Plotting confidence sets in rows below a time series
#'
#' @description Display a time series and aligned horizontal confidence sets,
#' with one row per method or break date. Gaps and singleton sets are preserved.
#'
#' @inheritParams plot_bubble_confsets
#' @param conf_sets Named list of binary or logical membership vectors or
#' one-column matrices. Names are used as row labels, in the supplied order.
#' @param break_points Optional estimate per row, or one common estimate for all
#' rows. Numeric values are row indices in `y`; Date/POSIXct values must match
#' `time`. A named vector must match `conf_sets`; NA omits the corresponding marker.
#' @param colours One colour or one colour per row; names may be used.
#' @param relative_heights Relative heights of the time-series and set panels.
#' @param x_label Optional label for the shared time axis.
#' @param empty_label Label printed in a row with no accepted date.
#'
#' @details
#' Computation of confidence sets is separate from plotting. Accepted consecutive
#' observations are drawn as segments between their first and last times;
#' a singleton is drawn as a point. Empty sets retain their labelled rows.
#' Estimates, when supplied, are shown as black diamonds on the corresponding rows.
#' Both panels have identical time limits. See [plot_bubble_confsets()] for the
#' explicit mapping of membership vectors from shorter windows to rows of `y`.
#'
#' @return A patchwork object with two aligned ggplot panels. Use `print()` to
#' display it in a loop or function, and `ggplot2::ggsave()` to save it.
#' The `confsets_data` attribute contains the exact components and estimates plotted.
#'
#' @references
#' Eo, Y. and Morley, J. (2015). Likelihood-ratio-based confidence sets for the
#' timing of structural breaks. Quantitative Economics, 6, 463-497.
#' \doi{10.3982/QE186}.
#' Kurozumi, E. and Skrobotov, A. (2018). Confidence Sets for the Break Date in
#' Cointegrating Regressions. Oxford Bulletin of Economics and Statistics, 80,
#' 514-535. \doi{10.1111/obes.12223}.
#'
#' @seealso [plot_bubble_confsets()], [coint.conf.sets()], [fcst.conf.set.EoMorley()]
#' @examples
#' t <- seq_len(100)
#' y <- sin(t / 8)
#' sets <- list("Break 1" = t %in% 20:35,
#'              "Break 2" = t %in% c(60:68, 75:78, 85))
#' plot_confsets(y, sets, break_points = c(28, 64))
#' @export
plot_confsets <- function(
  y,
  conf_sets,
  break_points = NULL,
  set_index = NULL,
  time = NULL,
  plot_type = c("paper", "presentation"),
  date_first = NULL,
  date_format = "%Y-%m",
  date_by = "month",
  colours = "#BD2A2A",
  relative_heights = c(2, 1),
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  empty_label = "Empty set"
) {
  plot_type <- match.arg(plot_type)
  if (!is.numeric(relative_heights) || length(relative_heights) != 2L ||
      any(!is.finite(relative_heights)) || any(relative_heights <= 0)) {
    stop("'relative_heights' must contain two positive finite numbers.")
  }
  if (!is.character(empty_label) || length(empty_label) != 1L || is.na(empty_label)) {
    stop("'empty_label' must be one string.")
  }
  d <- .confsets_prepare(y, conf_sets, break_points, set_index, time,
                         date_first, date_by)
  colours <- .confsets_colours(colours, d$set_names)
  x <- row <- xstart <- xend <- set <- label <- NULL
  K <- length(d$set_names)
  comp <- d$components
  comp$row <- K + 1L - match(comp$set, d$set_names)
  estimates <- d$estimates
  estimates$row <- K + 1L - match(estimates$set, d$set_names)
  empty <- setdiff(d$set_names, comp$set)
  empty_data <- data.frame(
    x = rep(.confsets_restore_time(as.numeric(d$limits[1L]) +
                                     0.015 * diff(as.numeric(d$limits)), d$time), length(empty)),
    row = K + 1L - match(empty, d$set_names), label = rep(empty_label, length(empty))
  )
  endpoints <- rbind(data.frame(set = comp$set, row = comp$row, x = comp$xstart),
                     data.frame(set = comp$set, row = comp$row, x = comp$xend))
  top <- .confsets_series(d, plot_type, date_format, y_label, title) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank())
  bottom <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = comp[comp$n > 1L, ],
      ggplot2::aes(x = xstart, xend = xend, y = row, yend = row, colour = set),
      linewidth = if (plot_type == "paper") 0.8 else 1.1
    ) + ggplot2::geom_segment(
      data = endpoints, ggplot2::aes(x = x, xend = x, y = row - 0.10,
                                    yend = row + 0.10, colour = set), linewidth = 0.45
    ) + ggplot2::geom_point(
      data = comp[comp$n == 1L, ], ggplot2::aes(x = xstart, y = row, colour = set),
      size = 2
    ) + ggplot2::geom_point(
      data = estimates, ggplot2::aes(x = x, y = row), shape = 18, size = 2.2,
      colour = "black"
    ) + ggplot2::geom_text(
      data = empty_data, ggplot2::aes(x = x, y = row, label = label), hjust = 0,
      colour = "grey45", fontface = "italic", size = if (plot_type == "paper") 3 else 4
    ) + .confsets_xscale(d, date_format) +
    ggplot2::scale_y_continuous(breaks = seq_len(K), labels = rev(d$set_names),
                                limits = c(0.5, K + 0.5), expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = colours, limits = d$set_names, guide = "none") +
    .confsets_theme(plot_type) +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                    panel.grid = ggplot2::element_blank(),
                    panel.background = ggplot2::element_blank(),
                    axis.line.x = ggplot2::element_line(colour = "black", linewidth = 0.3),
                    axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(x = x_label, y = NULL,
                  caption = if (nrow(estimates)) "Black diamonds: point estimates." else NULL)
  p <- patchwork::wrap_plots(top, bottom, ncol = 1, heights = relative_heights)
  d$mode <- "rows"
  attr(p, "confsets_data") <- d
  p
}
