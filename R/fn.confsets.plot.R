# Shared preparation for confidence-set plots. Membership is never interpolated.
.confsets_restore_time <- function(x, template) {
  if (inherits(template, "Date")) return(as.Date(x, origin = "1970-01-01"))
  if (inherits(template, "POSIXct")) {
    tz <- attr(template, "tzone")
    if (is.null(tz) || !length(tz)) tz <- ""
    return(as.POSIXct(x, origin = "1970-01-01", tz = tz[1L]))
  }
  as.numeric(x)
}

.confsets_prepare <- function(y, conf_sets, break_points, set_index, time,
                              date_first, date_by, bubble = FALSE) {
  if (!is.numeric(y) || (is.matrix(y) && ncol(y) != 1L) ||
      length(dim(y)) > 2L || length(y) < 2L ||
      any(is.infinite(y)) || sum(is.finite(y)) < 2L) {
    stop("'y' must be a numeric vector or one-column matrix with at least two finite values.")
  }
  N <- length(y)
  if (!is.null(time) && !is.null(date_first)) {
    stop("Supply either 'time' or 'date_first', not both.")
  }
  if (is.null(time)) {
    time <- if (!is.null(date_first)) {
      if (length(date_first) != 1L || is.na(as.Date(date_first))) {
        stop("'date_first' must specify one valid date.")
      }
      seq(as.Date(date_first), by = date_by, length.out = N)
    } else if (stats::is.ts(y)) {
      as.numeric(stats::time(y))
    } else seq_len(N)
  }
  if (!(is.numeric(time) || inherits(time, "Date") || inherits(time, "POSIXct")) ||
      !is.null(dim(time)) || length(time) != N ||
      any(!is.finite(as.numeric(time))) || any(diff(as.numeric(time)) <= 0)) {
    stop("'time' must contain one finite, strictly increasing numeric, Date or POSIXct value per observation.")
  }
  if (!is.list(conf_sets) || !length(conf_sets) || is.null(names(conf_sets)) ||
      anyNA(names(conf_sets)) || any(!nzchar(names(conf_sets))) ||
      anyDuplicated(names(conf_sets))) {
    stop("'conf_sets' must be a nonempty list with unique, nonempty names.")
  }
  if (bubble) {
    if (!setequal(names(conf_sets), c("Te", "Tc", "Tr"))) {
      stop("For a bubble, name the three confidence sets 'Te', 'Tc' and 'Tr'.")
    }
    conf_sets <- conf_sets[c("Te", "Tc", "Tr")]
  }
  set_names <- names(conf_sets)
  if (is.list(set_index)) {
    if (is.null(names(set_index)) || anyNA(names(set_index)) ||
        anyDuplicated(names(set_index)) || !setequal(names(set_index), set_names)) {
      stop("A list 'set_index' must have the same unique names as 'conf_sets'.")
    }
    set_index <- set_index[set_names]
  } else {
    set_index <- rep(list(set_index), length(conf_sets))
  }
  xx <- as.numeric(time)
  middle <- xx[-N] + diff(xx) / 2
  left <- c(xx[1L] - diff(xx)[1L] / 2, middle)
  right <- c(middle, xx[N] + diff(xx)[N - 1L] / 2)
  if (any(!is.finite(c(left, right)))) stop("The range of 'time' is too large.")
  components <- data.frame(
    set = character(), start = integer(), end = integer(), n = integer(),
    xmin = time[FALSE], xmax = time[FALSE],
    xstart = time[FALSE], xend = time[FALSE]
  )
  memberships <- data.frame(set = character(), index = integer(), accepted = logical())
  for (i in seq_along(conf_sets)) {
    z <- conf_sets[[i]]
    if (!(is.numeric(z) || is.logical(z)) ||
        (is.matrix(z) && ncol(z) != 1L) || length(dim(z)) > 2L ||
        anyNA(z) || any(!z %in% c(0, 1))) {
      stop("Each confidence set must be a 0/1 or logical vector (or one-column matrix), without NA.")
    }
    z <- as.logical(z)
    index <- set_index[[i]]
    if (is.null(index)) {
      if (length(z) != N) {
        stop("A shorter confidence set needs explicit 'set_index'; window offsets are never guessed.")
      }
      index <- seq_len(N)
    }
    if (!is.numeric(index) || !is.null(dim(index)) || length(index) != length(z) ||
        any(!is.finite(index)) || any(index != floor(index)) ||
        any(index < 1L | index > N) || any(diff(index) <= 0)) {
      stop("'set_index' must contain strictly increasing, unique row indices in 1:length(y), matching each set.")
    }
    index <- as.integer(index)
    memberships <- rbind(memberships, data.frame(
      set = rep(set_names[i], length(z)), index = index, accepted = z
    ))
    accepted <- index[z]
    if (!length(accepted)) next
    groups <- split(accepted, cumsum(c(TRUE, diff(accepted) != 1L)))
    for (g in groups) {
      first <- g[1L]
      last <- g[length(g)]
      components <- rbind(components, data.frame(
        set = set_names[i], start = first, end = last, n = length(g),
        xmin = .confsets_restore_time(left[first], time),
        xmax = .confsets_restore_time(right[last], time),
        xstart = time[first], xend = time[last]
      ))
    }
  }
  estimates <- data.frame(set = character(), index = integer(), x = time[FALSE])
  if (!is.null(break_points)) {
    if (is.logical(break_points) && all(is.na(break_points))) {
      break_points <- stats::setNames(as.numeric(break_points), names(break_points))
    }
    if (length(break_points) == 1L && !bubble) {
      break_points <- rep(unname(break_points), length(set_names))
    } else if (!is.null(names(break_points))) {
      if (anyNA(names(break_points)) || anyDuplicated(names(break_points)) ||
          !setequal(names(break_points), set_names)) {
        stop("Named 'break_points' must match the names of 'conf_sets'.")
      }
      break_points <- break_points[set_names]
    }
    if (length(break_points) != length(set_names)) {
      stop("Supply one 'break_points' value per set, or one common value for 'plot_confsets'.")
    }
    if (inherits(break_points, "Date") || inherits(break_points, "POSIXct")) {
      if (!identical(class(break_points), class(time))) {
        stop("Calendar 'break_points' and 'time' must have the same class.")
      }
      index <- match(break_points, time)
      if (any(!is.na(break_points) & is.na(index))) {
        stop("Calendar 'break_points' must match observed times exactly.")
      }
    } else {
      if (!is.numeric(break_points)) stop("'break_points' must be row indices or calendar dates.")
      index <- break_points
      valid <- !is.na(index)
      if (any(!is.finite(index[valid])) || any(index[valid] != floor(index[valid])) ||
          any(index[valid] < 1L | index[valid] > N)) {
        stop("Numeric 'break_points' are row indices in 1:length(y), not time coordinates.")
      }
    }
    keep <- !is.na(index)
    estimates <- data.frame(set = set_names[keep], index = as.integer(index[keep]),
                            x = time[index[keep]])
    if (bubble && any(diff(estimates$index) <= 0L)) {
      stop("Bubble estimates must be ordered Te < Tc < Tr.")
    }
  }
  list(time = time, y = as.numeric(y), set_names = set_names,
       memberships = memberships, components = components, estimates = estimates,
       limits = .confsets_restore_time(c(left[1L], right[N]), time),
       phases = data.frame(set = character(), xmin = time[FALSE], xmax = time[FALSE]))
}

.confsets_alpha <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value > 1) stop(paste0("'", name, "' must be between zero and one."))
}

.confsets_colours <- function(colours, set_names) {
  if (!is.character(colours) || !length(colours) || anyNA(colours)) {
    stop("'colours' must contain valid R colour strings.")
  }
  if (length(colours) == 1L) colours <- rep(unname(colours), length(set_names))
  if (!is.null(names(colours))) {
    if (anyDuplicated(names(colours)) || !all(set_names %in% names(colours))) {
      stop("Named 'colours' must cover every confidence set.")
    }
    colours <- colours[set_names]
  }
  if (length(colours) != length(set_names)) stop("Supply one colour per confidence set, or a single colour.")
  tryCatch(grDevices::col2rgb(colours), error = function(e) stop("Invalid value in 'colours'."))
  stats::setNames(colours, set_names)
}

.confsets_xscale <- function(d, date_format) {
  if (!is.character(date_format) || length(date_format) != 1L || is.na(date_format)) {
    stop("'date_format' must be one format string.")
  }
  if (inherits(d$time, "Date")) {
    ggplot2::scale_x_date(limits = d$limits, date_labels = date_format, expand = c(0, 0))
  } else if (inherits(d$time, "POSIXct")) {
    ggplot2::scale_x_datetime(limits = d$limits, date_labels = date_format, expand = c(0, 0))
  } else ggplot2::scale_x_continuous(limits = d$limits, expand = c(0, 0))
}

.confsets_theme <- function(plot_type) {
  # Match the existing plot_bubble() paper/presentation appearance.
  size <- if (plot_type == "paper") 10 else 14
  base <- if (plot_type == "paper") ggplot2::theme_bw(base_size = size) else
    ggplot2::theme_gray(base_size = size)
  base + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(10, 15, 8, 8),
    legend.position = "bottom"
  )
}

.confsets_series <- function(d, plot_type, date_format, y_label, title) {
  x <- y <- NULL
  ggplot2::ggplot(data.frame(x = d$time, y = d$y), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(linewidth = if (plot_type == "paper") 0.55 else 0.85,
                       colour = "black", na.rm = TRUE) +
    .confsets_xscale(d, date_format) + .confsets_theme(plot_type) +
    ggplot2::labs(x = NULL, y = y_label, title = title)
}
