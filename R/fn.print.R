#' @title
#' Custom functions for printing results in a nice way.
#'
#' @param x An object containing results.
#' @param ... Any additional arguments for [print] function.
#'
#' @keywords internal
#' @exportS3Method
print.bt_SADF <- function(x, ...) {
  test_name <- NULL

  if ("SADF.value" %in% names(x)) {
    test_statistic <- x$SADF.value
    test_name <- "SADF"
  } else if ("STADF.value" %in% names(x)) {
    test_statistic <- x$STADF.value
    test_name <- "STADF"
  } else if ("GSADF.value" %in% names(x)) {
    test_statistic <- x$GSADF.value
    test_name <- "GSADF"
  } else if ("GSTADF.value" %in% names(x)) {
    test_statistic <- x$GSTADF.value
    test_name <- "GSTADF"
  }

  if (!is.null(test_name)) {
    cat(test_name, "test statistic:", test_statistic, "\n")
  }

  if ("supBZ.value" %in% names(x)) {
    cat("supBZ statistic:", x$supBZ.value, "\n")
  }

  if ("p.value" %in% names(x)) {
    cat("P-value:", x$p.value, "\n")
  }

  if ("is.explosive" %in% names(x)) {
    cat(
      "Current process is",
      ifelse(x$is.explosive == 1, "", "not"),
      "explosive\n"
    )
  }
}


#' @rdname print.bt_SADF
#' @importFrom stringr str_split
#' @keywords internal
#' @exportS3Method
print.bt_mdfHLT <- function(obj, ...) {
  if (obj$const && !obj$trend) {
    cat("Model 0: Structural change in intercept\n")
    cat(
      "Y{t}=a0+a1*DU+b0*t+e{t}\n",
      "    where DU=1(t>TB)\n"
    )
  } else if (!obj$const && obj$trend) {
    cat("Model 1: Structural change in slope\n")
    cat(
      "Y{t}=a0+b0*t+b1*DT+e{t}\n",
      "    where DT=1(t>TB)*(t-TB)\n"
    )
  } else if (obj$const && obj$trend) {
    cat("Model 2: Structural change in both intercept and slope\n")
    cat(
      "Y{t}=a0+a1*DU+b0*t+b1*DT+e{t}\n",
      "    where DU=1(t>TB) and DT=1(t>TB)*(t-TB)\n"
    )
  }
  cat("\n")

  cat("The break date is estimated in", obj$break.time, "\n")
  cat("\n")

  cat(
    "Robust tests for the break with uncertaint over errors",
    "(integrated or stationary):\n\n",
    "\tstat\tc.v.\n"
  )
  for (v in c("HLT", "PY")) {
    cat(sprintf("%-7s\t%.4f\t%.4f\n", v, obj[[v]]$statistic, obj[[v]]$cv))
    if (obj[[v]]$statistic > obj[[v]]$cv) {
      cat("reject\n")
    } else {
      cat("fails to reject\n")
    }
    cat("\n")
  }
  cat("\n")

  cat(
    "Unit root tests:\n\n",
    "\tstat\tc.v.\n"
  )
  for (v in c("DF.GLS", "DF.OLS", "MDF.GLS", "MDF.OLS")) {
    cat(sprintf("%-7s\t%.4f\t%.4f\n", v, obj[[v]]$statistic, obj[[v]]$cv))
    if (obj[[v]]$statistic < obj[[v]]$cv) {
      cat("reject\n")
    } else {
      cat("fails to reject\n")
    }
    cat("\n")
  }
  cat("\n")

  cat("Testing strategies:\n\n")
  for (v in c("A.HLT", "A.PY", "UR.HLT", "UR.PY", "URR.HLT", "URR.PY")) {
    tmp_str <- stringr::str_split(v, "\\.")
    cat(tmp_str[[1]][1], "*(t_", tmp_str[[1]][2], ", s_alpha): ", sep = "")
    if (obj[[v]] == 1) {
      cat("reject\n")
    } else {
      cat("fails to reject\n")
    }
    cat("\n")
  }
}


#' @rdname print.bt_SADF
#' @importFrom stringr str_split
#' @keywords internal
#' @exportS3Method
print.bt_mdfHLTN <- function(obj, ...) {
  cat("\t\tstat\tc.v.\n\n")

  for (v in c(
    "MDF.GLS.1",
    "MDF.GLS.2",
    if (obj$breaks == 3) "MDF.GLS.3" else NULL,
    "MDF.OLS.1",
    "MDF.OLS.2",
    if (obj$breaks == 3) "MDF.OLS.3" else NULL
  )) {
    cat(sprintf("%-9s:\t%.4f\t%.4f\n", v, obj[[v]]$statistic, obj[[v]]$cv))
    if (obj[[v]]$statistic < obj[[v]]$cv) {
      cat("reject\n")
    } else {
      cat("fails to reject\n")
    }
    cat("\n")
  }
  cat("\n")

  cat(sprintf("UR^%d(s.alpha): ", obj$breaks))
  if (obj$UR1 == 1) {
    cat("reject\n")
  } else {
    cat("fails to reject\n")
  }

  cat(sprintf("UR^%d(s.alpha, %d): ", obj$breaks, obj$breaks.star))
  if (obj$UR == 1) {
    cat("reject\n")
  } else {
    cat("fails to reject\n")
  }
  cat("\n")
}


#' @rdname print.bt_SADF
#' @keywords internal
#' @exportS3Method
print.mdfCHLT <- function(obj, ...) {
  cat("\t\tstat\tc.v.\t wild c.v.\n\n")
  for (v in c("MZa", "MSB", "MZt", "ADF")) {
    cat(
      sprintf(
        "%s stat:\t%.4f\t%.4f\t%.4f\n",
        v,
        obj[[v]]$statistic,
        obj[[v]]$cv,
        obj[[v]]$cv.bootstrap
      )
    )
  }
}


#' @rdname print.bt_SADF
#' @keywords internal
#' @exportS3Method
print.cointGH <- function(x, ...) {
  cat("Gregory-Hansen tests\n")
  cat("\t\tstat\tc.v.\tasymptotic\n\n")
  for (v in c("Za", "Zt", "ADF")) {
    cat(
      sprintf(
        "%s stat:\t%.4f\t%.4f\t%.4f\n",
        v,
        x[[v]]$statistic,
        x[[v]]$cv,
        x[[v]]$asy.cv
      )
    )
  }
}


#' @rdname print.bt_SADF
#' @keywords internal
#' @exportS3Method
print.bt_robustUR <- function(obj, ...) {
  cat(sprintf("Estimated break moment: %d\n", obj$break.time))
  cat("\t\tstat\tc.v.\t\n\n")
  for (v in c(
    "HLT",
    "PY",
    "DF.GLS",
    "DF.OLS",
    "MDF.GLS",
    "MDF.OLS",
    "MDF.t"
  )) {
    cat(
      sprintf(
        "%s stat:\t%.4f\t%.4f\n",
        v,
        obj[[v]]$statistic,
        obj[[v]]$cv
      )
    )
  }
}


#' @rdname print.bt_SADF
#' @keywords internal
#' @exportS3Method
print.bt_robustURN <- function(obj, ...) {
  cat(sprintf("Estimated break moment: %d\n", obj$breaks.star))
  cat(sprintf(
    "Statistic F(%d,%d) = %.4f\n",
    obj$KP.sequential$breaks,
    obj$KP.sequential$breaks + 1,
    obj$KP.sequential$statistic
  ))
  cat("\t\tstat\tc.v.\n\n")
  for (v in c("MDF.GLS.1", "MDF.OLS.1", "MDF.GLS.2", "MDF.OLS.2")) {
    cat(
      sprintf(
        "%s stat:\t%.4f\t%.4f\n",
        v,
        obj[[v]]$statistic,
        obj[[v]]$cv
      )
    )
  }
}


#' @rdname print.bt_SADF
#' @keywords internal
#' @exportS3Method
print.bt_cointPR <- function(obj, ...) {
  mstr <- max(mapply(nchar, names(obj)))
  cat("Perron-Rodríguez testing procedure\n\n")
  cat(sprintf(paste0("%", mstr, "s \tstat.\tcr.value\n"), ""))
  for (v in names(obj)) {
    if (v == "lag") {
      cat(sprintf("\n\nInternal ADF lag length: %d", obj[[v]]))
    } else {
      cat(sprintf(
        paste0("%", mstr, "s:\t% 6.3f\t% 6.3f\n"),
        v,
        obj[[v]]$statistic,
        obj[[v]]$c.value
      ))
    }
  }
}

#' @rdname print.bt_SADF
#' @keywords internal
#' @exportS3Method
print.bt_confSet <- function(obj, ...) {
  cat(
    "Kurozumi-Skrobotov procedure to find confidence intervals",
    "for a structural break date",
    sep = "\n"
  )
  for (v in names(obj)) {
    if (v == "td") {
      cat(sprintf("Break date: %i\n\n", which(obj[[v]] == -1)))
    } else {
      Ts <- which(obj[[v]] == 1)
      if (!length(Ts) == 0) cat(sprintf("%8s: [%i, %i]\n", v, min(Ts), max(Ts)))
    }
  }
}

#' @rdname print.bt_SADF
#' @keywords internal
#' @exportS3Method
print.bt_mdfCHLT <- function(obj, ...) {
  if (identical(obj$params$method, "CHLT2014-MDF")) {
    cat("GLS minimum-DF test with", obj$params$breaks,
        "possible trend breaks (CHLT, January 2014)\n")
    cat("Minimizing break dates:", obj$params$break.points, "\n")
    cat(sprintf("ADF: %.6f\n", obj$ADF$statistic))
    if (!is.null(obj$ADF$boot.cv)) {
      cat(sprintf("Bootstrap 5%% critical value: %.6f; p-value: %.6f\n",
                  obj$ADF$boot.cv, obj$ADF$p.value))
    } else {
      cat("Use boot.cv = TRUE or bootstrap() for critical values and p-values.\n")
    }
    return(invisible(obj))
  }

  cat(
    "A modified DF test for a single break",
    "and possible heteroscedasticity",
    "by Cavaliere, Harvey, Leybourne, and Taylor (2011)",
    sep = "\n"
  )
  if (!"boot.cv" %in% names(obj[["ADF"]])) {
    cat("     stat.\tc.v\n")
    for (v in names(obj)) {
      if (v == "params") {
        next
      }
      cat(sprintf("%s: % 6.3f\t% 6.3f\n", v, obj[[v]]$statistic, obj[[v]]$cv))
    }
  } else {
    cat("     stat.\tc.v\tboot c.v.\n")
    for (v in names(obj)) {
      if (v == "params") {
        next
      }
      cat(sprintf(
        "%s: % 6.3f\t% 6.3f\t% 6.3f\n",
        v,
        obj[[v]]$statistic,
        obj[[v]]$cv,
        obj[[v]]$boot.cv
      ))
    }
  }
}
