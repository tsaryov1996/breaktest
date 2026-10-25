print.bt_mdfCHLT <- function(obj, ...) {
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
