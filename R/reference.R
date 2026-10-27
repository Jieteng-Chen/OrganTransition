#' Canonical organ names
#'
#' The nine organ systems used by the reference networks, in the order of the
#' final manuscript analysis script.
#'
#' @return A character vector of organ names.
#' @export
ot_organs <- function() {
  .organ_levels
}

#' Load the frozen reference object
#'
#' Reads the age-stratified reference networks (young: age <= 44, older:
#' age >= 65) and the reference transition-score distribution that were frozen
#' when the package was built. No network fitting is performed at load time.
#'
#' @param path Optional path to a `reference.rds` file. If `NULL`, the copy
#'   shipped in `inst/extdata` is used.
#' @param package Package name used to locate the shipped reference.
#'
#' @details
#' The reference object stores the global per-organ FI means and standard
#' deviations, the frozen network parameters, and aggregated quantile tables
#' for the reference scores. It does not contain individual-level FI values or
#' reference score rows. No model fitting is performed at run time.
#'
#' @return An object of class `ot_reference` with elements:
#'   \describe{
#'     \item{version}{Reference version string.}
#'     \item{organs}{Canonical organ names.}
#'     \item{standardization}{Global per-organ FI means and standard
#'       deviations.}
#'     \item{networks}{`young` and `old` lists, each containing the fitted
#'       Gaussian parameters (`intercept`, `coef`, `sd`, `parents`) per organ,
#'       plus the arc list, roots and out-degree.}
#'     \item{quantiles}{Aggregated percentile tables for the overall and
#'       per-organ scores, stratified by sex and age band.}
#'     \item{validation}{Synthetic FI values and their expected scores, used by
#'       the package tests.}
#'     \item{provenance}{Source file names, MD5 checksums and notes.}
#'   }
#'
#' @examples
#' ref <- load_reference()
#' ref
#' ot_organs()
#'
#' fi <- data.frame(
#'   cardiovascular = 0.05, hepatic = 0.02, urinary = 0.04,
#'   endocrine = 0.03, immune = 0.06, musculoskeletal = 0.06,
#'   neural = 0.04, gastrointestinal = 0.03, respiratory = 0.05
#' )
#' transition_percentile(fi, ref, age = 68, sex = "woman")
#' @export
load_reference <- function(path = NULL, package = "OrganTransition") {
  if (is.null(path)) {
    path <- system.file("extdata", "reference.rds", package = package)
    if (!nzchar(path)) {
      local <- file.path("inst", "extdata", "reference.rds")
      if (file.exists(local)) {
        path <- local
      } else {
        stop("Cannot locate the shipped reference. Pass `path` explicitly.",
             call. = FALSE)
      }
    }
  }
  if (!file.exists(path)) {
    stop("Reference file not found: ", path, call. = FALSE)
  }
  ref <- readRDS(path)
  if (!inherits(ref, "ot_reference")) {
    stop("Not an OrganTransition reference object: ", path, call. = FALSE)
  }
  ref
}

#' Print an OrganTransition reference
#'
#' @param x An `ot_reference` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.ot_reference <- function(x, ...) {
  cat("<ot_reference> version ", x$version, "\n", sep = "")
  cat("  organs (", length(x$organs), "): ", paste(x$organs, collapse = ", "),
      "\n", sep = "")
  cat("  young network: ", x$networks$young$label,
      " | n = ", x$networks$young$n, "\n", sep = "")
  cat("  older network: ", x$networks$old$label,
      " | n = ", x$networks$old$n, "\n", sep = "")
  cat("  reference rows (aggregated): ", x$quantiles$all$n, "\n", sep = "")
  cat("  standardisation: global mean/sd\n")
  cat("  quantile tables: overall, sex, age 5/10/20 years\n")
  cat("  individual-level rows stored: no\n")
  invisible(x)
}
