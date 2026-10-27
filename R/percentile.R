#' Percentiles of organ-specific transition scores
#'
#' Takes nine raw `[0, 1]` organ FI values for one or more individuals,
#' standardises them internally with frozen global means and standard
#' deviations, computes the overall and per-organ transition scores with the
#' reference networks, and returns their percentiles. Raw transition-score
#' values are not returned.
#'
#' Two percentiles are produced for the overall score and for each requested
#' organ: `percentile_all` ranks against the full reference cohort, and
#' `percentile_age_sex` ranks against reference participants of the same sex
#' and a nearby age band. Percentiles are defined as `100 * P(X <= x)` and are
#' obtained from aggregated quantile tables; no individual-level reference rows
#' are used at run time.
#'
#' @param data A data.frame, matrix, named numeric vector or list containing one
#'   row per individual and nine raw `[0, 1]` organ FI columns named exactly
#'   `cardiovascular`, `hepatic`, `urinary`, `endocrine`, `immune`,
#'   `musculoskeletal`, `neural`, `gastrointestinal`, `respiratory`.
#' @param ref A reference object from [load_reference()].
#' @param age Optional numeric age vector, one value per individual. Age is used
#'   only for age-sex percentile matching.
#' @param sex Optional sex vector, one value per individual. Accepts `"man"`,
#'   `"woman"`, `"male"`, `"female"`, `"m"`, `"f"`, `1` (man) or `2` (woman).
#' @param organ `"all"` (default) for the overall score plus all nine organs,
#'   a character vector of organ names, or `NULL` for the overall score only.
#' @param by `"all"` for the full-reference percentile only, or `"age_sex"`
#'   (default) for both. `"age_sex"` needs at least one of `age` and `sex`;
#'   otherwise the age-sex matched percentiles are `NA`.
#' @param min_n Minimum reference cell size required for an age-sex match. Cells
#'   smaller than this are skipped in favour of a wider band, a sex-only match,
#'   or the full reference.
#'
#' @return A data.frame with `percentile_all`, `percentile_age_sex`, `n_all`,
#'   `n_age_sex`, `matched_by`, `age_band`, `age`, `sex`, and, for each
#'   requested organ, `percentile_all_<organ>` and
#'   `percentile_age_sex_<organ>`. No raw transition-score values are returned.
#'
#' @details
#' Raw FI is standardised internally as `(FI - mean) / sd`, using global
#' per-organ means and standard deviations frozen from the reference cohort.
#' Percentiles come from aggregated quantile tables for the overall score and
#' for each organ, stratified by sex and by 5-, 10- and 20-year age bands.
#' Reference cells with fewer than 100 individuals are not stored. No
#' individual-level FI values or reference scores are included in the package.
#'
#' @examples
#' ref <- load_reference()
#'
#' fi <- data.frame(
#'   cardiovascular = 0.05, hepatic = 0.02, urinary = 0.04,
#'   endocrine = 0.03, immune = 0.06, musculoskeletal = 0.06,
#'   neural = 0.04, gastrointestinal = 0.03, respiratory = 0.05
#' )
#'
#' # Overall and all nine organ percentiles
#' transition_percentile(fi, ref, age = 68, sex = "woman")
#'
#' # Overall only
#' transition_percentile(fi, ref, age = 68, sex = "woman", organ = NULL)
#' @export
transition_percentile <- function(data, ref = load_reference(), age = NULL,
                                  sex = NULL, organ = "all",
                                  by = "age_sex", min_n = 100) {
  by <- match.arg(by, c("age_sex", "all"))
  if (!inherits(ref, "ot_reference")) {
    stop("`ref` must be an ot_reference object from load_reference().",
         call. = FALSE)
  }
  if (!is.numeric(min_n) || length(min_n) != 1L || min_n < 1) {
    stop("`min_n` must be a positive number.", call. = FALSE)
  }

  fi <- .prepare_fi(data)
  n <- nrow(fi)
  z <- .standardize_fi(fi, ref)

  age_used <- if (is.null(age)) rep(NA_real_, n) else as.numeric(age)
  if (length(age_used) != n) {
    stop("`age` must have one value per individual.", call. = FALSE)
  }
  sex_used <- if (is.null(sex)) rep(NA_character_, n) else .normalise_sex(sex)
  if (length(sex_used) != n) {
    stop("`sex` must have one value per individual.", call. = FALSE)
  }

  requested <- .resolve_requested_organs(organ)

  scores <- .score_from_fi(z, ref)
  values <- scores$transition_score
  organ_values <- list()
  for (org in requested) {
    organ_values[[org]] <- scores[[paste0("transition_score_", org)]]
  }

  quant <- ref$quantiles
  if (is.null(quant)) {
    stop("The reference does not contain quantile tables. Reinstall the ",
         "package or rebuild the reference object.", call. = FALSE)
  }

  all_cell <- quant$all
  percentile_all <- .pct_from_quantiles(
    values, all_cell$q[["overall"]], quant$probs_all
  )
  n_all <- all_cell$n
  organ_all <- list()
  for (org in requested) {
    organ_all[[org]] <- .pct_from_quantiles(
      organ_values[[org]], all_cell$q[[org]], quant$probs_all
    )
  }

  percentile_age_sex <- rep(NA_real_, n)
  organ_age_sex <- list()
  for (org in requested) organ_age_sex[[org]] <- rep(NA_real_, n)
  n_age_sex <- rep(NA_integer_, n)
  matched_by <- rep(NA_character_, n)
  age_band <- rep(NA_character_, n)

  if (by == "age_sex") {
    if (is.null(age) && is.null(sex)) {
      warning("`by = \"age_sex\"` needs at least one of `age` and `sex`; ",
              "age-sex matched percentiles are returned as NA.", call. = FALSE)
    } else {
      widths <- c(5, 10, 20)

      for (i in seq_len(n)) {
        cell <- NULL
        label <- NULL
        band_label <- NA_character_
        count <- NA_integer_
        is_all <- FALSE

        if (!is.na(age_used[i])) {
          sx <- if (is.na(sex_used[i])) "all" else sex_used[i]
          for (w in widths) {
            lo <- floor(age_used[i] / w) * w
            key <- as.character(lo)
            cand <- quant$age[[paste0("w", w)]][[sx]][[key]]
            if (!is.null(cand) && cand$n >= min_n) {
              cell <- cand
              label <- paste0(if (is.na(sex_used[i])) "" else "sex+", "age", w)
              band_label <- sprintf("[%g,%g)", lo, lo + w)
              count <- cand$n
              break
            }
          }
        }

        if (is.null(cell) && !is.na(sex_used[i])) {
          cand <- quant$sex[[sex_used[i]]]
          if (!is.null(cand) && cand$n >= min_n) {
            cell <- cand
            label <- "sex"
            count <- cand$n
          }
        }

        if (is.null(cell)) {
          cell <- all_cell
          label <- "all"
          count <- all_cell$n
          is_all <- TRUE
        }
        cell_probs <- if (is_all) quant$probs_all else quant$probs_cell

        percentile_age_sex[i] <- .pct_from_quantiles(
          values[i], cell$q[["overall"]], cell_probs
        )
        for (org in requested) {
          organ_age_sex[[org]][i] <- .pct_from_quantiles(
            organ_values[[org]][i], cell$q[[org]], cell_probs
          )
        }
        matched_by[i] <- label
        age_band[i] <- band_label
        n_age_sex[i] <- count
      }
    }
  }

  out <- data.frame(
    percentile_all     = percentile_all,
    percentile_age_sex = percentile_age_sex,
    n_all              = n_all,
    n_age_sex          = n_age_sex,
    matched_by         = matched_by,
    age_band           = age_band,
    age                = age_used,
    sex                = sex_used,
    check.names        = FALSE,
    stringsAsFactors   = FALSE
  )

  for (org in requested) {
    out[[paste0("percentile_all_", org)]] <- organ_all[[org]]
    out[[paste0("percentile_age_sex_", org)]] <- organ_age_sex[[org]]
  }

  out
}
