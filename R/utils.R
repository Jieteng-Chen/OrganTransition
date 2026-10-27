# Internal helpers -----------------------------------------------------------

.organ_levels <- c(
  "cardiovascular", "hepatic", "urinary", "endocrine", "immune",
  "musculoskeletal", "neural", "gastrointestinal", "respiratory"
)

# Coerce user input to a data.frame.
.as_data_frame <- function(data) {
  if (is.null(data)) stop("`data` is NULL.", call. = FALSE)
  if (is.numeric(data) && !is.null(names(data))) {
    data <- as.data.frame(as.list(data), stringsAsFactors = FALSE)
  } else if (is.list(data) && !is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame, matrix, named vector or list.",
         call. = FALSE)
  }
  data
}

# Resolve organ FI columns. Only the canonical names are accepted.
.resolve_fi_columns <- function(data_names) {
  nm <- trimws(data_names)
  out <- character(0)
  for (org in .organ_levels) {
    hit <- which(nm == org)
    if (!length(hit)) next
    if (length(hit) > 1L) {
      stop("Duplicate column for organ '", org, "'.", call. = FALSE)
    }
    out[org] <- nm[hit]
  }
  missing <- setdiff(.organ_levels, names(out))
  if (length(missing)) {
    stop("Missing organ FI column(s): ", paste(missing, collapse = ", "),
         ". Allowed names: ", paste(.organ_levels, collapse = ", "),
         call. = FALSE)
  }
  out
}

# Prepare nine raw [0, 1] organ FI columns.
.prepare_fi <- function(data) {
  data <- .as_data_frame(data)
  resolved <- .resolve_fi_columns(names(data))
  fi <- data[, unname(resolved), drop = FALSE]
  names(fi) <- names(resolved)

  for (org in .organ_levels) {
    raw <- fi[[org]]
    num <- suppressWarnings(as.numeric(raw))
    bad <- is.na(num) & !is.na(raw)
    if (any(bad)) {
      stop("Non-numeric value(s) in column '", resolved[[org]], "'.",
           call. = FALSE)
    }
    fi[[org]] <- num
  }

  if (anyNA(as.matrix(fi))) {
    stop("All nine organ FI values are required.", call. = FALSE)
  }
  mat <- as.matrix(fi)
  if (any(!is.finite(mat))) {
    stop("Organ FI values must be finite.", call. = FALSE)
  }
  if (any(mat < 0 | mat > 1)) {
    warning("Organ FI values should lie in [0, 1].", call. = FALSE)
  }

  fi
}

# Map raw [0, 1] FI to the scale of the reference networks.
.standardize_fi <- function(fi, ref) {
  st <- ref$standardization
  if (is.null(st)) {
    stop("The reference does not contain standardisation parameters. ",
         "Reinstall the package or rebuild the reference object.",
         call. = FALSE)
  }
  out <- fi
  for (org in .organ_levels) {
    out[[org]] <- (fi[[org]] - st$mean[[org]]) / st$sd[[org]]
  }
  out
}

# Invert a quantile table: return 100 * P(X <= x).
.pct_from_quantiles <- function(x, q, probs) {
  if (!length(x)) return(numeric(0))
  vapply(x, function(xi) {
    if (is.na(xi)) return(NA_real_)
    if (xi <= q[1]) return(0)
    if (xi >= q[length(q)]) return(100)
    y <- stats::approx(q, probs, xout = xi, rule = 2, ties = "ordered")$y
    as.numeric(y) * 100
  }, numeric(1))
}

# Gaussian conditional log-likelihood of one node, given its parents.
.node_loglik <- function(data, node, params) {
  if (length(params$parents)) {
    mu <- params$intercept +
      as.numeric(as.matrix(data[, params$parents, drop = FALSE]) %*%
                   params$coef)
  } else {
    mu <- rep(params$intercept, nrow(data))
  }
  stats::dnorm(data[[node]], mean = mu, sd = params$sd, log = TRUE)
}

# Compute per-organ and overall transition scores from corrected FI.
.score_from_fi <- function(fi, ref) {
  organs <- ref$organs
  ll_old <- sapply(organs, function(o) {
    .node_loglik(fi, o, ref$networks$old$nodes[[o]])
  })
  ll_old <- matrix(ll_old, nrow = nrow(fi), dimnames = list(NULL, organs))
  ll_young <- sapply(organs, function(o) {
    .node_loglik(fi, o, ref$networks$young$nodes[[o]])
  })
  ll_young <- matrix(ll_young, nrow = nrow(fi), dimnames = list(NULL, organs))

  per <- ll_old - ll_young
  colnames(per) <- paste0("transition_score_", organs)
  out <- as.data.frame(per, check.names = FALSE)
  out$transition_score <- rowSums(per)
  out
}

# Resolve the `organ` argument: NULL -> none, "all"/"*" -> all organs,
# otherwise canonical organ names only.
.resolve_requested_organs <- function(organ) {
  if (is.null(organ) || !length(organ)) return(character(0))
  organ <- as.character(organ)
  if (length(organ) == 1L && tolower(trimws(organ)) %in% c("all", "*")) {
    return(.organ_levels)
  }
  bad <- setdiff(organ, .organ_levels)
  if (length(bad)) {
    stop("Unknown organ name(s): ", paste(bad, collapse = ", "),
         ". Allowed names: ", paste(.organ_levels, collapse = ", "),
         call. = FALSE)
  }
  unique(organ)
}

# Normalise sex labels to "man" / "woman" (the reference table uses these).
.normalise_sex <- function(x) {
  if (is.null(x)) return(rep(NA_character_, 0L))
  key <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(key))
  out[key %in% c("man", "male", "m", "1", "1.0", "1.man", "1.male")] <- "man"
  out[key %in% c("woman", "female", "f", "2", "2.0",
                 "2.woman", "2.female")] <- "woman"
  out
}

.percentile <- function(values, reference) {
  if (!length(reference)) return(rep(NA_real_, length(values)))
  as.numeric(stats::ecdf(reference)(values)) * 100
}
