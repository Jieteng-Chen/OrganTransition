ref <- load_reference()
val <- ref$validation

test_that("standardisation maps raw [0, 1] FI to the reference scale", {
  z <- OrganTransition:::.standardize_fi(val$raw_fi, ref)
  expect_equal(as.matrix(z), as.matrix(val$z_fi), tolerance = 1e-12)
})

test_that("scoring reproduces expected scores from standardised FI", {
  scores <- OrganTransition:::.score_from_fi(val$z_fi, ref)
  expect_equal(unname(scores$transition_score),
               unname(val$expected$transition_score),
               tolerance = 1e-10)
  for (org in ref$organs) {
    col <- paste0("transition_score_", org)
    expect_equal(unname(scores[[col]]), unname(val$expected[[col]]),
                 tolerance = 1e-10)
  }
})

test_that("scoring output has one row per individual", {
  scores <- OrganTransition:::.score_from_fi(val$z_fi, ref)
  expect_equal(nrow(scores), nrow(val$z_fi))
  expect_equal(scores$transition_score,
               rowSums(as.matrix(scores[, paste0("transition_score_", ref$organs)])))
})

test_that("standardisation fails without parameters", {
  bad <- ref
  bad$standardization <- NULL
  expect_error(OrganTransition:::.standardize_fi(val$raw_fi, bad),
               "standardisation parameters")
})
