ref <- load_reference()
val <- ref$validation
q <- ref$quantiles

test_that("percentiles are returned without raw transition scores", {
  p <- transition_percentile(val$raw_fi[1:3, , drop = FALSE], ref,
                             age = rep(68, 3),
                             sex = c("woman", "man", "woman"))
  expect_false(any(grepl("^transition_score", names(p))))
  expect_true(all(c("percentile_all", "percentile_age_sex") %in% names(p)))
  expect_true(all(p$percentile_all >= 0 & p$percentile_all <= 100))
})

test_that("overall percentile matches the quantile table", {
  i <- 1
  p <- transition_percentile(val$raw_fi[i, , drop = FALSE], ref,
                             age = 68, sex = "woman", organ = NULL)
  z <- OrganTransition:::.standardize_fi(val$raw_fi[i, , drop = FALSE], ref)
  s <- OrganTransition:::.score_from_fi(z, ref)
  manual <- OrganTransition:::.pct_from_quantiles(
    s$transition_score, q$all$q$overall, q$probs_all
  )
  expect_equal(p$percentile_all, manual, tolerance = 1e-8)
})

test_that("per-organ percentiles match the quantile table by default", {
  i <- 2
  p <- transition_percentile(val$raw_fi[i, , drop = FALSE], ref,
                             age = 68, sex = "woman")
  z <- OrganTransition:::.standardize_fi(val$raw_fi[i, , drop = FALSE], ref)
  s <- OrganTransition:::.score_from_fi(z, ref)
  for (org in ref$organs) {
    col <- paste0("transition_score_", org)
    manual <- OrganTransition:::.pct_from_quantiles(
      s[[col]], q$all$q[[org]], q$probs_all
    )
    expect_equal(p[[paste0("percentile_all_", org)]], manual,
                 tolerance = 1e-8)
    expect_true(paste0("percentile_age_sex_", org) %in% names(p))
  }
})

test_that("organ = NULL returns the overall percentile only", {
  p <- transition_percentile(val$raw_fi[1, , drop = FALSE], ref,
                             age = 68, sex = "woman", organ = NULL)
  expect_false(any(grepl("^percentile_all_", names(p))))
})

test_that("age-sex matching returns a reference cell", {
  p <- transition_percentile(val$raw_fi[1:5, , drop = FALSE], ref,
                             age = c(40, 50, 60, 70, 80),
                             sex = c("man", "woman", "man", "woman", "man"))
  expect_true(all(p$n_age_sex > 0))
  expect_true(all(p$matched_by %in% c("sex+age5", "sex+age10", "sex+age20",
                                      "age5", "age10", "age20",
                                      "sex", "all")))
})

test_that("small cells fall back to the full reference", {
  p <- transition_percentile(val$raw_fi[1, , drop = FALSE], ref,
                             age = 68, sex = "woman", min_n = 1e9)
  expect_identical(p$matched_by, "all")
  expect_equal(p$n_age_sex, q$all$n)
})

test_that("sex labels are normalised", {
  a <- transition_percentile(val$raw_fi[1, , drop = FALSE], ref,
                             age = 70, sex = "male")
  b <- transition_percentile(val$raw_fi[1, , drop = FALSE], ref,
                             age = 70, sex = "man")
  expect_equal(a$percentile_age_sex, b$percentile_age_sex)
  expect_identical(a$sex, "man")
})

test_that("sex-only matching works without age", {
  p <- transition_percentile(val$raw_fi[1, , drop = FALSE], ref, sex = "woman")
  expect_identical(p$matched_by, "sex")
  expect_true(p$n_age_sex > 0)
})

test_that("missing age and sex warns and returns NA age-sex percentiles", {
  expect_warning(
    p <- transition_percentile(val$raw_fi[1, , drop = FALSE], ref),
    "needs at least one"
  )
  expect_true(is.na(p$percentile_age_sex))
})

test_that("strict organ names are enforced", {
  bad <- val$raw_fi[1, -1, drop = FALSE]
  expect_error(transition_percentile(bad, ref, age = 70, sex = "woman"),
               "Missing organ FI")
  bad2 <- val$raw_fi[1, , drop = FALSE]
  names(bad2)[names(bad2) == "gastrointestinal"] <- "digestive"
  expect_error(transition_percentile(bad2, ref, age = 70, sex = "woman"),
               "Missing organ FI")
})

test_that("FI outside [0, 1] triggers a warning", {
  bad <- val$raw_fi[1, , drop = FALSE]
  bad[1, 1] <- 1.5
  expect_warning(
    transition_percentile(bad, ref, age = 70, sex = "woman"),
    "\\[0, 1\\]"
  )
})

test_that("age must have one value per individual", {
  expect_error(
    transition_percentile(val$raw_fi[1:2, , drop = FALSE], ref,
                          age = 70, sex = c("man", "woman")),
    "one value per individual"
  )
})
