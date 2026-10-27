ref <- load_reference()

test_that("reference object has the expected structure", {
  expect_s3_class(ref, "ot_reference")
  expect_identical(ref$organs, ot_organs())
  expect_length(ref$organs, 9L)
  expect_true(all(c("young", "old") %in% names(ref$networks)))
  expect_true(all(ref$organs %in% names(ref$networks$young$nodes)))
  expect_true(all(ref$organs %in% names(ref$networks$old$nodes)))
  expect_false("scores" %in% names(ref))
})

test_that("network parameters are complete", {
  for (net in c("young", "old")) {
    for (org in ref$organs) {
      p <- ref$networks[[net]]$nodes[[org]]
      expect_true(is.numeric(p$intercept) && length(p$intercept) == 1L)
      expect_true(is.numeric(p$sd) && p$sd > 0)
      expect_identical(length(p$parents), length(p$coef))
      expect_true(all(p$parents %in% ref$organs))
    }
  }
})

test_that("global standardisation parameters are present", {
  st <- ref$standardization
  expect_true(all(c("mean", "sd") %in% names(st)))
  expect_true(all(ref$organs %in% names(st$mean)))
  expect_true(all(ref$organs %in% names(st$sd)))
  expect_true(all(is.finite(st$mean)))
  expect_true(all(is.finite(st$sd)))
  expect_true(all(st$sd > 0))
})

test_that("aggregated quantile tables are present", {
  q <- ref$quantiles
  expect_true(all(c("probs_all", "probs_cell", "scores", "all", "sex",
                    "age") %in% names(q)))
  expect_true("overall" %in% q$scores)
  expect_true(all(ref$organs %in% q$scores))
  expect_true(all(c("overall", ref$organs) %in% names(q$all$q)))
  expect_true(all(c("man", "woman") %in% names(q$sex)))
  expect_true(all(c("w5", "w10", "w20") %in% names(q$age)))
  expect_true(q$all$n > 10000)
  expect_true(all(diff(q$all$q$overall) >= 0))
})

test_that("no individual-level rows are stored", {
  expect_null(ref$scores)
  expect_true(all(c("raw_fi", "z_fi", "expected") %in% names(ref$validation)))
  expect_true(all(as.matrix(ref$validation$raw_fi) >= 0))
  expect_true(all(as.matrix(ref$validation$raw_fi) <= 1))
  expect_true("transition_score" %in% names(ref$validation$expected))
  expect_output(print(ref), "individual-level rows stored: no")
  expect_output(print(ref), "standardisation: global mean/sd")
})
