test_that("Does the function stop when there is no jags installation?", {
  # Mock `detect_jags` to return FALSE
  with_mocked_bindings(
    "detect_jags" = function() FALSE,
    {
      expect_error(bma(bayesian_meta_analysis_fake_data), regexp = "No installation of JAGS has been detected")
    }
  )
})

test_that("Does it correctly detect missing columns in the input data?", {
  # skip the test if JAGS is not installed
  if (!detect_jags()) {
    skip("JAGS software not detectable")
  }

  expect_error(bma(subset(bayesian_meta_analysis_fake_data, select = -c(index))),
    regexp = "data column names should include"
  )
})

test_that("Does it correctly detect a missing standard error column?", {
  # skip the test if JAGS is not installed
  if (!detect_jags()) {
    skip("JAGS software not detectable")
  }

  expect_error(bma(subset(bayesian_meta_analysis_fake_data, select = -c(se)), seFromData = TRUE),
    regexp = "Standard errors have not been provided"
  )
})

test_that("Does it correctly return an error if a deprecated model, in this case smooth_det, is chosen?", {
  # skip the test if JAGS is not installed
  if (!detect_jags()) {
    skip("JAGS software not detectable")
  }

  expect_error(bma(bayesian_meta_analysis_fake_data, model = "smooth_det"),
    regexp = "model has been deprecated"
  )
})

test_that("Does the function output, assuming jags is available, return a dataframe?", {
  # skip the test if JAGS is not installed
  if (!detect_jags()) {
    skip("JAGS software not detectable")
  }

  expect_equal(class(suppressWarnings(bma(bayesian_meta_analysis_fake_data, plot = FALSE))), "data.frame")
})

test_that("Does the function return a Markov Chain Monte Carlo (MCMC) output object?", {
  # skip the test if jagsUI is not installed
  if (!detect_jags()) {
    skip("JAGS software not detectable")
  }

  expect_equal(class(suppressWarnings(bma_objects(bayesian_meta_analysis_fake_data, plot = TRUE)$comb.samples)), "mcmc.list")
})
