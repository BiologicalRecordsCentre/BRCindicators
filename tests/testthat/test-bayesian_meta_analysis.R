# Test that function responds correctly when 'jagsUI' is not available
test_that("Does the function stop when 'jagsUI' is not available?", {

  # Mock `requireNamespace` to return FALSE
  mockery::stub(bma, "requireNamespace", FALSE)

  expect_error(bma(bayesian_meta_analysis_fake_data), regexp = "Package 'jagsUI' is needed")

})

test_that("Does it correctly detect missing columns in the input data?", {

expect_error(bma(subset(bayesian_meta_analysis_fake_data, select = -c(index))),
 regexp = "data column names should include")
})

test_that("Does it correctly detect a missing standard error column?", {

expect_error(bma(subset(bayesian_meta_analysis_fake_data, select = -c(se)), seFromData = TRUE),
 regexp = "Standard errors have not been provided")
})

test_that("Does it correctly return an error if a deprecated model, in this case smooth_det, is chosen?", {

expect_error(bma(bayesian_meta_analysis_fake_data, model = "smooth_det"),
 regexp = "model has been deprecated")
})

test_that("Does the function output, assuming 'jagsUI' is available, return a dataframe?", {

    # skip the test if jagsUI is not installed
    if ((!requireNamespace("jagsUI", quietly = TRUE))) {
    skip("jagsUI software not available")
  }

  expect_equal(class(suppressWarnings(bma(bayesian_meta_analysis_fake_data, plot = FALSE))), "data.frame")

})

test_that("Does the function return a Markov Chain Monte Carlo (MCMC) output object?", {

    # skip the test if jagsUI is not installed
    if ((!requireNamespace("jagsUI", quietly = TRUE))) {
    skip("jagsUI software not available")
  }

  expect_equal(class(suppressWarnings(bma_objects(bayesian_meta_analysis_fake_data, plot = TRUE)$comb.samples)), "mcmc.list")

})
