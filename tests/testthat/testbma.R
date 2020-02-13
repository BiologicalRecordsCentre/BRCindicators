context("Test bma")

set.seed(123)

# Create some example data in the format required
data <- data.frame(species = rep(letters, each = 50),
                   year = rep(1:50, length(letters)),
                   index = rnorm(n = 50 * length(letters), mean = 0, sd = 1),
                   se = runif(n = 50 * length(letters), min = 0.01, max = .1))

test_that("runs without error", {
  
  # Run the Bayesian meta-analysis
  bma_indicator <- bma(data, model="smooth", m.scale="logit", n.iter=100)
  
  expect_is(bma_indicator, 'data.frame')
  expect_equal(ncol(bma_indicator), 4)
  
})