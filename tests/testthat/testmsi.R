context("Test msi")

set.seed(123)

# Create some example data in the format required
nyr = 20
species = rep(letters, each = nyr)
year = rev(rep(1:nyr, length(letters)))

# Create an index value that increases with time
index = rep(seq(50, 100, length.out = nyr), length(letters))
# Add randomness to species
index = index * runif(n = length(index), 0.7, 1.3)
# Add correlated randomness acrosss species, to years
index = index * rep(runif(0.8, 1.2, n = nyr), length(letters))

se = runif(n = nyr * length(letters), min = 10, max = 20)

data <- data.frame(species, year, index, se)

# Species index values need to be 100 in the base year. Here I use
# the first year as my base year and rescale to 100. The standard error
# in the base year should be 0.
min_year <- min(data$year)

for(sp in unique(data$species)){

  subset_data <- data[data$species == sp, ]
  multi_factor <- 100 / subset_data$index[subset_data$year == min_year]
  data$index[data$species == sp] <- data$index[data$species == sp] * multi_factor
  data$se[data$species == sp] <- data$se[data$species == sp] * multi_factor
  data$se[data$species == sp][1] <- 0

}

test_that("runs without error", {
  
  # Run the MSI function
  msi_out <- msi(data, plot = FALSE)
  
  expect_is(msi_out, 'MSI')
  expect_equal(length(msi_out), 3)
  
})