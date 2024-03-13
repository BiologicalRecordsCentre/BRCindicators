test_that("Does function fail when given a non-matrix object?", {

df = as.data.frame(bootstrap_indicator_fake_data)

expect_error(bootstrap_indicator(df), "the Data parameter must be a matrix object")

})

test_that("Does function fail when given non-numeric values", {

letters_matrix <- matrix(sample(letters, 50 * length(letters), replace = TRUE), 
                                        nrow = 50, 
                                        ncol = length(letters))

# Assign the same column names
colnames(letters_matrix) <- letters

expect_error(bootstrap_indicator(letters_matrix), "Matrix values must all be numeric")

})

test_that("is the CI object ouput from sapply a matrix", {

expect_equal(class(bootstrap_indicator_objects(bootstrap_indicator_fake_data)$CIData), c("matrix", "array"))

})

test_that("is the CI object ouput from sapply a matrix", {

expect_equal(class(bootstrap_indicator(bootstrap_indicator_fake_data)), c("matrix", "array"))

})

test_that("is the CI object ouput from sapply a matrix", {

CI_limits = c(0.025, 0.975)

expect_equal(names(bootstrap_indicator(bootstrap_indicator_fake_data, CI_limits = CI_limits)), c("quant_025", "quant_975"))

})

test_that("is the CI object ouput from sapply a matrix", {

expect_equal(class(bootstrap_indicator(bootstrap_indicator_fake_data)), c("matrix", "array"))

})

test_that("Does the CI matrix output have the same column names as the CI_limits", {

CI_limits = c(0.025, 0.975)

expect_equal(colnames(bootstrap_indicator(bootstrap_indicator_fake_data, CI_limits = CI_limits)), c("quant_025", "quant_975"))

})