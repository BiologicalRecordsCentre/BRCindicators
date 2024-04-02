test_that("runs without error", {
  
  # Run the lambda_interpolation method on this data
  myIndicator <- lambda_indicator(myArray)
  expect_equal(class(myIndicator), 'list')
  expect_equal(length(myIndicator), 5)
  
})
