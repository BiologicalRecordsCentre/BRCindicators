test_that("Does function return correct error message when the user provides additional columns with character strings other than spp and iter", {
  fake_data$character_col_1 <- "some text"
  fake_data$character_col_2 <- "some text"

  # Create temporary file and dir
  input_path <- withr::local_tempfile(fileext = ".csv")
  output_path <- withr::local_tempdir()

  # Save the data to the temporary file
  write.csv(fake_data, file = input_path)

  expect_error(CompositeTrend(
    indata = input_path,
    output_path = output_path,
    trend_choice = "arithmetic_logit_occ",
    group_name = "TestGroup",
    save_iterations = "yes",
    TrendScale = 100,
    plot_output = FALSE), regexp = "Column\\(s\\) contains non-numeric fields")

})

test_that("Does the rbind of NULL value with a row work without error", {
  composite_trend_temp <- data.frame(
    X1985 = 1, X1986 = 1, X1987 = 1, X1988 = 1,
    X1989 = 1, X1990 = 1, X1991 = 1, X1992 = 1, X1993 = 1, X1994 = 1
  )

  expect_equal(nrow(rbind(NULL, composite_trend_temp)), nrow(composite_trend_temp))
})

test_that("Does the composite_trend object have the same number of rows as there are interations, and the same columns as there are years?", {
  
  # Create temporary file and dir
  input_path <- withr::local_tempfile(fileext = ".csv")
  output_path <- withr::local_tempdir()

  # Save the data to the temporary file
  write.csv(fake_data, file = input_path)

  composite_trend <- CompositeTrend_objects(
    indata = input_path,
    output_path = output_path,
    trend_choice = "arithmetic_logit_occ",
    group_name = "TestGroup",
    save_iterations = "yes",
    TrendScale = 100,
    plot_output = FALSE
  )$composite_trend

  expect_equal(dim(composite_trend)[1], n_iterations)
  expect_equal(dim(composite_trend)[2], n_years)
})

test_that("Is the composite_trend object a dataframe?", {

  # Create temporary file and dir
  input_path <- withr::local_tempfile(fileext = ".csv")
  output_path <- withr::local_tempdir()

  # Save the data to the temporary file
  write.csv(fake_data, file = input_path)

  composite_trend <- CompositeTrend_objects(
    indata = input_path,
    output_path = output_path,
    trend_choice = "arithmetic_logit_occ",
    group_name = "TestGroup",
    save_iterations = "yes",
    TrendScale = 100,
    plot_output = FALSE
  )$composite_trend

  expect_equal(class(composite_trend), c("matrix", "array"))
})

test_that("Is the composite_trend_summary object a dataframe?", {

  # Create temporary file and dir
  input_path <- withr::local_tempfile(fileext = ".csv")
  output_path <- withr::local_tempdir()

  # Save the data to the temporary file
  write.csv(fake_data, file = input_path)

  composite_trend_summary <- CompositeTrend_objects(
    indata = input_path,
    output_path = output_path,
    trend_choice = "arithmetic_logit_occ",
    group_name = "TestGroup",
    save_iterations = "yes",
    TrendScale = 100,
    plot_output = FALSE
  )$composite_trend_summary

  expect_equal(class(composite_trend_summary), "data.frame")
})


test_that("Is the composite_trend_summary object a ggplot?", {

  # Create temporary file and dir
  input_path <- withr::local_tempfile(fileext = ".csv")
  output_path <- withr::local_tempdir()

  # Save the data to the temporary file
  write.csv(fake_data, file = input_path)

  plot <- CompositeTrend_objects(
    indata = input_path,
    output_path = output_path,
    trend_choice = "arithmetic_logit_occ",
    group_name = "TestGroup",
    save_iterations = "yes",
    TrendScale = 100,
    plot_output = TRUE
  )$plot

  print(class(plot))

  expect_equal(class(plot), c("gg", "ggplot"))
})
