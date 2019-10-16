context("Test trend_summary_helper_functions")

test_that("Test trend_summary_helper_functions", {
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  set.seed(122)
  
  testdata<- rnorm(1000,0.5,0.1)
  
  gm_result <- geomean(testdata)
  
  expect_equal(gm_result, 0.4917808, tolerance = 1e-7)
  
  quanL_result <- quan_0.05(testdata)
  
  expect_equivalent(quanL_result, 0.3222446, tolerance = 1e-7)
  
  quanM_result <- quan_0.5(testdata)
  
  expect_equivalent(quanM_result, 0.5055763, tolerance = 1e-7)
  
  quanU_result <- quan_0.95(testdata)
  
  expect_equivalent(quanU_result, 0.669175, tolerance = 1e-7)
  
  cam_result <- cust_a_mean(testdata)
  
  expect_equivalent(cam_result, 0.5033326, tolerance = 1e-7)
  
  sem_result <- sem(testdata)
  
  expect_equivalent(sem_result, 0.003236797, tolerance = 1e-7)
  
  sink()
})