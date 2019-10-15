context("Test simulate_indicator")

test_that("Test simulate_indicator", {
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  set.seed(122)

  results<- simulate_indicator(nyears=30, nsp=8, mg= -0.02, sdg=0.2, sigma=0.1, SE=0.05)
  
  head_results<-structure(list(species = c(1L, 2L, 3L, 4L, 5L, 6L), year = c(1L, 1L, 1L, 1L, 1L, 1L),
                                index = c(-0.08350811, 0.02034243, 0.77392802, -0.29234623, 0.20268519, 0.16993784)),
                    .Names = c("species", "year", "index"), class = "data.frame")
  expect_equivalent(head(results), head_results)
  
  sink()
})
  