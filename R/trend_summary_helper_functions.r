### These functions are internally by a variety of the TrendSummary functions

geomean <- function(x){exp(mean(log(x)))}
quan_0.05 <- function(x) quantile(x, probs = 0.05, na.rm = TRUE)
quan_0.95 <- function(x) quantile(x, probs = 0.95, na.rm = TRUE)
quan_0.5 <- function(x) quantile(x, probs = 0.5, na.rm = TRUE)
cust_a_mean <-  function(x) mean(x, na.rm = T)
sem <- function(x) sd(x)/sqrt(length(x))