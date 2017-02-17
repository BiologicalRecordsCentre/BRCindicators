#' Bootstrapping indicator
#' 
#' This function takes in a dataframe of multi-species data that have been scaled
#' and adds confidence intervals to the indicator value.
#' 
#' @param Data A matrix where each named column give the species' values
#'        for each year in rows.
#' @param iterations The number of bootstrap iterations to use.
#' @param CI_limits The confidence limits to return as a vector of two
#'        numbers. This default to c(0.025, 0.975) the 95 percent conficence intervals.
#' @param verbose If \code{TRUE} then progress is printed to screen.
#' @return A matrix. In each row we have the year, the species and the
#'         scaled value. There is also an additional column, 'geomean' giving
#'         the geometric mean for each year.
#' @export

bootstrap_indicator <-  function(Data, iterations = 10000, CI_limits = c(0.025, 0.975),
                                 verbose = TRUE){
  geomean <- function(x) exp(mean(log(x), na.rm = T))
  
  ### bootstrap to estimate 95% CI around the geometric mean ###
  nSpecies <- ncol(Data)
  bootstrap_values <- matrix(data = NA, nrow = nrow(Data), ncol = iterations)
  bootstrap_data <- as.data.frame(Data)

  BS <- function(x){ # function to do bootstrapping, x is data
    
    # Randomly sample species 
    samp <- sample(x, size = ncol(x), replace = TRUE)    
    gm <- apply(samp, 1, geomean)
    
  }

  # Run bootstrapping
  if(verbose) cat('Running bootstrapping for', iterations, 'iterations...')
  CIData <- replicate(n = iterations, BS(bootstrap_data))
  if(verbose) cat('done\n')

  # Extract the 2.5 97.5% CI around the geometric mean (from the bootstrapped resamples)
  CIs <- as.data.frame(t(apply(X = CIData, MARGIN = 1, FUN = quantile,
                               probs = CI_limits, na.rm = TRUE)))
  names(CIs) <- c(paste('quant_', gsub('0\\.', '', as.character(CI_limits[1])), sep = ''),
                  paste('quant_', gsub('0\\.', '', as.character(CI_limits[2])), sep = ''))

  return(as.matrix(CIs))
}
  