#' Rescale species values for indicator using Lambda interpolation
#' 
#' @description  This function takes in the output from a sparta occupancy model
#' or a three dimensional array. The first year is set to an index value for all
#' species and then using the average change for one year to the next, across
#' species is used to calculate the indicator. A species only contributes to the
#' dataset in years where the standard deviation across iterations
#' meets the sd threshold. Missing years in the middle a a species dataset are
#' filled in using interpolation. 
#' 
#' @param input Either a string giving the path to occupancy model output files
#' produced by sparta, or a three dimensional array [species, year, iteration].
#' @param index The index value for the first year, defaults to 100.
#' @param threshold_sd The threshold standard deviation for a species-year value
#' to be included. If the standard deviation is above this value it is removed.  
#' @param threshold_yrs Numeric, the minimum number of years that a species must
#' fulfill the threshold_sd for it to be included.
#' @param upperQuantile The upper confidence interval to use (as a probability)
#' @param lowerQuantile The lower confidence interval to use (as a probability)
#' @param sample_size numeric, if not NULL then a subsample of the iterations are 
#' used, equal to the number given. This is useful when datasets are so large
#' that memory starts to become limiting.
#' @param year_range, numic vector of length 2 giving the start and end year of the 
#' data to be analysed.
#' @return A list with five elements: a summary (data.frame), the LogLambda values
#' , calculated after removing species that fail thresholds and including
#' interpolation, the raw
#' indicator value (a value for each iteration in each year), the average annual
#' percentage change for each species, and a table giving the 'good' years for
#' each species as defined by the thresholds.
#' @importFrom car logit
#' @export
#' @examples 
#' ### Running from an array ####
#' # number of species
#' nsp = 50
#' 
#' # number of years
#' nyr = 40
#' 
#' #number of iterations
#' iter = 500
#' 
#' # Build a random set of data
#' myArray <- array(rnorm(n = nsp*nyr*iter,
#'                  mean = 0.5,
#'                  sd = 0.1),
#'                  c(nsp, nyr, iter))
#' 
#' # Ensure values are bounded by 0 and 1
#' myArray[myArray > 1] <- 1
#' myArray[myArray < 0] <- 0
#' 
#' # Run the lambda_interpolation method on this data                
#' myIndicator <- lambda_indicator(myArray)
#' 
#' # Plot the indicator
#' plot_indicator(myIndicator$summary[,'indicator'],
#'                myIndicator$summary[,c('lower' ,'upper')])
#' 
#' ### Running from a directory of sparta ouput
#' # myIndicator <- lambda_indicator('myfilepath/directory')

lambda_indicator <-  function(input, 
                              index = 100,
                              threshold_sd = 0.2,
                              threshold_yrs = 20,
                              upperQuantile = 0.975,
                              lowerQuantile = 0.025,
                              sample_size = NULL,
                              year_range = NULL){
  
  # Load the data if path else return input if array
  Occ <- getData(input = input, sample_size = sample_size)
  
  # Subset to years
  if(!is.null(year_range)) Occ <- subset_years(Occ, year_range)
  
  # How many species?
  nsp1 <- dim(Occ)[1]
  
  # Remove the species that do not meet our threshold conditions
  Occ <- remove_bad_species(Occ = Occ,
                            threshold_sd = threshold_sd,
                            threshold_yrs = threshold_yrs)
  
  # Save the reliable years table to return 
  good_years <- attr(Occ, 'good_years')
  
  # Now how many?
  nsp2 <- dim(Occ)[1]
  
  if(dim(Occ)[1] == 0){
    
    stop('No species meet the threshold criteria')
    
  } else if((nsp1 - nsp2) > 0){
    
    message(paste(nsp1 - nsp2, "Species have been removed as they don't meet your thresholds"))
    
  }
  
  # Convert to Odds
  Occ <- car::logit(Occ, adjust = 0.001) # Is this the best option?

  # Calculate the lambda for all species-year combinations
  LogLambda <- apply(Occ, c(1,3), lambda_calc)
  LogLambda <- aperm(LogLambda, c(2,1,3))
  
  # Add back in the dimnames
  dimnames(LogLambda) <- dimnames(Occ)
  
  # The indicator changes each year by the mean of LogLambda
  Delta <- apply(LogLambda, c(2,3), mean, na.rm = T)
  Delta[1,] <- 0
  Theta <- apply(Delta, 2, cumsum)
  
  # Create summary data to return (i.e. the indicator is simply Lambda rescaled to start at 100) 
  Lambda <- exp(Theta)
  Indicator_data <- index * Lambda
  
  indicator <- apply(Indicator_data, 1, mean)
  Indicator_CI <- t(apply(Indicator_data, 1, quantile,
                          probs = c(lowerQuantile, upperQuantile), 
                          na.rm = TRUE))
  colnames(Indicator_CI) <- c('lower', 'upper')
  
  summary_table <- as.data.frame(cbind(indicator, Indicator_CI))
  summary_table$year <- as.numeric(row.names(summary_table))
  
  if(NA %in% summary_table) warning('Data not available for all years in output due to threshold data removal')
 
  sp_change <- species_assessment(LogLambda)
  
  return(list(summary = summary_table,
              LogLambda = LogLambda,
              ind_data = Indicator_data,
              species_change = sp_change,
              good_years = good_years))
  
} 