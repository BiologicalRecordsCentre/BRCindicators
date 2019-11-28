#' Assess trends for indicator
#' 
#' @description  This assessment covers both the indicator and the constituant species.
#' An assessment is made of the indicator over the time period given, examining whether
#' the initial indicator value falls within the credibal inverval of the final year.
#' Over the same time period the change in each species is assessed and reported.
#' 
#' @param lambda_output An object returned by lambda_interpolation
#' @param start_year (Optional) a numeric value, defaults to the first year
#' @param end_year (Optional) a numeric value, defaults to the last year
#' @param species_stat (Optional) character, the statistic used to average across a species
#' yearly change values to arrive at a single value per species. This can be 
#' either 'mean' (default) or 'median'.
#' @return Returns a list of two elements, a summary of the species and indicator
#' assessments. A plot of the species assessment is returned to the device. 
#' @export
#' @examples 
#' ### Running from an array ####
#' set.seed(123)
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
#' myArray <- array(data = rnorm(n = nsp*nyr*iter,
#'                               mean = 0.5,
#'                               sd = 0.1),
#'                  dim = c(nsp, nyr, iter),
#'                  dimnames = list(paste0('SP',1:nsp),
#'                                  1:nyr,
#'                                  1:iter))
#' 
#' # Ensure values are bounded by 0 and 1
#' myArray[myArray > 1] <- 1
#' myArray[myArray < 0] <- 0
#' 
#' # Run the lambda_indicator method on this data                
#' myIndicator <- lambda_indicator(myArray)
#' 
#' # Plot the trend stack
#' trend_assessment(myIndicator)

trend_assessment <- function(lambda_output,
                             start_year = min(lambda_output$summary$year),
                             end_year = max(lambda_output$summary$year),
                             species_stat = 'mean'){
  
  sp_assess <- species_assessment(LogLambda = lambda_output$LogLambda,
                                  start_year = start_year + 1,
                                  end_year = end_year,
                                  species_stat = species_stat,
                                  plot = TRUE)
  
  ind_assessment <- indicator_assessment(summary_table = lambda_output$summary,
                                         start_year = start_year,
                                         end_year = end_year)
  
  return(list(species_assessment = sp_assess,
              indicator_asssessment = ind_assessment))
  
}