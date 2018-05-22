#' Multi-Species Indicator
#' 
#' @description A simple wrapper for \code{msi_tool} to make it easier to use in R.
#' Multi-Species Indicators (MSI) are biodiversity indicators that combine the population
#' development of species into a single indicator. The MSI-tool calculates an MSI, confidence intervals
#' for the MSIs and linear and flexible (smoothed) trends. The trends are classified in terms like
#' "moderate increase", "strong decrease" or "stable". A number of additional analyses can be performed
#' like testing for changepoints, comparison of trends before and after a changepoint and the calculation
#' and testing of the total change in a time series.
#' 
#' @param data a data.frame with 4 columns in this order: 'species', 'year', 'index', 'se' (standard error) 
#' @param ... other parameters to pass to \code{msi_tool}
#' @return Returns a dataframe with 4 columns: Year, Index, lower2.5, upper97.5. The last two columns are the credible intervals
#' @import reshape2
#' @importFrom boot inv.logit
#' @export
#' @examples 
#' 
#' # Create some example data in the format required
#' data <- data.frame(species = rep(letters, each = 50),
#'                    year = rep(1:50, length(letters)),
#'                    index = runif(n = 50 * length(letters), min = 1, max = 10),
#'                    se = runif(n = 50 * length(letters), min = 0.01, max = .8))
#' 
#' # Run the MSI function
#' msi_out <- msi(data, plot = FALSE)
#' 
#' # Plot the resulting indicator
#' plot(msi_out)

msi <- function(data, ...){
  
  stopifnot(inherits(x = data, what = 'data.frame'))
  if(!all(colnames(data) == c("species", "year", "index", "se"))){
    stop('column names must be "species" "year" "index" "se"')
  }
  
  # order the data by year
  data <- data[order(data$year), ]
  
  # check column types
  stopifnot(inherits(data$year, 'integer') | inherits(data$year, 'numeric'))
  stopifnot(inherits(data$index, 'integer') | inherits(data$index, 'numeric'))
  stopifnot(inherits(data$se, 'integer') | inherits(data$se, 'numeric'))
  stopifnot(inherits(data$species, 'character') | inherits(data$species, 'factor'))
  if(inherits(data$species, 'factor')) data$species <- as.character(data$species)
  
  # data needs to be rescaled to 100
  #for(sp in unique(data$species)){
    
  #  multi_factor <- 100 / na.omit(data$index[data$species == sp])[1]
  #  data$index[data$species == sp] <- data$index[data$species == sp] * multi_factor
  #  data$se[data$species == sp] <- data$se[data$species == sp] * multi_factor
    
  #}
  
  dir <- tempdir()

  write.csv(data, file = file.path(dir, 'input.csv'), row.names = FALSE)
  
  msi_tool(wd = dir, inputFile = 'input.csv', jobname = 'MSI_job', ...)
  
  results <- read.table(file.path(dir, "MSI_job_RESULTS.csv"), sep = ';',
                        header = TRUE)
  # replace all commas with decimals and make them numbers
  # who does that?!
  for(i in 2:8){
    results[,i] <- as.numeric(gsub(',','.',as.character(results[,i])))
  }
  
  trends <- read.table(file.path(dir, "MSI_job_TRENDS.csv"), sep = ';',
                       header = TRUE)
  colnames(trends)[1] <- 'Measure'
  trends$value <- as.numeric(gsub(',','.',as.character(trends$value)))
  
  msi_out <- list(results = results, 
                  trends = trends)
  
  class(msi_out) <- 'MSI'
  
  return(msi_out)
}