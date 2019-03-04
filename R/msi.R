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
#' @param data a data.frame with 4 columns in this order: 'species', 'year', 'index', 'se' (standard error).
#' The index value in the base year (which need not be the first year), should be set to 100, with se of 0. 
#' @param ... other parameters to pass to \code{msi_tool}
#' @return Returns a dataframe with 4 columns: Year, Index, lower2.5, upper97.5. The last two columns are the credible intervals
#' @import reshape2
#' @importFrom boot inv.logit
#' @export
#' @examples 
#' 
#' # Create some example data in the format required
#' Create some example data in the format required
#' nyr = 20
#' species = rep(letters, each = nyr)
#' year = rev(rep(1:nyr, length(letters)))
#' 
#' # Create an index value that increases with time
#' index = rep(seq(50, 100, length.out = nyr), length(letters))
#' # Add randomness to species
#' index = index * runif(n = length(index), 0.7, 1.3)
#' # Add correlated randomness acrosss species, to years
#' index = index * rep(runif(0.8, 1.2, n = nyr), length(letters))
#' 
#' se = runif(n = nyr * length(letters), min = 10, max = 20)
#' 
#' data <- data.frame(species, year, index, se)
#' 
#' # Our species are decreasing
#' plot(data$year, data$index)
#' 
#' # Species index values need to be 100 in the base year. Here I use
#' # the first year as my base year and rescale to 100. The standard error
#' # in the base year should be 0.
#' min_year <- min(data$year)
#' 
#' for(sp in unique(data$species)){
#'   
#'   subset_data <- data[data$species == sp, ]
#'   multi_factor <- 100 / subset_data$index[subset_data$year == min_year]
#'   data$index[data$species == sp] <- data$index[data$species == sp] * multi_factor
#'   data$se[data$species == sp] <- data$se[data$species == sp] * multi_factor
#'   data$se[data$species == sp][1] <- 0
#'   
#' }
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
  
  # The base year index value should be 100 and se in this year should be 0
  # Here are a couple of warnings to try and pick up when the user has got this wrong
  all_species <- unique(data$species)
  species_with_100_index <- unique(data$species[round(data$index, digits = 5) == 100])
  species_with_0_se <- unique(data$species[data$se == 0])
  
  if(!(all(all_species %in% species_with_100_index))){
    warning('Species are expected to have an index value of 100 in their base year. Some of',
            ' your species do not have any index values of 100: ',
            paste(all_species[!all_species %in% species_with_100_index], collapse = ', '))
  }
  if(!(all(all_species %in% species_with_0_se))){
    warning('Species are expected to have an se value of 0 in their base year (where index is',
            ' set to 100).',
            ' Some of your species do not have any se values',
            ' of 0: ',
            paste(all_species[!all_species %in% species_with_0_se], collapse = ', '))
  }  
  
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

  CV <- read.table(file.path(dir, "species_CV_values.csv"), sep = ',',
                   header = TRUE)
    
  msi_out <- list(results = results, 
                  trends = trends,
                  CV = CV)
  
  class(msi_out) <- 'MSI'
  
  return(msi_out)
}