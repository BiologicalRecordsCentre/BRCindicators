#' Rescale species values for indicator
#' 
#' This function takes in a dataframe of multi-species data and rescales them
#' so that the value in the starting year is the same. This function accounts
#' for species that have no data at the beginning of the period or those that
#' have no data at teh end of the period.
#' 
#' @param input_dir the loction of occupancy model output files
#' @param subset_table dataframe with columns for years and rows for species. 
#' Column names must be character years and row names must match the file names
#' in input_dir
#' @param index The index value for the first year, defaults to 100.
#' @param max The upper limit allowed for scaled values. Values greater than 
#'        this will be set equal to \code{max}.
#' @param min The upper limit allowed for scaled values. Values greater than 
#'        this will be set equal to \code{min}.
#' @param year_limit The minimum length of 'good years' for a species to be
#'        included
#' @param upperCI The upper confidence interval to use (as a probability)
#' @param lowerCI The lower confidence interval to use (as a probability)
#' @param verbose If TRUE progress is written to console
#' @return A list with two elements, a summary and the rescaled data
#' @export
#' 
rescale_posterior <-  function(input_dir, subset_table = NULL,
                               index = 100, max = 10000,
                               min = 1, year_limit = 10, upperCI = 0.975,
                               lowerCI = 0.025, verbose = TRUE){
  
  # Read in all the data as a list
  library(reshape2)

  # get files from the input directory
  files <- list.files(path = paste(input_dir), ignore.case = TRUE, pattern = '\\.rdata$') # list of the files to loop through
  
  # sense check these file names
  if(length(files) == 0) stop('No .rdata files found in ', input_dir)
  if(length(files) < length(list.files(path = input_dir))) warning('Not all files in ', input_dir, ' are .rdata files, other file types have been ignored')
  
  if(verbose) cat('Loading data...')
  # Use lapply to run this function on all files
  list_summaries <- lapply(file.path(input_dir, files), read_posterior)
  names(list_summaries) <- files
  if(verbose) cat('done\n')
  
  # if we are useing a subset table get rid of all the data we dont need
  colnames(subset_table) <- gsub('^X', '', colnames(subset_table))
  
  if(!is.null(subset_table)){
    
    list_summaries <- lapply(names(list_summaries), remove_bad_years,
                             subset_table = subset_table, list_summaries = list_summaries)
    names(list_summaries) <- files
    
  }
  
  # Drop a species if it does not meet the number of years in time 
  # series criterea
  n_rows <- lapply(list_summaries, nrow)
  good_sp <- n_rows >= year_limit
  list_summaries <- list_summaries[good_sp]
  
  # Create a table for the first and last year of each spcies
  # NOTE this assumes there are no gaps in a species' time-series
  sp_periods <- as.data.frame(t(sapply(list_summaries, function(x) range(as.numeric(row.names(x))))))
  
  # Order by the entry year
  sp_periods <- sp_periods[order(sp_periods$V1),]
  
  # Take the 'entry species' and index to the index value
  first_year <- min(sp_periods[,1])
  initial_species <- row.names(sp_periods[sp_periods[,1] == first_year,])
  rescaled_list <- list_summaries[initial_species]
  rescaled_list <- entering_multiplier(data_list = rescaled_list,
                                       gm = index)
  # cap very high and very low values
  rescaled_list <- cap_index(rescaled_list, max = max, min = min)
  
  # Create summary output table
  summary_table <- data.frame(year = first_year:max(sp_periods[,2]),
                              indicator = NA, upperCI = NA,
                              lowerCI = NA, n_species = NA)
  summary_start <- list_geomean(data_list = rescaled_list, year = first_year, CI_min = lowerCI,
                                CI_max = upperCI)
  summary_table$indicator[summary_table$year == first_year] <- summary_start
  summary_table$upperCI[summary_table$year == first_year] <- attr(summary_start, 'CI_max')
  summary_table$lowerCI[summary_table$year == first_year] <- attr(summary_start, 'CI_min')
  summary_table$n_species[summary_table$year == first_year] <- attr(summary_start, 'n_species')
  
  # Now go through each year in turn.
  # First identify any species that leave and apply multiplier to remaining years
  # Second identfy new species and bring them in at the geomean
  for(i in (first_year + 1):max(sp_periods[,2])){
    
    # Find those species whose last years data was the year before
    leaving_species <- row.names(sp_periods[sp_periods[,2] == (i - 1),])
    
    if(length(leaving_species) > 0){
      # Get the geomean for the previous year with and without these species
      gm_with <- list_geomean(data_list = rescaled_list, year = i - 1)
      gm_without <- list_geomean(data_list = rescaled_list, year = i - 1,
                                 ignore_species = leaving_species)
      # Calculate the multiplication factor
      multiplier <- gm_with / gm_without
      
      # Apply this multiplier to all the species in all subsequent years (inc this one)
      rescaled_list <- leaving_multiplier(data_list = rescaled_list, year = i, multiplier = multiplier)
      
      # after multiplier cap very high and very low values
      rescaled_list <- cap_index(rescaled_list, max = max, min = min)
    }
    
    # identify all species entering the data set in this year
    entering_species <- row.names(sp_periods[sp_periods[,1] == i,])
    
    if(length(entering_species) > 0){
    
      # Get the geomean these species need to be set to
      gm <- list_geomean(data_list = rescaled_list, year = i)
      
      # rescale these species
      entering_data <- entering_multiplier(data_list = list_summaries[entering_species],
                                                       gm = gm)
      # cap very high and very low values
      entering_data <- cap_index(entering_data, max = max, min = min)
      
      # add these to the master data
      rescaled_list <- c(rescaled_list, entering_data)
    
    }
    
    summary <- list_geomean(data_list = rescaled_list, year = i, CI_min = lowerCI,
                            CI_max = upperCI)
    summary_table$indicator[summary_table$year == i] <- summary
    summary_table$upperCI[summary_table$year == i] <- attr(summary, 'CI_max')
    summary_table$lowerCI[summary_table$year == i] <- attr(summary, 'CI_min')
    summary_table$n_species[summary_table$year == i] <- attr(summary, 'n_species')
    
  }
  
return(list(summary = summary_table, data = rescaled_list))
  
}