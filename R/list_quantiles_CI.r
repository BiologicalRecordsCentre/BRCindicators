# A function to get the quantiles for each year
list_quantiles_CI <- function(rescaled_list, years, quantile_min = 0.025, quantile_max = 0.975){
  
  cat('Calculating quantiles and confidence intervals...')
  
  yr_val <- function(year, rescaled_list){
    
    year_vals <- lapply(rescaled_list, FUN = function(x){
          if(as.character(year) %in% row.names(x)){
            return(x[as.character(year),])
          } else {
            return(NULL)
          }
        })

    # Calculate the quantiles of all the posterior points from this year
    data_quantiles <- quantile(unlist(year_vals), probs = c(quantile_min, quantile_max), na.rm = TRUE)
    names(data_quantiles) <- paste('quantile',
                                   gsub('%', '', names(data_quantiles)),
                                   sep = '_')
    
    # rbind all the species from this year together, giving us a table, 
    # each row a species, each column an interation
    if(is.list(year_vals)){
      
      yr_iterations <- do.call(rbind, year_vals)
    
    } else { #if only one species is present
      
      yr_iterations <- t(as.data.frame(year_vals))
      
    }
    
    # apply geomean across columns, giving one geomean for each iteration
    yr_iteration_geomeans <- apply(X = yr_iterations, MARGIN = 2, FUN = geomean)

    # take the 95% CIs from these geomeans
    geomean_quantiles <- quantile(yr_iteration_geomeans, probs = c(0.025, 0.975), na.rm = TRUE)
    names(geomean_quantiles) <- paste('geomean_CI',
                                      gsub('%', '', names(geomean_quantiles)),
                                      sep = '_')
    # combine this data in a easy to handle format
    range_data <- cbind(t(as.data.frame(geomean_quantiles)), t(as.data.frame(data_quantiles)))
    row.names(range_data) <- year
    return(range_data)
    
  }
  
  range_summaries <- do.call(rbind, lapply(years, FUN = function(x) yr_val(year = x, rescaled_list = rescaled_list)))
  
  cat('done\n')
  
  return(range_summaries)
  
}