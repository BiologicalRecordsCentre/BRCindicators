# Create a function that will calculate the geomean from the list
# of species for a given year
list_geomean <- function(data_list, year, ignore_species = NULL){
  
  all_values <- NULL
  
  species <- names(data_list)
  
  if(!is.null(ignore_species)){
    
    species <- species[!species %in% ignore_species]
    
  }
  
  count <- 0
  
  for(sp in species){
    # if the species has a value for this year add its data
    if(as.character(year) %in% row.names(data_list[[sp]])){
      all_values <- c(all_values, median(data_list[[sp]][as.character(year),]))
      count <- count + 1
    }
  }
  
  if(!is.null(all_values)){
    year_mean <- geomean(all_values)
    attr(year_mean, 'n_species') <- count
    return(year_mean)
  } else {
    return(NA)    
  }
}