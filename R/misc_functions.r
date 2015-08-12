geomean <- function(x) exp(mean(log(x), na.rm = T))

# create a function to read in the data we want from these .rdata files
read_posterior <- function(file){
  
  load(file)
  
  # some old outputs dont have min year in which case make it == 1
  min_year <- ifelse(is.null(out$min_year), 1, out$min_year)
  
  # Extract the values for each iteration
  iteration_values <- t(out$BUGSoutput$sims.list$psi.fs)
  
  # Name the rows by year
  row.names(iteration_values) <- min_year:(min_year + nrow(iteration_values) - 1) 
  
  return(iteration_values)
}

remove_bad_years <-  function(species_name, subset_table, list_summaries){
  
  # Good years for this species
  good_years <- names(subset_table[species_name, subset_table[species_name, ] == TRUE])
  
  # Data for this species
  data <- list_summaries[[species_name]]
  
  # Now subset
  new_data <- data[good_years, ]
  
  return(new_data)
  
}  

# Create a function that will calculate the geomean from the list
# of species for a given year
list_geomean <- function(data_list, year, ignore_species = NULL,
                         CI_min = 0.025, CI_max = 0.975){
  
  all_values <- NULL
  
  species <- names(data_list)
  
  if(!is.null(ignore_species)){
    
    species <- species[!species %in% ignore_species]
    
  }
  
  count <- 0
  
  for(sp in species){
    # if the species has a value for this year add its data
    if(as.character(year) %in% row.names(data_list[[sp]])){
      all_values <- c(all_values, data_list[[sp]][as.character(year),])
      count <- count + 1
    }
  }
  
  year_mean <- geomean(all_values)
  attr(year_mean, 'n_species') <- count
  attr(year_mean, 'CI_min') <- quantile(x = all_values, probs = CI_min)
  attr(year_mean, 'CI_max') <- quantile(x = all_values, probs = CI_max)
  
  return(year_mean)
  
}

# Create a function to apply multiplier to years after a species leaves
# year - is the first of the years to be adjusted ie the first year in
# which the species is absent
leaving_multiplier <- function(data_list, year, multiplier){
  
  for(sp in names(data_list)){
    
    data_list[[sp]][as.numeric(row.names(data_list[[sp]])) >= year,] <- data_list[[sp]][as.numeric(row.names(data_list[[sp]])) >= year,] * multiplier
    
  }
  
  return(data_list)
  
}

# Create a function that rescales a set of data so that the first year geomean
# is equal to the geomean given
entering_multiplier <- function(data_list, gm){
  
  for(sp in names(data_list)){
    
    multiplier <- gm / geomean(data_list[[sp]][1,])
    
    data_list[[sp]] <- data_list[[sp]] * multiplier
    
  }
  
  return(data_list)
  
}

# A function to take in a list of tables and cap at the max and min index values
cap_index <- function(data_list, max = 10000, min = 1){
  
  for(sp in names(data_list)){
    
    data_list[[sp]][data_list[[sp]] > max] <- max
    data_list[[sp]][data_list[[sp]] < min] <- min
    
  }
  
  return(data_list)
  
}