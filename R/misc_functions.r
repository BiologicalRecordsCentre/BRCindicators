geomean <- function(x) exp(mean(log(x), na.rm = T))

# create a function to read in the data we want from these .rdata files
read_posterior <- function(file, sample_size = NULL){
  
  load(file)
  
  # some old outputs dont have min year in which case make it == 1
  min_year <- ifelse(is.null(out$min_year), 1, out$min_year)
  
  # Extract the values for each iteration
  iteration_values <- t(out$BUGSoutput$sims.list$psi.fs)
  
  if(!is.null(sample_size)){
    # Sample is needed
    iteration_values <- iteration_values[, sample(ncol(iteration_values),
                                                  size = sample_size,
                                                  replace = FALSE)]
  }
   
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
    
    multiplier <- gm / median(data_list[[sp]][1,])
    
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

# A function to bootstrap across a rescaled list to get confidence intervals
bootstrap_posteriors <- function(rescaled_list, iterations = 10000, years){
  
  pb <- txtProgressBar()
  cat(paste('Bootstrapping -', iterations, 'iterations', '\n'))
  
  rep <- function(rescaled_list, progress, years){
    
    setTxtProgressBar(pb, progress)
    
    sampled <- sample(rescaled_list, length(rescaled_list), replace = TRUE)
    
    return(sapply(years, FUN = function(x) list_geomean(sampled, year = x)))
    
  }
  
  bootstraps <- sapply(seq(from = 0, to = 1, length.out = iterations),
                       FUN = function(x) rep(rescaled_list, x, years))
  bootstrap_CIs <- t(apply(bootstraps, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  
  colnames(bootstrap_CIs) <- paste('bootstrap', colnames(bootstrap_CIs), sep = '_')
  
  cat('\n')
  
  return(bootstrap_CIs)
  
}

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


remove_bad_species <- function(Occ, threshold_sd, threshold_yrs){
  
  # Calculate standard deviations
  sds <- apply(Occ, c(1,2), sd) 
  
  # first define reliable estimates as those with standard deviations lower than the threshold
  reliable <- sds < threshold_sd
  
  # convert the FALSE elements to NA (for the maths to work below)
  reliable[!reliable] <- NA 
  
  # set the posteriors where sd >  0.2 to NA 
  OccRel <- sapply(1:dim(Occ)[3], function(j) Occ[,,j] * reliable, simplify='array')
  
  # OccRel is now identical to Occ except that unreliable estimates have been changed to NA
  # now strip out species with fewer reliable years than the desired number
  OccRel <- OccRel[rowSums(reliable, na.rm=T) >= threshold_yrs,,]

  # Save the good years table - for the species that pass
  # the thresholds - as an attribute of the data returned
  attr(OccRel, 'good_years') <- reliable[rowSums(reliable, na.rm = T) >= threshold_yrs, ]
  return(OccRel)
  
}

lambda_calc <- function(Occ, logOdds = TRUE){
  
  # first generate an empty array
  LogLambda <- array(NA, dim = dim(Occ), dimnames = dimnames(Occ))
  
  #loop over every year-species combo
  for(i in 1:dim(Occ)[1]){ 
    
    for(t in 2:dim(Occ)[2]){
      
      # which years in the past have reliable estimates
      relyrs <- (!is.na(Occ[i,1:(t-1),1]))
      
      if(any(relyrs)){
        
        # find the most recent year with reliable estimates (mrry)    
        mrry <- max(which(relyrs))
        
        # the compound interest over this period is given by the exponent of the difference in log odds between years
        if(logOdds){
          comp <- exp(Occ[i,t,] - Occ[i,mrry,]) # Odds ratio = proportional change in odds over the number of years
        } else {
          comp <- Occ[i,t,] / Occ[i,mrry,] # ratio
        }
        
        # the interpolated annual lamdba is the Nth root of comp, where N is the number of years between adjacent reliable estimates
        ann_change <- comp ^ (1/(t-mrry))
        
        # now interpolate all the years between mrry & t
        # NB if the gap is 0 years (i.e. mrry = t-1 then there is no interpolation)
        # if the gap is >0 years then the intervening years have already been set as NA
        LogLambda[i,(mrry+1):t,] <- log(ann_change)
        
      } else {
        
        LogLambda[i,t,] <- NA
        
      }
    }
  }
  
  return(LogLambda)
  
}


getData <- function(input, sample_size = NULL){
  
  if(!is.null(sample_size)) if(!is.numeric(sample_size)) stop('sample_size must be numeric')
  
  
  if(class(input) == 'character'){
    # get files from the input directory
    files <- list.files(path = paste(input), ignore.case = TRUE, pattern = '\\.rdata$') # list of the files to loop through
    
    # sense check these file names
    if(length(files) == 0) stop('No .rdata files found in ', input)
    if(length(files) < length(list.files(path = input))) warning('Not all files in ', input, ' are .rdata files, other file types have been ignored')
    
    cat('Loading data...')
    # Use lapply to run this function on all files
    org <- getwd()
    setwd(input)
    Occ <- sapply(files, read_posterior,
                  simplify = 'array', sample_size = sample_size)
    setwd(org)
    cat('done\n')
    
    # This does not come back as an array if the individual species
    # results do not have the same dimensions, this happens when the 
    # data come from different runs of the occupancy models
    if(class(Occ) != 'array'){
      
      Occ <- list_to_array(Occ)
      
    } else {
      
      # reorder dimensions into a more logical order
      # Species - Year - Iteration
      Occ <- aperm(Occ, c(3,1,2))
      
    }
    
    return(Occ)
  
  } else if(class(input) == "array"){
    
    if(!is.null(sample_size)){
      
      input <- input[ , , sample(ncol(input),
                                 size = sample_size,
                                 replace = FALSE)] 
      
    }
    
    return(input)
    
  } else {
    
    stop('Input should be either a file path or an array')
    
  }
  
}

list_to_array <- function(Occ){
  
  # Calculate dimensions
  nsp <- length(Occ)
  maxyr <- max(unlist(lapply(Occ, FUN = function(x) as.numeric(row.names(x)))))
  minyr <- min(unlist(lapply(Occ, FUN = function(x) as.numeric(row.names(x)))))
  nyr <- (maxyr - minyr) + 1
  
  # Throw an error if the number of iterations is not the same
  if(min(unlist(lapply(Occ, ncol))) != max(unlist(lapply(Occ, ncol)))){
    stop('Input data must have the same number of iterations')
  } else {
    iter <- as.numeric(unlist(lapply(Occ, ncol))[1])
  }
  
  # Build the array
  array_holder <- array(data = NA,
                        dim = c(nsp, nyr, iter),
                        dimnames = list(species = names(Occ),
                                        years = as.character(minyr:maxyr),
                                        iterations = as.character(1:iter)))
  
  # Fill the array
  for(i in 1:length(Occ)){ 
    
    array_holder[names(Occ[i]), row.names(Occ[[i]]), ] <- Occ[[i]]
    
  }
  
  return(array_holder)
  
}

subset_years <- function(Occ, year_range = NULL){
  
  # The second dimension of Occ is years
  year_names <- dimnames(Occ)[[2]]
  
  if(is.null(year_names)){
    
    if(max(year_range) > dim(Occ)[2]){
      
      stop(paste('Years are un-named in your data and the maximum year in your year range [',
                 max(year_range),
                 '] exceeds the number of years in the data [',
                 dim(Occ)[2],
                 ']', sep = '')
      )
      
    }
    
    Occ <- Occ[ , seq(from = min(year_range), to = max(year_range)), ]
    
  } else {
    
    Occ <- Occ[ , as.character(seq(from = min(year_range), to = max(year_range))), ]
    
  }
  
}

species_assessment <- function(LogLambda,
                               start_year = NULL,
                               end_year = NULL,
                               plot = FALSE){
  
  # If we are subsetting
  if(!is.null(start_year) | !is.null(end_year)){
  
    # If the assessment is only over a subset of the years do the subsetting first
    if(is.null(dimnames(LogLambda)[[2]])){
      
      if(is.null(start_year)) start_year <- 1
      if(is.null(end_year)) end_year <- dim(LogLambda)[2]
    
      if(end_year > dim(LogLambda)[2]){
        
        stop(paste('Years are un-named in your data and the specified end_year [',
                   end_year,
                   '] exceeds the number of years in the data [',
                   dim(LogLambda)[2],
                   ']', sep = ''))
        
      } else {
        
        LogLambda <- LogLambda[ , start_year:end_year, ]
        
      }
      
    } else {
      
      if(is.null(start_year)) start_year <- as.character(min(as.numeric(dimnames(LogLambda)[[2]])))
      if(is.null(end_year)) end_year <- as.character(max(as.numeric(dimnames(LogLambda)[[2]])))
      
      LogLambda <- LogLambda[ , as.character(start_year:end_year), ]
      
    }
    
  }
  
  # Calculate the average change across this time period
  spLogLamda <- rowMeans(apply(LogLambda, c(1,2), mean, na.rm = T), na.rm = T) # one value per species
  
  # Remove NAs
  spLogLamda <- spLogLamda[!is.na(spLogLamda)]
  
  # this last value is simple, but conflates uncertainty with interannual variation
  # of one value per species, convert to a percentage change per year
  sp_pcpy <- 100 * (exp(spLogLamda) - 1)
  
  # Assign to cats
  sp_cat <- cut(sp_pcpy, 
                breaks = c(-Inf,-2.73,-1.14,1.16,2.81,Inf),
                labels = c('strong decrease','decrease', 'no change', 'increase','strong increase'),
                ordered = T)
  
  # build DD
  sp_change <- data.frame(percent_change_year = sp_pcpy, category = sp_cat)
  
  # Plot is desired
  if(plot) plot_trend_stack(species_change = sp_change$category)
  
  # Test change
  
  return(sp_change)
}

indicator_assessment <- function(summary_table,
                                 start_year = NULL,
                                 end_year = NULL){
  
  if(is.null(start_year)) start_year <- min(summary_table$year)
  if(is.null(end_year)) end_year <- max(summary_table$year)
  
  # Get the indicator value at the start
  start_ind <- summary_table$indicator[summary_table$year == start_year]
  
  # Get the confidence intervals for rhe end of the period
  end_CIs <- summary_table[summary_table$year == end_year, c('lower', 'upper')]
  
  # assess
  if(start_ind < end_CIs$upper & start_ind > end_CIs$lower){
    assessment <- 'stable'
  } else if(start_ind < end_CIs$lower){
    assessment <- 'decreasing'
  } else if(start_ind > end_CIs$upper){
    assessment <- 'increasing'
  }
  
  # Return the assessment
  assessment <- data.frame(start_index = start_ind,
                           end_lower = end_CIs$lower,
                           end_upper = end_CIs$upper,
                           assessment = assessment)
  
}