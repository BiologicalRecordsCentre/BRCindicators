#' Rescale species values for indicator
#' 
#' This function takes in a dataframe of multi-species data and rescales them
#' so that the value in the starting year is the same. This function accounts
#' for species that have no data at the beginning of the period or those that
#' have no data at teh end of the period.
#' 
#' @param Data A matrix, the first column, named year, gives the year. Subsequent named
#'        columns give the species. Values in the table are the yearly values
#'        to be rescaled.
#' @param index The index value for the first year, defaults to 100.
#' @param max The upper limit allowed for scaled values. Values greater than 
#'        this will be set equal to \code{max}.
#' @param min The upper limit allowed for scaled values. Values greater than 
#'        this will be set equal to \code{min}.
#' @param weights an optional dataframe of weights for the species when calculating the geometric mean.
#'        Must contain species names in first column and a weight in the second column. 
#'        Returns an error if not all species are included.      
#' @return A matrix. In each row we have the year, the species and the
#'         scaled value. There is also an additional column, 'geomean' giving
#'         the geometric mean for each year.
#' @export

rescale_species <-  function(Data, index = 100, max = 10000,
                             min = 1, weights=NULL){
  
  if(!is.null(weights)) {
    # first check the species names match up & weights don't have negative numbers or NAs
    if(!all(dimnames(Data)[[2]][-1] == weights$Species)) 
      stop("names in weights don't match species names in data")
    else if(any(weights$weight < 0) | any(is.na(weights$weight)))
      stop("weights include negative or NA values")
    else wt <- weights$weight/mean(weights$weight[weights$weight>0]) # scale to mean of 1, ignoring zeros
  } else wt <- rep(1, dim(Data)[2]-1)
  
  # weighted geomean function
  geomean <- function(x) exp(mean(log(x)*wt, na.rm = T))
  
  if(!inherits(Data, 'matrix')){
    if(inherits(Data, 'data.frame')){
      Data <- data.matrix(Data)
    } else {
      stop('Data must be a matrix')  
    }
  } 
  
  # Get the multipliers needed to achieve the index value. 
  #multipliers <- as.numeric(index / Data[1,2:ncol(Data)])
  # change to scaling by 1, so that starting value is zero on log scale. This makes it easy to weight.
  # This is the value in the first NA year
  multipliers <- 1 / apply(Data[,2:ncol(Data)], 2, function(x) x[!is.na(x)][1])
  
  # Apply these multipliers to scale the data
  #indicator_scaled <- t(t(Data[,2:ncol(Data)]) * multipliers) # returns only a vector, not 2D data
  indicator_scaled <- t(apply(Data[,2:ncol(Data)], 1, function(x) x * multipliers))

  # Make values over max == max, and < min == min
  indicator_scaled[indicator_scaled < min & !is.na(indicator_scaled)] <- min
  indicator_scaled[indicator_scaled > max & !is.na(indicator_scaled)] <- max
  
  # Species with missing values are now in here with all NA rows
  # We want to calculate the geomean from this data first
  geomean_vals <- apply(X = indicator_scaled, MARGIN = 1, FUN = geomean)
  
  indicator_scaled <- cbind(indicator_scaled, geomean_vals)
  colnames(indicator_scaled)[ncol(indicator_scaled)] <- 'geomean'
  
  # Now we need to add back in the NA species and scale them all so that
  # their first value is the same as the geomean for that year
  # These have NAs at the beggining
  NAtop <- colnames(Data)[is.na(Data[1,])]

  # These need to be sorted into the order in which they enter the dataset for
  # this method to work (else assigning geomean values goes wrong)
  firstYear <- function(x) min(which(!is.na(x)))
  # re-order if there is more than one
  if(length(NAtop) > 1) NAtop <- names(sort(apply(X = Data[,NAtop], MARGIN = 2, FUN = firstYear)))
  
  if(length(NAtop) > 0){
    # Deal with ones at the begining first
    for(i in 1:length(NAtop)){# Create a column of T/F if NA or not
      
      # Create a temporary dataframe for this species
      temp_col <- data.frame(species = !is.na(Data[,NAtop[i]]), row = 1:nrow(indicator_scaled))
      
      # index of first good year
      first_year <- min(temp_col$row[temp_col$species])
      
      # get geomean for this year 
      temp_gm <- indicator_scaled[first_year,'geomean']
      
      # Calculate multiplier needed
      multi <- temp_gm / Data[first_year, NAtop[i]]
      
      # Apply this...
      d <- Data[,NAtop[i]]*multi
      
      # ...make sure non are under 1 or over 10000...
      d[d < min & !is.na(d)] <- min
      d[d > max & !is.na(d)] <- max
      
      # ...and put it in the table
      indicator_scaled[,NAtop[i]] <- d
      
      # Recalculate the geomean before the next species
      indicator_scaled[,'geomean'] <- apply(X = indicator_scaled[,!colnames(indicator_scaled) %in% 'geomean'], MARGIN = 1, FUN = geomean)
      
    }
  }
  
  # If a species drops out hold at last value
  fillTailNAs <- function(x){
    
    # Get trues and falses for locations of NAs
    na_true_false <- is.na(x)
    
    # Get the position of all falses
    na_position <- grep(FALSE, na_true_false)
    
    # If the max false is hte last year dont do anything...
    if(!max(na_position) == length(x)){
      
      # else give all the last years the value at the last false
      x[(max(na_position)+1):length(x)] <- x[max(na_position)]
    }
    return(x)    
  }

  # apply tail function
  temp_indicator_scaled <- apply(X = indicator_scaled[,-ncol(indicator_scaled)],
                            MARGIN = 2, FUN = fillTailNAs)

  
  # Recalculate geomean
  indicator <- index * apply(X = temp_indicator_scaled, MARGIN = 1, FUN = geomean)
  
  # bind to the species indicies
  indicator_scaled <- cbind(temp_indicator_scaled, indicator)
  
  # Format the columns
  colnames(indicator_scaled)[ncol(indicator_scaled)] <- "indicator"
  indicator_scaled <- cbind(Data[, "year"], indicator_scaled)
  colnames(indicator_scaled)[1] <- "year"
  
  return(indicator_scaled)
}