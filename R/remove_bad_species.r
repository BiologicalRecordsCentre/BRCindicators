remove_bad_species <- function(Occ, threshold_sd, threshold_yrs, threshold_Rhat){
  
  # If we have sd and Rhat from the input file use those
  if('sd' %in% dimnames(Occ)[[3]] & 'Rhat' %in% dimnames(Occ)[[3]]){
    
    reliable <- Occ[ , ,'sd', drop = FALSE] < threshold_sd & Occ[ , ,'Rhat', drop = FALSE] < threshold_Rhat
    
    # drop the third dimension (we want to keep the 3rd dimension even if it has
    # length 1, so this is a little hacky)
    reliable <- apply(reliable, c(1,2), mean)
    reliable[reliable == 0] <- NA
    
    # set the posteriors where sd or Rhat not met to NA 
    # remove the sd and Rhat columns
    Occ <- Occ[ , , !dimnames(Occ)[[3]] %in% c('sd','Rhat'), drop = FALSE]

  # else if it is an array calc sd on the fly
  } else {
    
    # Calculate standard deviations
    sds <- apply(Occ, c(1,2), sd, na.rm = TRUE) 
    
    # first define reliable estimates as those with standard deviations lower than the threshold
    reliable <- sds < threshold_sd
    
    # convert the FALSE elements to NA (for the maths to work below)
    reliable[!reliable] <- NA 

  }
  
  OccRel <- sapply(1:dim(Occ)[3], function(j) Occ[ , ,j] * reliable,
                   simplify='array')
  
  # OccRel is now identical to Occ except that unreliable estimates have been changed to NA
  # now strip out species with fewer reliable years than the desired number
  OccRel <- OccRel[rowSums(reliable, na.rm=T) >= threshold_yrs, , , drop = FALSE]

  if(dim(OccRel)[1] == 0) stop('None of your species meet the thresholds')
  
  # Save the good years table - for the species that pass
  # the thresholds - as an attribute of the data returned
  attr(OccRel, 'good_years') <- reliable[rowSums(reliable, na.rm = T) >= threshold_yrs, ]
  return(OccRel)
  
}