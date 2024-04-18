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
    iter <- colnames(Occ[[1]])
  }
  
  # Build the array
  array_holder <- array(data = NA,
                        dim = c(nsp, nyr, length(iter)),
                        dimnames = list(species = names(Occ),
                                        years = as.character(minyr:maxyr),
                                        iterations = iter))
  
  # Fill the array
  for(i in 1:length(Occ)){ 
    
    array_holder[names(Occ[i]), row.names(Occ[[i]]), ] <- Occ[[i]]
    
  }
  
  return(array_holder)
  
}