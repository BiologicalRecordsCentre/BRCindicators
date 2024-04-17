lambda_calc <- function(x){
  
  # simpler than before. It assumes data are on the Log scale
  # takes just a vector of occupancy scores, not a 3D array
  LogLambda <- rep(NA, length(x))
  
  # All values should be NA then we start after the first year
  # with data. The first year with data is set to 0
  yearOne <- min(which(!is.na(x)))
  LogLambda[yearOne] <- 0 
  
  for (t in (yearOne + 1):length(x)){
    
    relyrs <- (!is.na(x[1:(t - 1)]))
    
    if (any(relyrs)){
      
      mrry <- max(which(relyrs)) # most recent reliable year
      LogLambda[(mrry + 1):t]  <- (x[t] - x[mrry])/(t - mrry)
      
    } else {
      
      LogLambda[t] <- NA
      
    }
  }
  
  return(LogLambda)
  
}