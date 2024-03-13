bootstrap_indicator_fake_data <- matrix(runif(50 * length(letters), max = 100), 
                     nrow = 50, 
                     ncol = length(letters))

# Assign the same column names
colnames(bootstrap_indicator_fake_data) <- letters


bootstrap_indicator_objects <-  function(Data, iterations = 10000, CI_limits = c(0.025, 0.975),
                                 verbose = TRUE){

    # To capture function objects
    output = list()

  if(!is.matrix(Data)){
    stop("the Data parameter must be a matrix object.")
  }

  if(!all(is.numeric(Data))){
   stop("Matrix values must all be numeric.") 
  }
  
  ### bootstrap to estimate 95% CI around the geometric mean ###
  nSpecies <- ncol(Data)
  bootstrap_values <- matrix(data = NA, nrow = nrow(Data), ncol = iterations)
  bootstrap_data <- as.data.frame(Data)

  pb <- txtProgressBar(min = 0, max = iterations, style = 3)

  # Run bootstrapping
  CIData <- sapply(1:iterations, simplify = TRUE, function(iteration) {

  if(verbose){setTxtProgressBar(pb, iteration)}

  # Randomly sample species 
  samp <- sample(bootstrap_data, size = ncol(bootstrap_data), replace = TRUE)
 
  apply(samp, 1, function(x){
    exp(mean(log(x), na.rm = T))
    })
})

output$CIData = CIData

  close(pb)

  # Extract the 2.5 97.5% CI around the geometric mean (from the bootstrapped resamples)
  CIs <- as.data.frame(t(apply(X = CIData, MARGIN = 1, FUN = quantile,
                               probs = CI_limits, na.rm = TRUE)))
  names(CIs) <- c(paste('quant_', gsub('0\\.', '', as.character(CI_limits[1])), sep = ''),
                  paste('quant_', gsub('0\\.', '', as.character(CI_limits[2])), sep = ''))

    output$CI_matrix = as.matrix(CIs)

  return(output)
}