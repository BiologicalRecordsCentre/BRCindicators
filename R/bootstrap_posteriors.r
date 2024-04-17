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