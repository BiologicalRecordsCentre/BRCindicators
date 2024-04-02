# create a function to read in the data we want from these .rdata files
read_posterior <- function(file, sample_size = NULL, region = NULL){
  
  load(file)
  
  # some old outputs dont have min year in which case make it == 1
  min_year <- ifelse(is.null(out$min_year), 1, out$min_year)
  
  # Extract the values for each iteration
  if(!is.null(region)){
    iteration_values <- t(out$BUGSoutput$sims.list[paste0("psi.fs.r_",region)][[1]])
  } else {
    iteration_values <- t(out$BUGSoutput$sims.list$psi.fs)
    warning('No region specified therefore defaulting to the full dataset (psi.fs)')
  }
  
  if(!is.null(sample_size)){
    # Sample is needed
    iteration_values <- iteration_values[, sample(ncol(iteration_values),
                                                  size = sample_size,
                                                  replace = FALSE)]
  }
   
  # Name the rows by year
  row.names(iteration_values) <- min_year:(min_year + nrow(iteration_values) - 1) 
  
  # Name the columns by iteration
  colnames(iteration_values) <- paste('i', 1:ncol(iteration_values), sep = '') 
  
  # Add rhat & sd
  if(!is.null(region)){
    sum_dat <- out$BUGSoutput$summary[grep(paste0("psi.fs.r_", region), row.names(out$BUGSoutput$summary)), c('sd', 'Rhat')]
  } else {
    sum_dat <- out$BUGSoutput$summary[grep("psi.fs\\[", row.names(out$BUGSoutput$summary)), c('sd', 'Rhat')]
  }
    
#    sum_dat <- out$BUGSoutput$summary[grepl("psi.fs",row.names(out$BUGSoutput$summary)),
#                                    c('sd', 'Rhat')]
  
  iteration_values <- cbind(iteration_values, sum_dat)
  
  return(iteration_values)
}