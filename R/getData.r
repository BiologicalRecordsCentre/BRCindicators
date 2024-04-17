getData <- function(input, sample_size = NULL, region = NULL){
  
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
                  simplify = 'array', sample_size = sample_size, region = region)
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