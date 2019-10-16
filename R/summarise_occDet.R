#' Summarise Occupancy model output
#' 
#' This function reads in and summarises all the .rdata files that are output
#' by the occDetModel in sparta.
#' 
#' @param input_dir The directory where all the .rdata files are located.
#'        This will be the same as the \code{output_dir} arguement given
#'        to occDetModel in sparta when creating the output.
#'        
#' @param region If the model output contains regional estimates specify the region
#'        name here to extract estimates from that region only (e.g."ENGLAND"). If NULL
#'        the estimate based on all sites will be used instead.
#'        
#' @param verbose If \code{TRUE} then progress is printed to screen.
#' 
#' @return A data.frame. In each row we have the year, the species and the
#'         mean of the posterior for the predicted proportion of sites occupied.
#' @importFrom reshape2 dcast   
#' @export

summarise_occDet <-  function(input_dir, region = NULL, verbose = TRUE){
    
    library(reshape2)

    # get files from the input directory
    files <- list.files(path = paste(input_dir), ignore.case = TRUE, pattern = '\\.rdata$') # list of the files to loop through
    
    # sense check these file names
    if(length(files) == 0) stop('No .rdata files found in ', input_dir)
    if(length(files) < length(list.files(path = input_dir))) warning('Not all files in ', input_dir, ' are .rdata files, other file types have been ignored')
    
    loadRData <- function(fileName){
      #loads an RData file, and returns it
      load(fileName)
      if(length(ls()[ls() != "fileName"]) > 1){
        stop('The rdata file(s) contain more than one object. They should contain a single object, the JAGS model output from the sparta package')
      }
      get(ls()[ls() != "fileName"])
    }
    
    # create a function to read in the data we want from these .rdata files
    read_bayes <- function(file, in_region = region){
      
      out <- loadRData(file) 
      
      if(!(is.null(in_region))){
        if(!(in_region %in% out$regions)) stop("Requested region is not listed in model output")
      }
      
      out <- loadRData(file) 

      # some old outputs dont have min year in which case make it == 1
      min_year <- ifelse(is.null(out$min_year), 1, out$min_year)
      #Get the summary output for the rows and columns that we are interested in
      temp_out <- as.data.frame(out$BUGSoutput$summary)
      if(is.null(in_region)){
        rows <- grep("^(psi.fs[^.r])", row.names(temp_out))
        summary_out <- data.frame(year = (min_year - 1) + as.numeric(gsub("psi.fs", "", gsub("\\[|\\]","", row.names(temp_out)[rows]))),
                                  mean_proportion_sites = temp_out[rows, c("mean")],
                                  species_name = gsub('.rdata', '', basename(file)))
      }else{
        rows <- grep(paste0("psi.fs.r_",in_region), row.names(temp_out))
        summary_out <- data.frame(year = (min_year - 1) + as.numeric(gsub(paste0("psi.fs.r_",in_region), "", gsub("\\[|\\]","", row.names(temp_out)[rows]))),
                                mean_proportion_sites = temp_out[rows, c("mean")],
                                species_name = gsub('.rdata', '', basename(file)))
      }
      ### Replace the name bit with something better ###
      return(summary_out)
    }
    
    if(verbose) cat('Loading data...')
    # Use lapply to run this function on all files
    list_summaries <- lapply(file.path(input_dir, files), read_bayes)
    if(verbose) cat('done\n')
    
    # Unlist these and bind them together
    spp_data <- do.call(rbind, list_summaries)

    # Cast this into a wide format (this is a format people may have there data in)
    spp_data <- as.matrix(dcast(data = spp_data,
                                formula = year ~ species_name, value.var = 'mean_proportion_sites'))
    
    return(spp_data)
}