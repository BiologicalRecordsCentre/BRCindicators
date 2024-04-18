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
    
    Occ <- Occ[ , seq(from = min(year_range), to = max(year_range)), , drop = FALSE]
    
  } else {
    
    Occ <- Occ[ , as.character(seq(from = min(year_range), to = max(year_range))), , drop = FALSE]
    
  }
  
}