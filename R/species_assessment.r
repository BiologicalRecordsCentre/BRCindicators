species_assessment <- function(dat,
                               method = "lambda",
                               start_year = NULL,
                               end_year = NULL,
                               species_stat = 'mean',
                               plot = FALSE) {
  # Sense checks
  if(!method %in% c("lambda", "bma")) stop("Method must be one of 'lambda' or 'bma'")
  if(!species_stat %in% c('mean', 'median')) stop("species_stat must be either 'mean' or 'median'")
  
  if(method == "lambda") {
    
    # lambda method
    LogLambda <- dat
    
    # If we are subsetting
    if(!is.null(start_year) | !is.null(end_year)) {
      
      # If the assessment is only over a subset of the years do the subsetting first
      if(is.null(dimnames(LogLambda)[[2]])) {
        
        if(is.null(start_year)) start_year <- 1
        if(is.null(end_year)) end_year <- dim(LogLambda)[2]
        
        if(end_year > dim(LogLambda)[2]) {
          
          stop(paste('Years are un-named in your data and the specified end_year [',
                     end_year,
                     '] exceeds the number of years in the data [',
                     dim(LogLambda)[2],
                     ']', sep = ''))
          
        } else {
          
          LogLambda <- LogLambda[ , start_year:end_year, ]
          
        }
        
      } else {
        
        if(is.null(start_year)) start_year <- as.character(min(as.numeric(dimnames(LogLambda)[[2]])))
        if(is.null(end_year)) end_year <- as.character(max(as.numeric(dimnames(LogLambda)[[2]])))
        
        LogLambda <- LogLambda[ , as.character(start_year:end_year), ]
        
      }
      
    }
    
    # Calculate the average change across this time period
    if(species_stat == 'mean') spLogLamda <- rowMeans(apply(LogLambda, c(1,2), mean, na.rm = T), na.rm = T) # one value per species
    if(species_stat == 'median') spLogLamda <- apply(apply(LogLambda, c(1,2), mean, na.rm = T), 1, FUN = median, na.rm = T)
    
    # Remove NAs
    spLogLamda <- spLogLamda[!is.na(spLogLamda)]
    
    # this last value is simple, but conflates uncertainty with interannual variation
    # of one value per species, convert to a percentage change per year
    sp_pcpy <- 100 * (exp(spLogLamda) - 1)
    
    # Assign to cats
    sp_cat <- cut(sp_pcpy, 
                  breaks = c(-Inf,-2.73,-1.14,1.16,2.81,Inf),
                  labels = c('strong decrease','decrease', 'no change', 'increase','strong increase'),
                  ordered = T)
    
    # build DD
    sp_change <- data.frame(percent_change_year = sp_pcpy, category = sp_cat)
    
    # Plot is desired
    if(plot) plot_trend_stack(species_change = sp_change$category)
    
    return(sp_change)
    
  } else {
    
    # bma method
    bma_df <- dat
    
    # mean growth rate between years per species
    spgrowth <- attr(bma_df, "model")$mean$spgrowth
    
    yr_df <- data.frame(startyr = seq(min(bma_df$year), (max(bma_df$year) - 1), 1), 
                        endyr = seq(min(bma_df$year) + 1, max(bma_df$year), 1),
                        yr_col = seq(1, ncol(spgrowth), 1))
    
    if(is.null(start_year)) start_year <- 1 
    # convert years from CE to numeric (e.g., 1970 to 1)
    else start_year <- yr_df[yr_df$startyr == start_year,]$yr_col
    
    if(is.null(end_year)) end_year <- ncol(spgrowth)
    # convert years from CE to numeric (e.g., 2016 to 46)
    else end_year <- yr_df[yr_df$endyr == end_year,]$yr_col
    
    # subset to focal years
    spgrowth_sub <- spgrowth[ , start_year:end_year]
    
    # average growth rate across years per species
    if(species_stat == 'mean') spgrowth_av <- rowMeans(spgrowth_sub)
    if(species_stat == 'median') spgrowth_av <- apply(spgrowth_sub, 1, median, na.rm = TRUE)
    
    # back-transform from log-scale
    spg_pcpy <- 100*(exp(spgrowth_av)-1)
    
    # assign to categories
    spg_cat <- cut(spg_pcpy, 
                   breaks = c(-Inf,-2.73,-1.14,1.16,2.81,Inf),
                   labels = c('Strong decrease','Decrease', 'No change', 'Increase','Strong increase'),
                   ordered = TRUE)
    
    # build dataframe
    sp_change <- data.frame(percent_change_year = spg_pcpy, category = spg_cat)
    
    # Plot is desired
    if(plot) plot_trend_stack(species_change = sp_change$category)
    
    return(sp_change)
    
  }
  
}