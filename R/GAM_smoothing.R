#' GAM smoothing
#' 
#' This function uses a GAM to smooth a numeric vector
#'  
#' @param values A numeric vector that will be smoothed
#' @param years Optional, the years that \code{value} corresponds
#' to as a numeric vector.
#' @param plot If \code{TRUE} a plot is generated showing the smoothed
#' spline.
#' @param ... other arguements passed to \code{mgcv::gam}
#' @return The smoothed/predicted values are returned as well as the 
#' model object used for making predictions.
#' @import mgcv
#' @export

GAM_smoothing <- function(values, years = 1:length(values), plot = FALSE, ...){
  
  if(!is.numeric(years) | !is.vector(years)) 
    stop('years must be a numeric vector')
  if(!is.numeric(values) | !is.vector(values)) 
    stop('values must be a numeric vector')
  if(length(years) != length(values)) 
    stop('years and values must be the same length')
  
  mod_gam1 <- gam(values ~ s(years), ...)
  
  predicted_values <- predict(mod_gam1, newdata = data.frame(years = years))
  
  attr(predicted_values, which = 'model') <- mod_gam1
  
  if(plot == TRUE){
    
    plot(mod_gam1, residuals = T, pch = 19, cex = 0.9,
         scheme = 1, col = 'black', shade = T, shade.col = 'gray90',
         xlab = 'Year', ylab = 'Predicted Value')
    
  }
  
  dimnames(predicted_values)[[1]] <- years
  return(predicted_values)
  
}