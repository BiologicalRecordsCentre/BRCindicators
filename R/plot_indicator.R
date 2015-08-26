#' Plot indicator
#' 
#' This function plots the indicator using ggplot2.
#' 
#' @param indicator A numeric vector, one value for each year
#' @param CIs A matrix with two columns, the first column is the lower bound
#'        and the second is the upper bound.
#' @param smoothed_line A numeric vector giving the smoothed values to plot
#'        see \code{\link{GAM_smoothing}}
#' @param year A numeric vector as long as \code{indicator} giving the years
#' @param index Numeric, the index value, by default it is taken to be the 
#'        first value of \code{indicator}
#' @param main Title given to the plot
#' @return The plot object is returned
#' @import ggplot2
#' @export

plot_indicator <-  function(indicator, CIs, year = 1:length(indicator), index = indicator[1],
                            smoothed_line = NULL, main = ''){
  
    combo <- as.data.frame(cbind(year, indicator, CIs))
    colnames(combo) <- c('year','indicator','lower','upper')
    
    if(!is.null(smoothed_line)){
      if(length(indicator) != length(smoothed_line)) stop('indicator and smoothed_line must be the same length')
    }
    
    if(!is.null(smoothed_line)){
      combo$smoothed <- smoothed_line
    }
      
    # Plot
    plot <- ggplot(combo, aes(x = year)) + 
      geom_line(aes(y = indicator), size = 1) +
      scale_y_continuous(limits = c(0, max(combo$upper))) +
      geom_hline(yintercept = index, linetype = "dashed") +
      labs(x = "Year", y = "Scaled value") + 
      geom_point(aes(y = indicator), size = 3) +
      theme(text = element_text(size = 15)) +
      geom_ribbon(data = combo, aes(ymin = lower, ymax = upper), alpha = 0.2) +
      ggtitle(main)
    
    if(!is.null(smoothed_line)){
      
      plot <- plot +
              geom_line(aes(y = smoothed, colour = 'red'), size = 1) +
              theme(legend.position = "none")
    }
    
    print(plot)
    
    #return(NULL)
}
