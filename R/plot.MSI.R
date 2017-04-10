#' A function to plots msi results
#' 
#' Plots a summary of the MSI analysis
#' 
#' @param x object of class MSI
#' @param title Character, the title for the plot
#' @param \dots currently ignored
#'
#' @name plot.MSI
#' @method plot MSI
#' @import ggplot2
#' @export

plot.MSI <- function(x, title = 'MSI results', ...) {
 
  ggplot(data = x$results, aes(x = year, y = MSI)) +
    geom_ribbon(aes(ymin = x$results$lower_CL_trend, 
                    ymax = x$results$upper_CL_trend),
                colour = 'grey', fill = 'blue', alpha = 0.4) +
    geom_line(aes(y = x$results$Trend), colour = 'blue') +
    geom_point() +
    ylim(c(0, max(x$results$MSI) + (max(x$results$MSI)/100)*10)) +
    geom_errorbar(aes(ymin = x$results$MSI - (x$results$sd_MSI*1.96),
                      ymax = x$results$MSI + (x$results$sd_MSI*1.96))) +
    ggtitle(title) +
    annotate('text', y = 20, x = min(x$results$year), 
             label = paste(paste0(x$trends$Measure[1], ':'),
                           x$trends$significance[1]), hjust = 0) +
    annotate('text', y = 10, x = min(x$results$year), 
             label = paste(paste0(x$trends$Measure[3], ':'),
                           x$trends$significance[3]), hjust = 0)
   
}