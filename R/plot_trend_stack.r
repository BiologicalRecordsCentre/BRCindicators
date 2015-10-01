#' Plot species trend stack plots
#' 
#' @description  This function is used for reporting on species trends generated
#' in indicators. It takes a vector of categories and plots then in a stacked
#' column plot
#' 
#' @param species_change factor vector of change categories. This should have the
#' following levels "strong increase", "increase", "no change", "decrease",
#' "strong decrease".
#' @return Returns a data.frame of the data in the plot. A plot is sent to the device. 
#' @import ggplot2 RColorBrewer
#' @export
#' @examples 
#' ### Running from an array ####
#' set.seed(123)
#' # number of species
#' nsp = 50
#' 
#' # number of years
#' nyr = 40
#' 
#' #number of iterations
#' iter = 3
#' 
#' # Build a random set of data
#' myArray <- array(rnorm(n = nsp*nyr*iter,
#'                  mean = 0.5,
#'                  sd = 0.2),
#'                  c(nsp, nyr, iter))
#' 
#' # Ensure values are bounded by 0 and 1
#' myArray[myArray > 1] <- 1
#' myArray[myArray < 0] <- 0
#' 
#' # Run the lambda_interpolation method on this data                
#' myIndicator <- lambda_interpolation(myArray)
#' 
#' # Plot the trend stack
#' plot_trend_stack(myIndicator$species_change[,'category'])

plot_trend_stack <- function(species_change){
  
  # prep data to plot
  plot_data <- as.data.frame(table(species_change))
  plot_data$prop_spp <- plot_data$Freq/sum(plot_data$Freq)
  plot_data$species_change <- factor(plot_data$species_change, levels = c("strong increase", "increase", "no change", "decrease", "strong decrease"))
  plot_data$group <- "group"
  plot_data <- plot_data[rev(order(plot_data$species_change)),]
  
  # set manual colour brewer.pal
  cbPalette <- rev(brewer.pal(5, 'RdYlBu'))

  stackPlot <- ggplot(plot_data, aes(x = group, y = prop_spp, fill = species_change)) + 
    geom_bar(stat = "identity", width=.5) + 
    theme_classic() +
    scale_fill_manual(values=cbPalette) +
    ylab(paste("Proportion of species (n=", length(species_change),")", sep = '')) +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
          legend.position = "right", legend.title = element_blank(),
          axis.title.y = element_text(size = 14, vjust=1)) +
    scale_x_discrete(name = "") 
  
  # print plot
  print(stackPlot)
  
  return(plot_data[,1:3])
  
}