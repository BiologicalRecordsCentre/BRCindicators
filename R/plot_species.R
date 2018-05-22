#' @import ggplot2
#' @export

plot_species <- function(spindex, upper, lower, data_in, year = 1:ncol(spindex), main=''){
  
  spindex <- melt(spindex, value.name='spindex')
  upper <- melt(upper, value.name='upper')
  lower <- melt(lower, value.name='lower')
  combo <- merge(merge(spindex, upper), lower)
  names(combo)[1:2] <- c('species', 'year')
  combo <- merge(combo, data_in)
#  combo <- subset(combo, year > 1)# optional, since spindex is not estimated in year 1
    
  # Plot
  plot <- ggplot(combo, aes(x = year)) + 
    geom_line(aes(y = spindex), size = 1) +
    scale_y_continuous(limits = c(min(combo$lower), max(combo$upper))) +
    #geom_hline(yintercept = index, linetype = "dashed") +
    labs(x = "Year", y = "Index") + 
    #geom_point(aes(y = spindex), size = 3) +
    theme(text = element_text(size = 15)) +
    geom_ribbon(data = combo, aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_wrap(~species) +
    geom_point(aes(y = index), col = 'red') +
    ggtitle(main)

  print(plot)
  
}
