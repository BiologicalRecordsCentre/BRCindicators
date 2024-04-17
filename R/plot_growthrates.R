#' @import ggplot2
plot_growthrates <- function(logLambda, CIs, year = 1:length(logLambda), main=''){
  combo <- as.data.frame(cbind(year, logLambda, CIs))
  
  # Plot
  plot <- ggplot(combo, aes(x = year)) + 
    geom_line(aes(y = logLambda), size = 1) +
    scale_y_continuous(limits = c(min(combo$lower), max(combo$upper))) +
    #geom_hline(yintercept = index, linetype = "dashed") +
    labs(x = "Year", y = "LogLambda") + 
    geom_point(aes(y = logLambda), size = 3) +
    theme(text = element_text(size = 15)) +
    geom_ribbon(data = combo, aes(ymin = lower, ymax = upper), alpha = 0.2) +
    ggtitle(main)
  
  print(plot)
  return(NULL)
}

#' @import ggplot2
plot_spgrowthrates <- function(logLambda, upper, lower, year = 1:ncol(logLambda), main=''){
  
  logLambda <- melt(logLambda, value.name='logLambda')
  upper <- melt(upper, value.name='upper')
  lower <- melt(lower, value.name='lower')
  combo <- merge(merge(logLambda, upper), lower)
  names(combo)[1:2] <- c('species', 'year')

  # Plot
  plot <- ggplot(combo, aes(x = year)) + 
    geom_line(aes(y = logLambda), size = 1) +
    scale_y_continuous(limits = c(min(combo$lower), max(combo$upper))) +
    #geom_hline(yintercept = index, linetype = "dashed") +
    labs(x = "Year", y = "Growth rate") + 
    geom_point(aes(y = logLambda), size = 3) +
    theme(text = element_text(size = 15)) +
    geom_ribbon(data = combo, aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_wrap(~species) +
    ggtitle(main)
  
  print(plot)
  return(NULL)
}

#' @import ggplot2
wrap_plot <- function(model, name="", data_in=NULL){
  
  plot_growthrates(logLambda = model$mean$logLambda,
                   CIs = with(model, data.frame(lower=q2.5$logLambda, upper=q97.5$logLambda)),
                   main = name)
  plot_spgrowthrates(logLambda = model$mean$spgrowth,
                     lower = model$q2.5$spgrowth,
                     upper = model$q97.5$spgrowth,
                     main = name)
  plot_species(spindex = model$mean$spindex,
               lower = model$q2.5$spindex,
               upper = model$q97.5$spindex,
               data_in = data_in,
               main = name)
  return(NULL)
}