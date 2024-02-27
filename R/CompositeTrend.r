#' CompositeTrend function
#' 
#' @description This function can be used to produce composite metrics of change 
#' (indicators), whilst propagating uncertainty in the individual 
#' species trend estimates through to the final composite trend 
#' metric. This function takes in a dataframe of sampled annual
#' occupancy estimates across multiple species and returns a 
#' a single composite trend metric with uncertainty. This approach
#' is only suitable for species without missing years.
#'
#' @param indata The file path to an .rdata file containing a dataframe
#'		  (called samp_post - GP to change this) which contains 
#'		  year columns (prefixed with "X", e.g "X1985"), a species
#'		  column ("spp"), and an iteration identifier ("iter"). The 
#'		  year columns contain the annual occupancy estimates for 
#' 		  the species-year-iteration combination in question.
#' @param output_path The location where the outputs should be saved.
#' @param trend_choice The approach used to combine the individual species 
#'		  estimates into a single composite trend. See details.
#' @param group_name The name of the species group we are running, used for
#'  	  naming output files.
#' @param save_iterations Do we want to save the composite trend estimates
#'		  for each individual iteration, these are generally used for 
#'		  estimating temporal trends with uncertainty.
#' @param TrendScale Traditionally some indicators are scaled so the first 
#'		  year is set to a given number, 100 in the case of the UK 
#'		  biodiversity indicators. This value can be chosen here, with no
#'		  scaling as the default.
#' @param plot_output plot the resulting composite indicator: TRUE or FALSE.
#  @details There are a number of model to choose from:
#' \itemize{
#'  \item{\code{"arithmetic_logit_occ"}}{ - The raw occupancy values are 
#'  converted to the log odds scale (using the logit function). The
#'  arithmetic mean across species is used to create a composite trend for each
#'  iteration. These means are then converted back to the odds scale (exp)  }
#'  \item{\code{"geometric_raw_occ"}}{ - Take the geometric mean across species 
#'  raw occupancy estimates }
#'  \item{\code{"arithmetic_raw_occ"}}{ - potentially drop this, as not used.}
#' }
#' @return A summary file. This .csv is saved in the output_path location and 
#' 		   contains the annual composite indicator estimate (summarized across 
#' 		   the iterations as wither the mean or median). The unique number of 
#' 		   species contributing to the indicator is shown in the "spp_num" column.
#'		   Various forms of uncertainty (estimated across the iterations) for the 
#'		   annual composite trend are presented, including the upper and lower 
#'		   95% credible intervals, SD of the mean and standard error of the mean.
#' @import ggplot2
#' @import reshape2
#' @import car
#' @export

CompositeTrend <- function(indata, output_path, trend_choice = "arithmetic_logit_occ", group_name,
                           save_iterations = "yes", TrendScale = NULL, plot_output = TRUE){
  
  load(indata)
  
  number_of_spp <- length(unique(as.character(samp_post$spp))) # How many species contribute to the indicator?
  
  # loop through iterations - later convert to array and apply across array, should be quicker #
  composite_trend <- NULL
  for (j in 1:length(unique(samp_post$iter))){
    print(j)
    temp_table <- NULL
    temp_table <- samp_post[samp_post$iter == j,]
    t_table <- temp_table[,(1:(ncol(temp_table)-2))] # convert shape of the table
    
    # arithmean on the occ scale #
    logit_temp_table <- subset(t_table, select = -c(spp, iter))
    logit_temp_table <- as.data.frame(car::logit(as.matrix(logit_temp_table)))
    
    # geomean on the occ scale #
    log_temp_table <- subset(t_table, select = -c(spp, iter))
    log_temp_table <- log(log_temp_table)
    
    # geometric mean raw occupancy #
    if(trend_choice == "geometric_raw_occ"){
      composite_trend_temp <- apply(subset(t_table, select = -c(spp, iter)), 2, geomean)
      composite_trend <- rbind(composite_trend, composite_trend_temp)
    }
    
    # arithmetic mean raw occupancy #
    if(trend_choice == "arithmetic_raw_occ"){
      composite_trend_temp <- apply(subset(t_table, select = -c(spp, iter)), 2, mean)
      composite_trend <- rbind(composite_trend, composite_trend_temp)
    }
    
    # arithmetic log odds (logit) occupancy back converted to occupancy scale with inv.logit
    if(trend_choice == "arithmetic_logit_occ"){
      composite_trend_temp <- apply(logit_temp_table, 2, mean)
      composite_trend <- rbind(composite_trend, composite_trend_temp)
    }
    
  }
  
  # if the trend is based on logit, back convert to odds (rather than occupancy) following Steve Freeman's comments #
  if(trend_choice == "arithmetic_logit_occ"){
    composite_trend <- exp(composite_trend)
  }
  
  # Are we scaling the indicator? - Scale to 100 for UK biodiversity indicators
  if(!is.null(TrendScale)){
    multiplier <- TrendScale/mean(composite_trend[,1])  # identify multiplier 
    composite_trend <- composite_trend * multiplier # scale logit arithmetic mean so mean value in year 1 = the input value for TrendScale #
  }
  
  if(save_iterations == "yes"){
    write.csv(composite_trend, file = paste(output_path, group_name, "_", trend_choice, "_composite_trend_iterations.csv", sep = "") , row.names = FALSE)
  }
    
  # save the summarised iterations #
  composite_trend_summary <- data.frame(
    year = as.numeric(gsub("X", "", colnames(composite_trend))),
    mean_occupancy = apply(composite_trend, 2, mean),
    median_occupancy = apply(composite_trend, 2, median),
    lower_5_perc_CI_occ = apply(composite_trend, 2, quan_0.05),
    upper_95_perc_CI_occ = apply(composite_trend, 2, quan_0.95),
    sd_occupancy = apply(composite_trend, 2, sd),
    sem_occupancy = apply(composite_trend, 2, sem)
  )
  
  # add species number column #
  composite_trend_summary$spp_num <- number_of_spp
  write.csv(composite_trend_summary, file = paste(output_path, group_name, "_", trend_choice, "_composite_trend_summary.csv", sep = ""), row.names = FALSE)
  
  if(plot_output == TRUE){
    
    ggplot(composite_trend_summary, aes_string(x = "year", y = "mean_occupancy")) + 
      theme_bw() +
      geom_ribbon(data = composite_trend_summary, aes_string(group = 1, ymin = "lower_5_perc_CI_occ", ymax = "upper_95_perc_CI_occ"), alpha = 0.2) +
      geom_line(size = 1, col = "black") +
      geom_point(size = 2) +
      geom_hline(yintercept = 100, linetype = "dashed") +
      ylab("Index") +
      xlab("Year") +
      scale_y_continuous(limits = c(0, max(composite_trend_summary$upper_95_perc_CI_occ)))
    
    ggsave(paste(group_name, "_", trend_choice, "_composite_trend.png", sep = ""), plot = last_plot(), path = output_path, width=6, height=4, units="in", dpi = 300)	
    
  }
  return(composite_trend_summary)
}

