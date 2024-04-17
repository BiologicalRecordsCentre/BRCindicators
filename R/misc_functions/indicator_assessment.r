indicator_assessment <- function(summary_table,
                                 start_year = NULL,
                                 end_year = NULL){
  
  if(is.null(start_year)) start_year <- min(summary_table$year)
  if(is.null(end_year)) end_year <- max(summary_table$year)
  
  # Get the indicator value at the start
  start_ind <- summary_table$indicator[summary_table$year == start_year]
  
  # Get the confidence intervals for the end of the period
  end_CIs <- summary_table[summary_table$year == end_year, c('lower', 'upper')]
  
  # assess
  if(start_ind < end_CIs$upper & start_ind > end_CIs$lower){
    assessment <- 'stable'
  } else if(start_ind < end_CIs$lower){
    assessment <- 'increasing'
  } else if(start_ind > end_CIs$upper){
    assessment <- 'decreasing'
  }
  
  # Return the assessment
  assessment <- data.frame(start_index = start_ind,
                           end_lower = end_CIs$lower,
                           end_upper = end_CIs$upper,
                           assessment = assessment)
  
}