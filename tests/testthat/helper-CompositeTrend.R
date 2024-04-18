# Fake data  
n_species <- 5
n_years <- 10
n_iterations <- 3

# Create a dataframe
species <- paste0("Species_", 1:n_species)
years <- paste0("X", 1985:(1985 + n_years - 1))
iterations <- 1:n_iterations

# Expand grid to create all combinations
data_expanded <- expand.grid(spp = species, year = years, iter = iterations)

# Generate random occupancy estimates
data_expanded$occupancy <- runif(nrow(data_expanded), 0.001, 0.999)

# Pivot wider to get the desired format
fake_data <- reshape2::dcast(data_expanded, spp + iter ~ year, value.var = "occupancy")

# stemmed CompositeTrend_stemmed to test that the dimensions are correct
CompositeTrend_objects <- function(indata, output_path, trend_choice = "arithmetic_logit_occ", group_name,
                           save_iterations = "yes", TrendScale = NULL, plot_output = TRUE) {
  
  # Create an empty list to save outputs
  output = list()

  # read in a csv file
  samp_post = read.csv(indata)[, -1]

  number_of_spp <- length(unique(as.character(samp_post$spp))) # How many species contribute to the indicator?

  # loop through iterations - later convert to array and apply across array, should be quicker #
  composite_trend <- NULL
  for (i in 1:length(unique(samp_post$iter))) {
    print(paste0("Running iteration: ", unique(samp_post$iter)[i]))

    temp_table <- NULL
    temp_table <- samp_post[samp_post$iter == i,]

    t_table <- subset(temp_table, select = -c(spp, iter))

  # Check if any column is non-numeric
  non_numeric_check <- !sapply(t_table, is.numeric)
  if (any(non_numeric_check)) {
    non_numeric_cols <- names(t_table)[non_numeric_check]
    
    # Stop the operation and show the non-numeric columns found
    stop(paste0("Column(s) contains non-numeric fields: ", paste(non_numeric_cols, collapse = ", ")))
  }

    # arithmean on the occ scale
    logit_temp_table <- as.data.frame(car::logit(as.matrix(t_table)))

    # geomean on the occ scale
    log_temp_table <- t_table
    log_temp_table <- log(log_temp_table)

    # geometric mean raw occupancy
    if (trend_choice == "geometric_raw_occ") {
      composite_trend_temp <- apply(t_table, 2, geomean)
      composite_trend <- rbind(composite_trend, composite_trend_temp)
    }

    # arithmetic mean raw occupancy #
    if (trend_choice == "arithmetic_raw_occ") {
      composite_trend_temp <- apply(t_table, 2, mean)
      composite_trend <- rbind(composite_trend, composite_trend_temp)
    }

    # arithmetic log odds (logit) occupancy back converted to occupancy scale with inv.logit
    if (trend_choice == "arithmetic_logit_occ") {
      composite_trend_temp <- apply(logit_temp_table, 2, mean)
      composite_trend <- rbind(composite_trend, composite_trend_temp)
    }
  }

  # if the trend is based on logit, back convert to odds (rather than occupancy) following Steve Freeman's comments #
  if (trend_choice == "arithmetic_logit_occ") {
    composite_trend <- exp(composite_trend)
  }

  # Are we scaling the indicator? - Scale to 100 for UK biodiversity indicators
  if (!is.null(TrendScale)) {
    multiplier <- TrendScale / mean(composite_trend[, 1]) # identify multiplier
    composite_trend <- composite_trend * multiplier # scale logit arithmetic mean so mean value in year 1 = the input value for TrendScale #
  }

  if (save_iterations == "yes") {
    write.csv(composite_trend, file = paste(output_path, group_name, "_", trend_choice, "_composite_trend_iterations.csv", sep = ""), row.names = FALSE)
  }

  output$composite_trend = composite_trend

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

  output$composite_trend_summary = composite_trend_summary

  if (plot_output == TRUE) {
    plot = ggplot(composite_trend_summary, aes_string(x = "year", y = "mean_occupancy")) +
      theme_bw() +
      geom_ribbon(data = composite_trend_summary, aes_string(group = 1, ymin = "lower_5_perc_CI_occ", ymax = "upper_95_perc_CI_occ"), alpha = 0.2) +
      geom_line(size = 1, col = "black") +
      geom_point(size = 2) +
      geom_hline(yintercept = 100, linetype = "dashed") +
      ylab("Index") +
      xlab("Year") +
      scale_y_continuous(limits = c(0, max(composite_trend_summary$upper_95_perc_CI_occ)))

    ggsave(paste(group_name, "_", trend_choice, "_composite_trend.png", sep = ""), plot = plot, path = output_path, width = 6, height = 4, units = "in", dpi = 300)
  
  output$plot = plot

  }
  return(output)
}
