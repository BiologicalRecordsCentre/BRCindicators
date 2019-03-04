#' Multi-Species Indicator tool
#' 
#' @description A functionalised version of code developed by Statistics Netherlands (see source). See `msi` for a 
#' wrapper that removes the need to write files to a folder first.
#' Multi-Species Indicators (MSI) are biodiversity indicators that combine the population
#' development of species into a single indicator. The MSI-tool calculates an MSI, confidence intervals
#' for the MSIs and linear and flexible (smoothed) trends. The trends are classified in terms like
#' "moderate increase", "strong decrease" or "stable". A number of additional analyses can be performed
#' like testing for changepoints, comparison of trends before and after a changepoint and the calculation
#' and testing of the total change in a time series.
#' 
#' @param wd The path for input and ouput files
#' @param inputFile The name of the input file. This must have columns species, year, index, se
#' in that order. The index value in the base year (which need not be the first year), should be set to 100
#' with se of 0.
#' @param jobname Generic name for output files
#' @param nsim Number of Monte Carlo simulations
#' @param SEbaseyear Desired year to set MSI to 100 and SE to 0; usually the first year of the time series
#' @param plotbaseyear desired year to set to 100 in plots
#' @param index_smooth A character, either 'INDEX' or 'SMOOTH'. "INDEX" will cause the index
#' in plotbaseyear set to 100; "SMOOTH" will set the smoothed trend value in the plotbaseyear to 100.
#' @param span Span is the proportion of points to use in the weighted estimation of the smoothed line.
#' Values close to 1 are more smoothed than values close to 0
#' @param lastyears last X years of time series for which separate trend (short-term trends) should be calculated
#' @param maxCV maximum allowed mean Coefficient of Variation (CV) of species indices (0.5 = 50%).
#'Species with higher mean CV are excluded.
#' @param changepoint compare trends before and after this year
#' @param truncfac truncation factor (=max year-to-year index ratio). Default for Living Planet Index = 10.
#' @param TRUNC set all indices below TRUNC to this value and their SE to 0. TRUNC = 0 means no truncation.
#' @param plot logical, should plots be created
#' @source https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim--/msi-tool
#' @export

msi_tool <- function(wd = getwd(),
                     inputFile,
                     jobname = 'MSI_job',
                     nsim = 1000,
                     SEbaseyear = NULL,
                     plotbaseyear = NULL,
                     index_smooth = "SMOOTH",
                     span = 0.75,
                     lastyears = 10,
                     maxCV = 3,
                     changepoint = NULL,
                     truncfac = 10,
                     TRUNC = 1,
                     plot = TRUE){
  
  org_wd <- getwd()
  on.exit({setwd(org_wd)})
  setwd(wd)
  
  rdata <- read.csv(inputFile, stringsAsFactors = FALSE)
  
  # catch changepoint error
  if(is.null(changepoint)) changepoint <- floor(max(rdata[,"year"]) - ((max(rdata[,"year"]) - min(rdata[,"year"]))/2))
  if(!changepoint %in% unique(rdata$year)) stop('changepoint year is not in range of years in data')
  
  # check SEbaseyear
  if(is.null(SEbaseyear)) SEbaseyear <- min(rdata[,"year"])
  if(!SEbaseyear %in% unique(rdata$year)) stop('SEbaseyear year is not in range of years in data')
  
  # Check plotbaseyear
  if(is.null(plotbaseyear)) plotbaseyear <- min(rdata[,"year"])
  if(!plotbaseyear %in% unique(rdata$year)) stop('plotbaseyear year is not in range of years in data')
  
  # check span
  if(span > 1 | span < 0) stop('span must be between 0 and 1')
  
  # set names of output files
  jobnameRESULTS <- paste (jobname, "_RESULTS.csv", sep = "")
  jobnameTRENDS <- paste (jobname, "_TRENDS.csv", sep = "")
  jobnameSIMTRENDS <- paste (jobname, "_SIMTRENDS.csv", sep = "")
  jobnameGRAPH <- paste (jobname, "_GRAPH.jpg", sep = "")
  
  # Define and calculate model parameters
  species <- rdata[,"species"]
  year <- rdata[,"year"]
  index <- rdata[,"index"]
  se <- rdata[,"se"]
  uspecies <- sort(unique(species))
  nspecies <- length(uspecies)
  uyear <- sort(unique(year))
  nyear <- length(unique(year))
  meanindex <- tapply(index, species, mean, na.rm=TRUE)/100
  mnindex1 <- as.data.frame(rep(meanindex,each=nyear))
  minyear <- min(year)
  maxyear <- max(year)
  plotbaseyear <- plotbaseyear-minyear+1
  baseyear <- max(1,SEbaseyear-minyear+1)
  lastyears <- min(lastyears,nyear)
  # INP <- data.frame(species, year, index, se)
  SPEC <- as.matrix(species)
  species <- sort(rep(uspecies,each=nyear))
  year <- rep(uyear,nspecies)
  # INP1 <- data.frame(species, year)
  # INP2 <- merge(INP, INP1, by=c("species","year"), sort=TRUE, all=TRUE)
  INP2 <- data.frame(species = rdata[,"species"],
                     year = rdata[,"year"],
                     index = rdata[,"index"],
                     se = rdata[,"se"])[order(rdata[,"species"],
                                         rdata[,"year"]), ]
  
  # Calculate and plot mean CV for indices per species
  CVtemp <- INP2[INP2$index >= 10, ] # select records with index >= 10 for CV calculation
  CVtemp <- CVtemp[!is.na(CVtemp$index), ] # reject missing values
  CV1 <- CVtemp$se/CVtemp$index
  CV1[CV1== 0] <- NA
  CV1[CV1== Inf] <- NA
  species2 <- CVtemp$species
  mnCV <- tapply(CV1, species2, mean, na.rm=TRUE)
  
  if(plot){
    plot.new()
    plot(1:length(uspecies), mnCV, type="p", pch=19, col="black",
         main=jobname, xlab="species code", ylab="meanCV")
  }
  
  write.csv(data.frame(species = uspecies, mean_CV = mnCV),
            file = 'species_CV_values.csv', row.names = FALSE)
  
  CV <- as.data.frame(rep(mnCV, each=nyear))
  
  # replace small indices by TRUNC and their SE by zero
  INP3 <- data.frame(cbind(species, year, CV, mnindex1))
  colnames(INP3) <- c("species","year", "CV", "mnindex1")
  INP4 <- merge(INP2, INP3, by=c("species","year"), sort=FALSE, all=TRUE)
  INP5 <- subset(INP4, CV < maxCV, select = c(species, year, index, se, mnindex1))
  INP5$index[INP5$index < TRUNC & !is.na(INP5$index)] <- TRUNC 
  INP5$se[INP5$index == TRUNC & !is.na(INP5$index)] <- 0
  
  # reset parameters
  index <- as.vector(INP5["index"])
  se <- as.vector(INP5["se"])
  nobs <- NROW(INP5)
  uspecies <- sort(unique(INP5$species))
  nspecies <- length(uspecies)
  year <- rep(uyear, nspecies)
  
  # Transform indices and standard deviations to log scale (Delta method)
  LNindex <- as.vector(log(index))
  LNse <- as.vector(se/index)
  
  # Monte Carlo simulations of species indices
  MC <- matrix(NA, nobs, nsim)
  for (s in 1:nsim) { 
    for (o in 1:nobs) {
      MC[o,s] <- rnorm(1, LNindex[o,1], LNse[o,1])
    } 
  }
  MC[MC < log(TRUNC)] <- log(TRUNC)
  
  # impute missing values using chain method
  CHAIN <- matrix(NA, nobs, nsim)
  for (s in 1:nsim ){
    for (o in 1:nobs-1) {
      CHAIN[o,s] <- MC[o+1,s]-MC[o,s]
    }
    for (sp in 1:nspecies) {
      CHAIN[sp*nyear,s] <- NA
    }
  }
  
  CHAIN[CHAIN > log(truncfac)] <- log(truncfac)
  CHAIN[CHAIN < log(1/truncfac)] <- log(1/truncfac)
  
  mnCHAIN <- matrix(NA, nyear, nsim)
  for (s in 1:nsim) {
    mnCHAIN[,s] <- tapply(CHAIN[,s], year, mean, na.rm=TRUE)
  }
  
  simMSI <- matrix(NA, nyear, nsim)
  bmin <- baseyear-1
  bplus <- baseyear+1
  for (s in 1:nsim){
    simMSI[baseyear,s] <- log(100)
    for (y in bmin:min(bmin,1)){
      simMSI[y,s] <- simMSI[y+1,s] - mnCHAIN[y,s]
    }
    for (y in min(bplus,nyear):nyear){
      simMSI[y,s] <- simMSI[y-1,s] + mnCHAIN[y-1,s]
    }
  }
  
  # calculate MSI and SE
  meanMSI <- array(NA, dim=c(1, nyear))
  stdevMSI <- array(NA, dim=c(1, nyear))   
  for (y in 1:nyear) {
    meanMSI[y] <- mean(simMSI[y,])
    stdevMSI[y] <- sd(simMSI[y,])
  }
  
  # Monte Carlo simulations of MSI (for trend calculation)
  simMSI <- matrix(NA, nyear, nsim)
  for (s in 1:nsim) { 
    for (y in 1:nyear) {
      simMSI[y,s] <- rnorm(1, meanMSI[y], stdevMSI[y])
    } 
  }
  
  # Back-transformation of MSI to index scale (Delta method)
  meanMSI <- array(NA, dim=c(1, nyear))
  stdevMSI <- array(NA, dim=c(1, nyear))   
  for (y in 1:nyear) {
    meanMSI[y] <- round(exp(mean(simMSI[y,])), digits=2)
    stdevMSI[y] <- round(sd(simMSI[y,])*meanMSI[y], digits=2)
  }
  
  # Confidence interval of MSI on index-scale
  CI <- array(NA, dim=c(nyear, 2))
  for (y in 1:nyear){
    CI[y,1] <- exp(mean(simMSI[y,])-1.96*sd(simMSI[y,]))
    CI[y,2] <- exp(mean(simMSI[y,])+1.96*sd(simMSI[y,]))
  }
  
  # Overall linear trend per simulation
  lrMSI <- array(NA, dim=c(2, nsim))
  for (s in 1:nsim) {
    summ <- summary(lm(simMSI[,s] ~ uyear, na.action=na.exclude))
    lrMSI[1,s] <- round(summ$coefficients[2], digits=4)	# coeff[2] = slope
    lrMSI[2,s] <- round(summ$coefficients[2,2], digits=4) 	# coeff[2,2] = se of slope 
  }
  
  SLOPE <- mean(lrMSI[1,])
  sdSLOPE <- sd(as.vector(lrMSI[1,]))
  
  # Back transform trends to index scale:
  SLOPE_mult <- round(exp(SLOPE), digits=4)
  sdSLOPE_mult <- round(sdSLOPE*SLOPE_mult, digits=4)
  
  SLOPE
  sdSLOPE
  SLOPE_mult
  sdSLOPE_mult
  
  # Overall trend classification (see Soldaat et al. 2007 (J. Ornithol. DOI 10.1007/s10336-007-0176-7))
  TrendClass <- "X"
  if (SLOPE_mult - 1.96*sdSLOPE_mult > 1.05) {TrendClass <- "strong increase" } else
    if (SLOPE_mult + 1.96*sdSLOPE_mult < 0.95) { TrendClass <- "steep decline" } else
      if (SLOPE_mult - 1.96*sdSLOPE_mult > 1.00) { TrendClass <- "moderate increase"} else
        if (SLOPE_mult + 1.96*sdSLOPE_mult < 1.00) { TrendClass <- "moderate decline" } else
          if ((SLOPE_mult - 1.96*sdSLOPE_mult - 0.95)*(1.05 - SLOPE_mult + 1.96*sdSLOPE_mult) < 0.00) { TrendClass <- "uncertain"} else
            if ((SLOPE_mult + 1.96*sdSLOPE_mult) - (SLOPE_mult - 1.96*sdSLOPE_mult) > 0.10) { TrendClass <- "uncertain" } else
            {TrendClass <- "stable"}
  TrendClass
  
  # Short term linear trend per simulation
  lrMSI_short <- array(NA, dim=c(2, nsim))
  lyear <- c((maxyear-lastyears+1):maxyear)
  simMSI_short <- array(NA, dim=c(lastyears, nsim))
  for (s in 1:nsim) {
    for (y in 1:lastyears) {
      simMSI_short[y, s] <- simMSI[(nyear-lastyears+y),s]
    }
  }
  for (s in 1:nsim) {
    summ <- summary(lm(simMSI_short[,s] ~ lyear, na.action=na.exclude))
    lrMSI_short[1,s] <- round(summ$coefficients[2], digits=4)	# coeff[2] = slope
    lrMSI_short[2,s] <- round(summ$coefficients[2,2], digits=4) 	# coeff[2,2] = sd of slope 
  }
  
  SLOPE_short <- mean(lrMSI_short[1,])
  sdSLOPE_short <- sd(as.vector(lrMSI_short[1,]))
  
  # Back-transform short-term trends to index scale
  SLOPE_short_mult <- round(exp(SLOPE_short), digits=4)
  sdSLOPE_short_mult <- round(sdSLOPE_short*SLOPE_short_mult, digits=4)
  
  SLOPE_short
  sdSLOPE_short
  SLOPE_short_mult
  sdSLOPE_short_mult
  
  # Short term trend classification (Soldaat et al. 2007 (J. Ornithol. DOI 10.1007/s10336-007-0176-7))
  TrendClass_short <- "X"
  if (SLOPE_short_mult - 1.96*sdSLOPE_short_mult > 1.05) {TrendClass_short <- "strong increase" } else
    if (SLOPE_short_mult + 1.96*sdSLOPE_short_mult < 0.95) { TrendClass_short <- "steep decline" } else
      if (SLOPE_short_mult - 1.96*sdSLOPE_short_mult > 1.00) { TrendClass_short <- "moderate increase"} else
        if (SLOPE_short_mult + 1.96*sdSLOPE_short_mult < 1.00) { TrendClass_short <- "moderate decline" } else
          if ((SLOPE_short_mult - 1.96*sdSLOPE_short_mult - 0.95)*(1.05 - SLOPE_short_mult + 1.96*sdSLOPE_short_mult) < 0.00) { TrendClass_short <- "uncertain"} else
            if ((SLOPE_short_mult + 1.96*sdSLOPE_short_mult) - (SLOPE_short_mult - 1.96*sdSLOPE_short_mult) > 0.10) { TrendClass_short <- "uncertain" } else
            {TrendClass_short <- "stable"}
  TrendClass_short
  
  # linear trend before changepoint per simulation
  lrMSI_before <- array(NA, dim=c(2, nsim))
  byear <- c(minyear:changepoint)
  byears <- length(byear)
  simMSI_before <- array(NA, dim=c(byears, nsim))
  for (s in 1:nsim) {
    for (y in 1:byears) {
      simMSI_before[y, s] <- simMSI[y,s]
    }
  }
  for (s in 1:nsim) {
    summb <- summary(lm(simMSI_before[,s] ~ byear, na.action=na.exclude))
    lrMSI_before[1,s] <- round(summb$coefficients[2], digits=4)  # coeff[2] = slope
    lrMSI_before[2,s] <- round(summb$coefficients[2,2], digits=4) 	# coeff[2,2] = sd of slope 
  }
  SLOPE_before <- round(mean(lrMSI_before[1,]), digits=5)
  sdSLOPE_before <- round(sd(as.vector(lrMSI_before[1,])), digits=5)
  
  
  # Back-transform trends before changepoint to index scale
  SLOPE_before_mult <- round(exp(SLOPE_before), digits=4)
  sdSLOPE_before_mult <- round(sdSLOPE_before*SLOPE_before_mult, digits=4)
  
  SLOPE_before
  sdSLOPE_before
  SLOPE_before_mult
  sdSLOPE_before_mult
  
  # Trend classification before changepoint
  TrendClass_before <- "X"
  if (SLOPE_before_mult - 1.96*sdSLOPE_before_mult > 1.05) {TrendClass_before <- "strong increase" } else
    if (SLOPE_before_mult + 1.96*sdSLOPE_before_mult < 0.95) { TrendClass_before <- "steep decline" } else
      if (SLOPE_before_mult - 1.96*sdSLOPE_before_mult > 1.00) { TrendClass_before <- "moderate increase"} else
        if (SLOPE_before_mult + 1.96*sdSLOPE_before_mult < 1.00) { TrendClass_before <- "moderate decline" } else
          if ((SLOPE_before_mult - 1.96*sdSLOPE_before_mult - 0.95)*(1.05 - SLOPE_before_mult + 1.96*sdSLOPE_before_mult) < 0.00) { TrendClass_before <- "uncertain"} else
            if ((SLOPE_before_mult + 1.96*sdSLOPE_before_mult) - (SLOPE_before_mult - 1.96*sdSLOPE_before_mult) > 0.10) { TrendClass_before <- "uncertain" } else
            {TrendClass_before <- "stable"}
  TrendClass_before
  
  # linear trend after changepoint per simulation
  lrMSI_after <- array(NA, dim=c(2, nsim))
  ayear <- c(changepoint:maxyear)
  ayears <- length(ayear)
  simMSI_after <- array(NA, dim=c(ayears, nsim))
  for (s in 1:nsim) {
    for (y in 1:ayears) {
      simMSI_after[y, s] <- simMSI[(byears-1+y),s]
    }
  }
  for (s in 1:nsim) {
    summa <- summary(lm(simMSI_after[,s] ~ ayear, na.action=na.exclude))
    lrMSI_after[1,s] <- round(summa$coefficients[2], digits=4)  # coeff[2] = slope
    lrMSI_after[2,s] <- round(summa$coefficients[2,2], digits=4)   # coeff[2,2] = sd of slope 
  }
  SLOPE_after <- round(mean(lrMSI_after[1,]), digits=5)
  sdSLOPE_after <- round(sd(as.vector(lrMSI_after[1,])),digits=5)
  
  # Trend classification before changepoint
  SLOPE_after_mult <- round(exp(SLOPE_after), digits=4)
  sdSLOPE_after_mult <- round(sdSLOPE_after*SLOPE_after_mult, digits=4)
  
  SLOPE_after
  sdSLOPE_after
  SLOPE_after_mult
  sdSLOPE_after_mult
  
  # Trend classification after changepoint
  TrendClass_after <- "X"
  if (SLOPE_after_mult - 1.96*sdSLOPE_after_mult > 1.05) {TrendClass_after <- "strong increase" } else
    if (SLOPE_after_mult + 1.96*sdSLOPE_after_mult < 0.95) { TrendClass_after <- "steep decline" } else
      if (SLOPE_after_mult - 1.96*sdSLOPE_after_mult > 1.00) { TrendClass_after <- "moderate increase"} else
        if (SLOPE_after_mult + 1.96*sdSLOPE_after_mult < 1.00) { TrendClass_after <- "moderate decline" } else
          if ((SLOPE_after_mult - 1.96*sdSLOPE_after_mult - 0.95)*(1.05 - SLOPE_after_mult + 1.96*sdSLOPE_after_mult) < 0.00) { TrendClass_after <- "uncertain"} else
            if ((SLOPE_after_mult + 1.96*sdSLOPE_after_mult) - (SLOPE_after_mult - 1.96*sdSLOPE_after_mult) > 0.10) { TrendClass_after <- "uncertain" } else
            {TrendClass_after <- "stable"}
  TrendClass_after
  
  # compare linear trends before and after changepoint
  compare <- array(NA, dim=c(nsim, 1))
  for (s in 1:nsim){
    compare[s,1] <- lrMSI_after[1,s]-lrMSI_before[1,s]
  }
  meanTrendDiff <- mean(compare)
  seTrendDiff <- sd(as.vector(compare))
  significance <- NA
  if (abs(meanTrendDiff)-2.58*seTrendDiff > 0) { significance <- "p<0.01" } else
    if (abs(meanTrendDiff)-1.96*seTrendDiff > 0) { significance <- "p<0.05" } else
    {significance <- "n.s."}
  
  # Smoothing
  # span = 0.75 is default in R and usually fits well. But trying other values may give better results, depending on your data
  loessMSI <- array(NA, dim=c(nyear, nsim))
  Diff <- array(NA, dim=c(nyear, nsim))
  for (s in 1:nsim) {
    # Span is the proportion of points to use in the weighted estimation
    smooth <- predict(loess(simMSI[,s]~uyear, span=span, degree=2, na.action=na.exclude),
                      data.frame(uyear), se=TRUE)
    loessMSI[,s] <- round(smooth$fit, digits=4)
    
    for (y in 1:nyear) {
      Diff[y,s] <- loessMSI[nyear,s] - loessMSI[y,s]
    } 
  }
  
  # calculate overall change (first year set to 100) and change in last years
  pct_CHANGE <- array(NA, dim=c(nsim, 2))
  for (s in 1:nsim) {
    pct_CHANGE[s,1] <- exp(loessMSI[nyear,s])*100/exp(loessMSI[1,s])-100
    pct_CHANGE[s,2] <- exp(loessMSI[nyear,s])*100/exp(loessMSI[1,s])-exp(loessMSI[(nyear-lastyears+1),s])*100/exp(loessMSI[1,s])
  }
  pct_CHANGE_long <- round(mean(pct_CHANGE[,1]),digits=3)
  sd_pct_CHANGE_long <- round(sd(pct_CHANGE[,1]), digits=3)
  pct_CHANGE_long
  sd_pct_CHANGE_long
  pct_CHANGE_short <- round(mean(pct_CHANGE[,2]),digits=3)
  sd_pct_CHANGE_short <- round(sd(pct_CHANGE[,2]), digits=3)
  pct_CHANGE_short
  sd_pct_CHANGE_short
  
  # significance of percentage change
  significance_PCT_long <- NA
  if (abs(pct_CHANGE_long)-2.58*sd_pct_CHANGE_long > 0) { significance_PCT_long <- "p<0.01" } else
    if (abs(pct_CHANGE_long)-1.96*sd_pct_CHANGE_long > 0) { significance_PCT_long <- "p<0.05" } else
    {significance_PCT_long <- "n.s."}
  
  significance_PCT_short <- NA
  if (abs(pct_CHANGE_short)-2.58*sd_pct_CHANGE_short > 0) { significance_PCT_short <- "p<0.01" } else
    if (abs(pct_CHANGE_short)-1.96*sd_pct_CHANGE_short > 0) { significance_PCT_short <- "p<0.05" } else
    {significance_PCT_short <- "n.s."}
  
  # create output for flexible trend estimates
  smoothMSI <- array(NA, dim=c(nyear, 13))
  for (y in 1:nyear) {
    smoothMSI[y,1] <- round(mean(loessMSI[y,]), digits=4) # smooth MSI on log scale
    smoothMSI[y,2] <- round(sd(loessMSI[y,]), digits=4) # SE of smooth MSI on log scale
    smoothMSI[y,3] <- round(mean(Diff[y,]), digits=4) # Difference smooth MSI with last year on log scale
    smoothMSI[y,4] <- round(sd(Diff[y,]), digits=4) # SE of difference with ast year
    smoothMSI[y,12] <- round(smoothMSI[y,1]-1.96*smoothMSI[y,2], digits=2) # lower CI of smooth MSI on log scale
    smoothMSI[y,13] <- round(smoothMSI[y,1]+1.96*smoothMSI[y,2], digits=2) # upper CI of smooth MSI on log scale
  }
  
  # Trendclassification, based on smoothed trends (Soldaat et al. 2007 (J. Ornithol. DOI 10.1007/s10336-007-0176-7))
  for (y in 1:nyear) {
    smoothMSI[y,5] <- round(smoothMSI[nyear,1]/smoothMSI[y,1], digits=4)		# = TCR
    smoothMSI[y,6] <- round(smoothMSI[y,5]-1.96*(smoothMSI[y,4]/smoothMSI[y,1]), digits=4)	# = CI- TCR
    smoothMSI[y,7] <- round(smoothMSI[y,5]+1.96*(smoothMSI[y,4]/smoothMSI[y,1]), digits=4)	# = CI+ TCR
    smoothMSI[y,8] <- round(exp(log(smoothMSI[y,5])/(nyear-y)), digits=4)		# = YCR
    smoothMSI[y,9] <- round(exp(log(smoothMSI[y,6])/(nyear-y)), digits=4)		# = CI- YCR
    smoothMSI[y,10] <- round(exp(log(smoothMSI[y,7])/(nyear-y)), digits=4)		# = CI+ YCR
  }
  
  for (y in 1:(nyear-1)) {
    if (smoothMSI[y,9] > 1.05) {smoothMSI[y,11] <- 1} else
      if (smoothMSI[y,10] < 0.95) {smoothMSI[y,11] <- 6} else
        if (smoothMSI[y,9] > 1.00) {smoothMSI[y,11] <- 2} else
          if (smoothMSI[y,10] < 1.00) {smoothMSI[y,11] <- 5} else
            if ((smoothMSI[y,9] - 0.95)*(1.05 - smoothMSI[y,10]) < 0.00) {smoothMSI[y,11] <- 3} else
              if ((smoothMSI[y,10]) - (smoothMSI[y,9]) > 0.10) {smoothMSI[y,11] <- 3} else
              {smoothMSI[y,11] <- 4}
  } 
  
  TrendClass_flex <- matrix(NA, nrow= nyear, ncol=1)
  for (y in 1:(nyear-1)) {
    if (smoothMSI[y,9] > 1.05) {TrendClass_flex[y] <- "strong_increase" } else
      if (smoothMSI[y,10] < 0.95) { TrendClass_flex[y] <- "steep_decline" } else
        if (smoothMSI[y,9] > 1.00) { TrendClass_flex[y] <- "moderate_increase"} else
          if (smoothMSI[y,10] < 1.00) { TrendClass_flex[y] <- "moderate_decline" } else
            if ((smoothMSI[y,9] - 0.95)*(1.05 - smoothMSI[y,10]) < 0.00) { TrendClass_flex[y] <- "uncertain"} else
              if ((smoothMSI[y,10]) - (smoothMSI[y,9]) > 0.10) { TrendClass_flex[y] <- "uncertain" } else
              {TrendClass_flex[y] <- "stable"}
  } 
  
  # rescaling (for presentation of plot)
  rescale <- NA
  if (index_smooth =="INDEX") {rescale <- 100/meanMSI[plotbaseyear]} else
    if (index_smooth =="SMOOTH") {rescale <- 100/exp(smoothMSI[plotbaseyear,1])} else
    {rescale <- NA}
  simMSImean <- round(as.vector(rescale*meanMSI), digits=2)
  simMSIsd <- round(as.vector(rescale*stdevMSI), digits=2)
  uppCI_MSI <- round(rescale*CI[,2], digits=2)
  lowCI_MSI <- round(rescale*CI[,1], digits=2)
  trend_flex <- round(rescale*exp(smoothMSI[,1]), digits=2)
  lowCI_trend_flex <- round(rescale*exp(smoothMSI[,12]), digits=2)
  uppCI_trend_flex <- round(rescale* exp(smoothMSI[,13]), digits=2)
  
  # create output file for MSI + smoothed trend
  RES <- as.data.frame(cbind(uyear, simMSImean, simMSIsd, lowCI_MSI, uppCI_MSI, trend_flex, lowCI_trend_flex, uppCI_trend_flex))
  RES$trend_class <- TrendClass_flex
  names(RES) <- c("year", "MSI", "sd_MSI", "lower_CL_MSI", "upper_CL_MSI", "Trend", "lower_CL_trend", "upper_CL_trend", "trend_class")
  write.csv2(RES, file=jobnameRESULTS, row.names=FALSE, quote=FALSE)
  
  # create output file with all linear trends
  SIMTRENDS <- t(rbind(lrMSI[1,],lrMSI_short[1,]))
  write.csv2(SIMTRENDS, file=jobnameSIMTRENDS, row.names=FALSE)
  
  # create output for linear trend estimates and % change
  rownames <- rbind("overall trend","SE overall trend",
                    paste("trend last", lastyears, "years"),
                    paste("SE trend last", lastyears, "years"),
                    paste("changepoint", paste0('(', changepoint, ')')),
                    paste("trend before changepoint", paste0('(', changepoint, ')')),
                    paste("SE trend before changepoint", paste0('(', changepoint, ')')),
                    paste("trend after changepoint", paste0('(', changepoint, ')')),
                    paste("SE trend after changepoint", paste0('(', changepoint, ')')),
                    "% change", "SE % change",
                    paste("% change last", lastyears, "years"),
                    paste("SE % change last", lastyears, "years"))
  output <- cbind(SLOPE_mult, sdSLOPE_mult, SLOPE_short_mult, sdSLOPE_short_mult,
                  changepoint, SLOPE_before_mult, sdSLOPE_before_mult, SLOPE_after_mult, sdSLOPE_after_mult,
                  pct_CHANGE_long, sd_pct_CHANGE_long, pct_CHANGE_short, sd_pct_CHANGE_short)
  TRENDS <- as.data.frame(t(output), row.names=rownames)
  TRENDS$significance <- ""
  TRENDS["changepoint", "significance"] <- significance
  TRENDS[1, "significance"] <- TrendClass
  TRENDS[3, "significance"] <- TrendClass_short
  TRENDS[6, "significance"] <- TrendClass_before
  TRENDS[8, "significance"] <- TrendClass_after
  TRENDS[10, "significance"] <- significance_PCT_long
  TRENDS[12, "significance"] <- significance_PCT_short
  names(TRENDS) <- c("value", "significance")
  write.csv2(TRENDS, file=jobnameTRENDS)
  
  # create output file with plot
  if(plot){
    maxy <- max(uppCI_MSI) + 10
    jpeg(filename = jobnameGRAPH, width = 480, height = 480,units = "px", pointsize = 12, bg = "white", res = NA, restoreConsole = TRUE)
    legend1 <- paste(nspecies, "species")
    legend2 <- paste(TrendClass)
    legend3 <- paste("last", lastyears, "years:", TrendClass_short)
    plot(uyear, simMSImean, type="p", pch=19, col="black", main=jobname, xlab="", ylab="Index", xlim = c(minyear, maxyear), ylim = c(0,maxy))
    text(minyear, 40, legend1, adj=0, cex=0.8)
    text(minyear, 25, legend2, adj=0, cex=0.8)
    text(minyear, 10, legend3, adj=0, cex=0.8)
    lines(uyear, rescale*exp(smoothMSI[,1]), lty=1, col="black")
    lines(uyear, rescale*exp(smoothMSI[,12]), lty=3, col="black")
    lines(uyear, rescale*exp(smoothMSI[,13]), lty=3, col="black")
    #arrows(uyear, simMSImean, uyear, (simMSImean-simMSIsd), angle = 90, code = 3, length=0)
    #arrows(uyear, simMSImean, uyear, (simMSImean+simMSIsd), angle = 90, code = 3, length=0)
    dev.off()
    
    # plot to screen
    maxy <- max(uppCI_MSI) + 10
    plot(uyear, simMSImean, type="p", pch=19, col="black", main=jobname, xlab="", ylab="Index", xlim = c(minyear, maxyear), ylim = c(0,maxy))
    text(minyear, 40, legend1, adj=0, cex=0.8)
    text(minyear, 25, legend2, adj=0, cex=0.8)
    text(minyear, 10, legend3, adj=0, cex=0.8)
    lines(uyear, rescale*exp(smoothMSI[,1]), lty=1, col="black")
    lines(uyear, rescale*exp(smoothMSI[,12]), lty=3, col="black")
    lines(uyear, rescale*exp(smoothMSI[,13]), lty=3, col="black")
    arrows(uyear, simMSImean, uyear, (simMSImean-simMSIsd), angle = 90, code = 3, length=0)
    arrows(uyear, simMSImean, uyear, (simMSImean+simMSIsd), angle = 90, code = 3, length=0)
  }
}