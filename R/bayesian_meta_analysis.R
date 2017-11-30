#' Bayesian Meta-analysis
#' 
#' @description Use a bayesian meta-analysis to create an indicator from species index values and standard error.
#' 
#' @param data a data.frame with 4 columns in this order: species, year, index, se (standard error) 
#' @param plot Logical, should a trace plot be plotted?
#' @param model The type of model to be used. See details.
#' @param parallel if \code{TRUE} the model chains will be run in parallel using one fewer cores than
#' are availbale on the computer.
#' @param incl.model if \code{TRUE} the model is added as an attribute of the object returned
#' @param n.iter The number of iterations of the model to run. Defaults to 10,000 to avoid long run times
#' though much longer runs are usually required to reach convergence
#' @param m.scale The measurement scale of the data. The scale of the data is assumed to be logarithmic.
#' Here you specify which log scale the data is on ('loge', 'log10', or 'logit'). Defaults to 'loge'.
#' @param num.knots If using either of the smooth models this specifies the number of knots.
#' @param rescaleYr1 Logical, should all iterations be scaled so that the first year is equal? If TRUE
#' year one will have 0 error.
#' @details There are a number of model to choose from:
#' \itemize{
#'  \item{\code{"random_walk"}}{ - Also known as BMA3, strongly recommended.}
#'  \item{\code{"uniform"}}{ - Also known as BMA2.}
#'  \item{\code{"uniform_noeta"}}{ - Also known as BMA1.}
#'  \item{\code{"FNgr"}}{ - Indicator defined by Growth rates, as in Freeman & Newson.}
#'  \item{\code{"smooth_stoch"}}{ - Indicator defined by Growth rates, with Ruppert smoother (stochastic version).}
#'  \item{\code{"smooth_det"}}{ - Indicator defined by Growth rates, with Ruppert smoother (deterministic version).}
#' }
#' @return Returns a dataframe with 4 columns: Year, Index, lower2.5, upper97.5. The last two columns are the credible intervals
#' @import reshape2
#' @import jagsUI
#' @importFrom boot inv.logit
#' @importFrom coda mcmc.list as.mcmc
#' @export
#' @examples 
#' 
#' # Create some example data in the format required
#' data <- data.frame(species = rep(letters, each = 50),
#'                    year = rep(1:50, length(letters)),
#'                    index = runif(n = 50 * length(letters), min = 0, max = 1),
#'                    se = runif(n = 50 * length(letters), min = 0.01, max = .1))
#' 
#' # Run the Bayesian meta-analysis
#' bma_indicator <- bma(data)
#' 
#' # Plot the resulting indicator
#' plot_indicator(indicator = bma_indicator[,'Index'],
#'                CIs = bma_indicator[,c(3,4)])


bma <- function (data,
                 plot = TRUE,
                 model = 'random_walk',
                 parallel = FALSE,
                 incl.model = TRUE,
                 n.iter = 1e4,
                 m.scale = 'loge',
                 num.knots = 10,
                 rescaleYr1 = TRUE){
  
  if (!identical(colnames(data), c("species", "year", "index", 
                                   "se"))) {
    stop('data column names should be: "species", "year", "index", "se"')
  }
  
  # This is not my preferrred behaviour
  if(!m.scale %in% c('loge', 'log10', 'logit')) stop("m.scale must be 'loge', 'log10', or 'logit'")
  
  # pick the correct model
  switch(tolower(model),
         random_walk = {bugs_path <- bma_model_ranwalk()},
         uniform = {bugs_path <- bma_model_uniform()},
         uniform_noeta = {bugs_path <- bma_model_uniform_noeta()},
         fngr = {bugs_path <- bma_model_FNgr()},
         smooth_stoch = {bugs_path <- bma_model_smooth_stoch()},
         smooth_det = {bugs_path <- bma_model_smooth_det()},
         {stop(paste("model type not know. Must be one of 'random_walk',",
                     "'uniform', 'uniform_noeta', 'FNgr', 'smooth_stoch',",
                     "'smooth_det'"))})
  
  # 24 Feb - the index is already on the log scale (for butteflies at least)
  #index <- log(acast(data, species ~ year, value.var = "index"))
  index <- (acast(data, species ~ year, value.var = "index"))
  
  se <- acast(data, species ~ year, value.var = "se")
  
  # Setup BUGS data
  bugs_data <- list(nsp = nrow(index),
                    nyears = ncol(index),
                    estimate = index, 
                    sigma.obs = se, 
                    max_se = ifelse(test = all(is.na(se)),
                                    yes = 10,
                                    no = max(se, na.rm = TRUE))) 
  
  if(model %in% c('smooth_stoch', 'smooth_det')){
    ZX <- makeZX(num.knots = num.knots,
                 covariate = seq(min(data$year),
                                 max(data$year)))
    bugs_data[['Z']] <- ZX[['Z']]
    bugs_data[['X']] <- ZX[['X']]
    bugs_data[['num.knots']] <- num.knots
  }
  
  # Setup parameters to monitor
  params = c("tau.spi", "logI", "spindex", "sigma.obs")
  if(model %in% c('smooth_stoch', 'smooth_det')) params <- c(params, "logI.raw")
  if(model %in% c('random_walk', 'uniform', 'uniform_noeta')) params <- c(params, "tau.eta")
  if(model %in% c('random_walk')) params <- c(params, "tau.I")
  if(model %in% c('smooth_stoch', 'smooth_det', 'FNgr')) params <- c(params, "spgrowth")
  
  model <- jagsUI::jags(data = bugs_data,
                        inits = NULL,
                        param = params,
                        parallel = parallel,
                        n.cores = parallel::detectCores()-1,
                        model.file = bugs_path,
                        n.chains = 3,
                        n.thin = 2,
                        n.iter = n.iter,
                        n.burnin = floor(n.iter/2))
  
  if (plot) {
    array_sim <- model$samples
    comb.samples <- mcmc.list(lapply(1:3, FUN = function(x, 
                                                         array_sim) {
      year_ests <- colnames(array_sim[[x]])[grepl("^logI\\[",  # changed from I
                                                  colnames(array_sim[[x]]))]
      ar_temp <- array_sim[[x]][,c(head(year_ests, 1), tail(year_ests,1))]
      colnames(ar_temp) <- c("First year", "Last year")
      as.mcmc(ar_temp)
    }, array_sim = array_sim))
    plot(comb.samples)
  }
  
  # rescale year 1/ or not
  if(!rescaleYr1){ 
    pd <- data.frame(mean = unlist(model$mean),
                     q2.5 = unlist(model$q2.5),
                     q97.5 = unlist(model$q97.5))
  } else {
    logI_rescaled <- t(apply(model$sims.list$logI, 1, function(x) x - x[1]))
    pd <- data.frame(mean = apply(logI_rescaled, 2, mean),
                     q2.5 = apply(logI_rescaled, 2, quantile, probs = 0.025),
                     q97.5 = apply(logI_rescaled, 2, quantile, probs = 0.975),
                     row.names = paste0('logI', 1:ncol(logI_rescaled)))
  }
  # convert the logI back to the measurement scale  
  unb2b <- function(x, m.scale){
    switch(m.scale, 
           loge = x <- exp(x),
           log10 = x <- 10^x,
           logit = {x <- boot::inv.logit(as.matrix(x))},
           warning(paste(m.scale, 'unknown, no back-transformation applied')))
    return(x)
  }
  pd <- unb2b(pd[grepl("^logI[[:digit:]]+",dimnames(pd)[[1]]),], m.scale)
  
  rescale_bayesian_indicator <- function(x, centering = "firstyr") {
    if (centering == "firstyr") {
      x <- x/x[1, 1]
    }
    100 * x
  }
  pd <- rescale_bayesian_indicator(pd, centering = "firstyr")
  pd <- dcast(melt(as.matrix(pd)), Var1 ~ 
                Var2)
  names(pd) <- c("Year", "Index", "lower2.5", "upper97.5")
  pd$Year <- as.numeric(pd$Year)
  
  if(incl.model) attr(pd, 'model') <- model
  
  return(pd)
}