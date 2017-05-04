#' Bayesian Meta-analysis
#' 
#' @description Use a bayesian meta-analysis to create an indicator from species index values and standard error.
#' 
#' @param data a data.frame with 4 columns in this order: species, year, index, se (standard error) 
#' @param plot Logical, should a trace plot be plotted?
#' @param model The type of model to be used. The default of 'random_walk' is strongly recommended.
#' Two other models 'uniform' and 'uniform_noeta' exist but are almost always inferior.
#' @param parallel if \code{TRUE} the model chains will be run in parallel using one fewer cores than
#' are availbale on the computer.
#' @param incl.model if \code{TRUE} the model is added as an attribute of the object returned
#' @param n.iter The number of iterations of the model to run. Defaults to 10,000 to avoid long run times
#' though much longer runs are usually required to reach convergence
#' @param m.scale The measurement scale of the data. The scale of the data is assumed to be logarithmic.
#' Here you specify which log scale the data is on ('loge', 'log10', or 'logit'). Defaults to 'loge'.
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
                 m.scale = 'loge'){
  
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
         {stop("model type not know. Must be one of 'random_walk', 'uniform', 'uniform_noeta'")})
  
  #index <- log(acast(na.omit(data), species ~ year, value.var = "index"))
  #se <- acast(na.omit(data), species ~ year, value.var = "se")
  
  # 24 Feb - the index is already on the log scale (for butteflies at least)
  #index <- log(acast(data, species ~ year, value.var = "index"))
  index <- (acast(data, species ~ year, value.var = "index"))
  
  se <- acast(data, species ~ year, value.var = "se")
  
  bugs_data <- list(Nsp = nrow(index), Nyr = ncol(index), estimate = index, 
                    se = se, max_se = max(se)) 
  
  params = c("tau.b0", "tau.eta", "tau.I", "logI", "alpha")
  
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
  
  pd <- data.frame(mean = unlist(model$mean),
                   q2.5 = unlist(model$q2.5),
                   q97.5 = unlist(model$q97.5))
  
  # convert the logI back to the measurement scale  
  unb2b <- function(x, m.scale){
    switch(m.scale, 
           loge = x <- exp(x),
           log10 = x <- 10^x,
           logit = {x <- boot::inv.logit(as.matrix(x))},
           warning(paste(m.scale, 'unknown, no back-transformation applied')))
    return(x)
  }
  pd <- unb2b(pd[grepl("^logI",dimnames(pd)[[1]]),], m.scale)
  
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