#' Bayesian Meta-analysis
#' 
#' @description Use a Bayesian meta-analysis to create an indicator from species index values, optionally incorporating standard error.
#' 
#' @param data a data.frame with 4 columns in this order: species, year, index, se (standard error) 
#' @param plot Logical, should a trace plot be plotted?
#' @param model The type of model to be used. See details.
#' @param parallel if \code{TRUE} the model chains will be run in parallel using one fewer cores than
#' are available on the computer. NOTE: this will typically not work for parallel use on cluster PCs.
#' @param incl.model if \code{TRUE} the model is added as an attribute of the object returned
#' @param n.iter The number of iterations of the model to run. Defaults to 10,000 to avoid long run times
#' though much longer runs are usually required to reach convergence
#' @param m.scale The measurement scale of the data. The scale of the data is assumed to be logarithmic.
#' Here you specify which log scale the data is on ('loge', 'log10', or 'logit'). Defaults to 'loge'.
#' @param num.knots If using either of the smooth models this specifies the number of knots.
#' @param rescaleYr1 Logical, should all iterations be scaled so that the first year is equal? If TRUE
#' year one will have 0 error.
#' @param n.thin Thinning rate for the Markov chains. Defaults to 5.
#' @param save.sppars Logical. Should the species-specific parameters be monitored? Defaults to TRUE 
#' @details There are a number of model to choose from:
#' \itemize{
#'  \item{\code{"random_walk"}}{ - Also known as BMA3.}
#'  \item{\code{"uniform"}}{ - Also known as BMA2.}
#'  \item{\code{"uniform_noeta"}}{ - Also known as BMA1.}
#'  \item{\code{"FNgr"}}{ - Indicator defined by Growth rates, as in Freeman & Newson.}
#'  \item{\code{"smooth_stoch"}}{ - Indicator defined by Growth rates, with Ruppert smoother (stochastic version).}
#'  \item{\code{"smooth_det"}}{ - Indicator defined by Growth rates, with Ruppert smoother (deterministic version).}
#'  \item{\code{"FNgr2"}}{ - Variant where species can join the series late and error on the first year is 0 (check with Nick and Steve).}
#'  \item{\code{"smooth_stoch2"}}{ - Variant where species can join the series late and error on the first year is 0 (check with Nick and Steve).}
#'  \item{\code{"smooth_det2"}}{ - Variant where species can join the series late and error on the first year is 0 (check with Nick and Steve).}
#'  \item{\code{"smooth_det_sigtheta"}}{ - Variant of det2 in which standard errors are assumed constant (check with Nick and Steve).}
#' }
#' @return Returns a dataframe with 4 columns: Year, Index, lower2.5, upper97.5. The last two columns are the credible intervals
#' @import reshape2
#' @import jagsUI
#' @importFrom boot inv.logit
#' @importFrom coda mcmc.list as.mcmc
#' @references Freeman, S.N., Isaac, N.J.B., Besbeas, P.T., Dennis, E.B. & Morgan, B.J.T. (2019) 
#'             A generic method for estimating and smoothing multispecies biodiversity indices, robust to intermittent data. 
#'             \emph{JABES}, in revision.
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
                 num.knots = 12,
                 rescaleYr = 1,
                 n.thin = 5,
                 #save.spindex = TRUE){
                 save.sppars = TRUE){
  
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
         fngr2 = {bugs_path <- bma_model_FNgr2()},
         smooth_stoch2 = {bugs_path <- bma_model_smooth_stoch2()},
         smooth_det2 = {bugs_path <- bma_model_smooth_det2()},
         smooth_det_sigtheta = {bugs_path <- bma_model_smooth_det_sigtheta()},
         {stop(paste("model type not know. Must be one of 'random_walk',",
                     "'uniform', 'uniform_noeta', 'FNgr', 'smooth_stoch',",
                     "'smooth_det', 'smooth_stoch2', 'smooth_det2', 'FNgr2', 'smooth_det_sigtheta'"))})

  
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
  
  # if(model %in% c('random_walk', 'uniform', 'uniform_noeta', 'FNgr', 'smooth_stoch', 'smooth_det')){
  #   bugs_data[['max_se']] <- ifelse(test = all(is.na(se)),
  #                                   yes = 10,
  #                                   no = max(se, na.rm = TRUE))
  # }
  
  if(model %in% c('smooth_stoch', 'smooth_det', 'smooth_det_sigtheta',
                  'smooth_stoch2', 'smooth_det2')){
    ZX <- makeZX(num.knots = num.knots,
                 covariate = seq(min(data$year),
                                 max(data$year)))
    bugs_data[['Z']] <- ZX[['Z']]
    bugs_data[['X']] <- ZX[['X']]
    bugs_data[['num.knots']] <- num.knots
  }

  
  if(model %in% c('smooth_stoch2', 'smooth_det2', 'smooth_det_sigtheta', 'FNgr2')){
    # using row.names should ensure the same order in the bugs data
    FY <- sapply(row.names(index), FUN = function(x){
      min(data$year[!is.na(data$index) & data$species == x])
    })

    bugs_data[['FY']] <- FY
  }
  
  # Setup parameters to monitor
  params = c("tau.spi", "logI", "sigma.obs")
  if(model %in% c('smooth_stoch', 'smooth_det', 'smooth_det_sigtheta')) params <- c(params, "logI.raw")
  if(model %in% c('random_walk', 'uniform', 'uniform_noeta')) params <- c(params, "tau.eta")
  if(model %in% c('random_walk')) params <- c(params, "tau.I")
  if(model %in% c('smooth_stoch', 'smooth_det', 'FNgr','smooth_det_sigtheta',
                  'smooth_stoch2', 'smooth_det2', 'FNgr2')) params <- c(params, "logLambda", "spgrowth", "logI2")
  if(model %in% c('smooth_stoch', 'smooth_det', 'FNgr', 'smooth_det_sigtheta')) params <- c(params, "tau.sg")
  if(model %in% c('smooth_stoch', 'smooth_det','smooth_stoch2', 'smooth_det2','smooth_det_sigtheta')) params <- c(params, "beta", "taub")
  if(save.sppars) {
    params <- c(params, "spindex")
  } else {
    params <- params[!params %in% c("spgrowth", "sigma.obs")]
  }

  model <- jagsUI::jags(data = bugs_data,
                        inits = NULL,
                        param = params,
                        parallel = parallel,
                        n.cores = parallel::detectCores()-1,
                        model.file = bugs_path,
                        store.data = TRUE,
                        n.chains = 3,
                        n.thin = n.thin,
                        n.iter = n.iter,
                        n.burnin = floor(n.iter/2))
  
  if (plot==TRUE) {
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
  if(rescaleYr == 0){ 
    pd <- data.frame(mean = unlist(model$mean),
                     q2.5 = unlist(model$q2.5),
                     q97.5 = unlist(model$q97.5))
  } else {
    logI_rescaled <- t(apply(model$sims.list$logI, 1, function(x) x - x[rescaleYr]))
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