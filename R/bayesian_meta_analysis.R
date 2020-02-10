#' Bayesian Meta-analysis
#' 
#' @description Use a Bayesian meta-analysis to create an indicator from species index values, optionally incorporating standard error.
#' @param data a data.frame with 3-4 columns in this order: `species`, `year`, `index`, `se` (standard error). The `se` column is optional 
#' NB: Index values are assumed to be on the unbounded (logarithmic scale)
#' @param plot Logical, should a trace plot be plotted to diagnose the model output?
#' @param model The type of model to be used. See details.
#' @param parallel if \code{TRUE} the model chains will be run in parallel using one fewer cores than
#' are available on the computer. NOTE: this will typically not work for parallel use on cluster PCs.
#' @param incl.model if \code{TRUE} the model is added as an attribute of the object returned
#' @param n.iter The number of iterations of the model to run. Defaults to 10,000 to avoid long run times
#' though much longer runs are usually required to reach convergence
#' @param m.scale The measurement scale of the data. The scale of the data is assumed to be logarithmic.
#' Here you specify which log scale the data is on ('loge', 'log10', or 'logit'). Defaults to 'loge'.
#' @param num.knots If using either of the smooth models this specifies the number of knots.
#' @param rescaleYr Logical, should all iterations be scaled so that the first year is equal? If TRUE
#' year one will have 0 error.
#' @param seFromData Logical. Should the standard errors be read in from data (`TRUE`) or estimated (`FALSE`)?  Defaults to `FALSE` 
#' @param Y1perfect Logical. Should the first year of a species' index be assumed known without error (`TRUE`)? Defaults to `TRUE` 
#' @param incl.2deriv Logical. Option to include estimation of second derivatives on the indicator (`TRUE`)? Defaults to `FALSE` 
#' @param n.thin Thinning rate for the Markov chains. Defaults to 5.
#' @param save.sppars Logical. Should the species-specific parameters be monitored? Defaults to TRUE 
#' @param q defines the quantiles of the posterior distribution to report. Defaults to c(0.025, 0.975), i.e. the 95th percentile credible intervals
#' @details There are a number of model to choose from:
#' \itemize{
#'  \item{\code{"smooth"}}{ The default option. Indicator defined by Growth rates, with Ruppert smoother, allowing for species to join late. Error on the first year of each species' time-series is assumed to be zero. The indicator is the expected value of the geometric mean across species (with missing species imputed). 
#'  Includes three options: `seFromData` `Y1perfect` and `incl.2deriv`. See bayesian_meta_analysis for mode details.}
#'  \item{\code{"smooth_JABES"}}{ Equivalent to smooth with `seFromData = TRUE` and `Y1perfect = TRUE`. This is the version implemented in the JABES paper. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  \item{\code{"smooth_det2"}}{ Equivalent to smooth with `seFromData = TRUE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  \item{\code{"smooth_det_sigtheta"}}{ Equivalent to smooth with `seFromData = FALSE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  \item{\code{"smooth_det"}}{ Specific variant of smooth_det2 - under review. Likely to be deprecated}
#'  }
#' @return Returns a dataframe with 4 columns: Year, Index, lower2.5, upper97.5. The last two columns are the credible intervals
#' @import reshape2
#' @importFrom boot inv.logit
#' @importFrom coda mcmc.list as.mcmc
#' @references Freeman, S.N., Isaac, N.J.B., Besbeas, P.T., Dennis, E.B. & Morgan, B.J.T. (2019) 
#'             A generic method for estimating and smoothing multispecies biodiversity indices, robust to intermittent data. 
#'             \emph{JABES}, in revision.
#' @export
#' @examples 
#' Create some example data in the format required
#' data <- data.frame(species = rep(letters, each = 50),
#'                    year = rep(1:50, length(letters)),
#'                    index = rnorm(n = 50 * length(letters), mean = 0, sd = 1),
#'                    se = runif(n = 50 * length(letters), min = 0.01, max = .1))
#' 
#' Run the Bayesian meta-analysis
#' bma_indicator <- bma(data, model="smooth", m.scale="logit")
#' 
#' Plot the resulting indicator
#' plot_indicator(indicator = bma_indicator[,'Index'],
#'                CIs = bma_indicator[,c(3,4)])


bma <- function (data,
                 plot = TRUE,
                 model = 'smooth',
                 parallel = FALSE,
                 incl.model = TRUE,
                 n.iter = 1e4,
                 m.scale = 'loge',
                 num.knots = 12,
                 rescaleYr = 1,
                 seFromData = FALSE,
                 Y1perfect = TRUE,
                 incl.2deriv = FALSE,
                 save.sppars = TRUE,
                 n.thin = 5,
                 q = c(0.025, 0.975)){
  
  # Check if jagsUI is installed
  if (!requireNamespace("jagsUI", quietly = TRUE)) {
    stop("Package 'jagsUI' is needed for the 'bma' function to work. Please insatll this from CRAN. You will also be required to install JAGS, which you can download from https://sourceforge.net/projects/mcmc-jags/files/JAGS/",
         call. = FALSE)
  }
  
  if (!identical(colnames(data)[1:3], c("species", "year", "index"))) {
    stop('data column names should be: "species", "year", "index"')
  }
  
  if(colnames(data)[4] != "se" | ncol(data) < 4) {# add a set of NAs
    data$se <- NA
    if(seFromData) {
      stop("Error: Standard errors have not been provided")
    }
  }
  
  switch(tolower(model),
         smooth = {}, # the default
         smooth_jabes = {# this is the version SF ran for the paper
           model = "smooth"
           seFromData = FALSE
           Y1perfect = TRUE},
         smooth_det2 = {# this is the version NI tested for ISEC 
           model = "smooth"
           seFromData = TRUE
           Y1perfect = FALSE},
         smooth_det_sigtheta = {# this is the version NI ran for Scottish indicators
           model = "smooth"
           seFromData=FALSE
           Y1perfect = FALSE},
         smooth_det = {
           seFromData = TRUE
           Y1perfect = FALSE
         },
         random_walk = stop("Random walk model has been deprecated"),
         uniform = stop("Uniform model has been deprecated"),
         uniform_noeta = stop("Uniform model has been deprecated"),
         fngr = stop("This model option has been deprecated"),
         smooth_stoch = stop("This model option has been deprecated"),
         fngr2 = stop("This model option has been deprecated"),
         smooth_stoch2 = stop("This model option has been deprecated"),
         stop(paste("Model type not known. Check the help file for details"))
         )
  
  # do a quick check for whether the index values have been transformed
  if(min(data$index, na.rm = T) >= 0)
    print("Warning: No negative index values detected. Are you sure you transformed the data?")
  
  # check whether the data contain any infinite values
  if(any(is.infinite(data$index)))
    stop('Dataset contains Infinite values. Fix this before proceeding')
  
  # This is not my preferrred behaviour
  if(!m.scale %in% c('loge', 'log10', 'logit')) stop("m.scale must be 'loge', 'log10', or 'logit'")
  
  # pick the correct model
  model_code <- get_bmaBUGScode(option = model, 
                                incl.2deriv=incl.2deriv,
                                seFromData = seFromData, 
                                Y1perfect = Y1perfect)

  # save it to a temp file 
  bugs_path <- tempfile()
  writeLines(text = model_code, con = bugs_path)

  # include an option here to standardise the data to some value in year 1 
  
  # we assume that the index values are already on the unbounded (log) scale (we checked for negative values above)  
  index <- (acast(data, species ~ year, value.var = "index"))
  
  # Setup BUGS data
  bugs_data <- list(nsp = nrow(index),
                    nyears = ncol(index),
                    estimate = index
                    ) 
  
  if(seFromData) {
    se <- acast(data, species ~ year, value.var = "se")
    bugs_data$sigma.obs <- se
    bugs_data$max_se = ifelse(test = all(is.na(se)),
                    yes = 10,
                    no = max(se, na.rm = TRUE))
  }
  
  if(model %in% c('smooth', 'smooth_stoch', 'smooth_det', 'smooth_stoch2')){
    ZX <- makeZX(num.knots = num.knots,
                 covariate = seq(min(data$year),
                                 max(data$year)))
    bugs_data[['Z']] <- ZX[['Z']]
    bugs_data[['X']] <- ZX[['X']]
    bugs_data[['num.knots']] <- num.knots
  }

  
  if(model %in% c('smooth', 'smooth_stoch2', 'FNgr2')){
    # using row.names should ensure the same order in the bugs data
    FY <- sapply(row.names(index), FUN = function(x){
      min(data$year[!is.na(data$index) & data$species == x])
    })
    bugs_data[['FY']] <- FY - min(FY) + 1 # set lowest value to 1
    
  }
  
  # Setup parameters to monitor
  params = c("tau.spi", "logI", "sigma.obs")
  if(model %in% c('smooth', 'smooth_stoch', 'smooth_det')) params <- c(params, "logI.raw")
  if(model %in% c('random_walk', 'uniform', 'uniform_noeta')) params <- c(params, "tau.eta")
  if(model %in% c('random_walk')) params <- c(params, "tau.I")
  if(model %in% c('smooth', 'smooth_stoch', 'smooth_det', 'FNgr',
                  'smooth_stoch2', 'FNgr2')) params <- c(params, "logLambda", "spgrowth", "logI2")
  if(model %in% c('smooth', 'smooth_stoch', 'smooth_det', 'FNgr')) params <- c(params, "tau.sg")
  if(model %in% c('smooth', 'smooth_stoch', 'smooth_det','smooth_stoch2')) params <- c(params, "beta", "taub")
  if(!seFromData) params <- c(params, "theta")
  if(save.sppars) {
    params <- c(params, "spindex")
  } else {
    params <- params[!params %in% c("spgrowth", "sigma.obs")]
  }
  if(incl.2deriv) params <- c(params, "t2dash")
  
  model.out <- jagsUI::jags(data = bugs_data,
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
    array_sim <- model.out$samples
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
    pd <- data.frame(mean = unlist(model.out$mean),
                     q2.5 = unlist(model.out$q2.5),
                     q97.5 = unlist(model.out$q97.5))
  } else {
    logI_rescaled <- t(apply(model.out$sims.list$logI, 1, function(x) x - x[rescaleYr]))
    pd <- data.frame(mean = apply(logI_rescaled, 2, mean),
                     lowerCI = apply(logI_rescaled, 2, quantile, probs = q[1]),
                     upperCI = apply(logI_rescaled, 2, quantile, probs = q[2]),
                     row.names = paste0('logI', 1:ncol(logI_rescaled)))
    if(q[1] == 0.025 & q[2] == 0.975) names(pd)[2:3] <- c("q2.5", "q97.5")
  }
  # convert the logI back to the measurement scale  
  unb2b <- function(x, m.scale){
    switch(m.scale, 
           loge = x <- exp(x),
           log10 = x <- 10^x,
           logit = x <- exp(x), # Counter-intuitively, since we want geometric mean odds
           #logit = {x <- boot::inv.logit(as.matrix(x))}, # this would give occupancy, which is hard to interpret
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
  pd$Year <- as.numeric(pd$Year) + min(data$year) - 1 
  
  if(incl.model) attr(pd, 'model') <- model.out
  
  return(pd)
}