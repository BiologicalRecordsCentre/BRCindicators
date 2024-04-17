#' Bayesian Meta-analysis
#' 
#' @description Use a Bayesian meta-analysis to create an indicator from species index values, optionally incorporating standard error.
#' @param data a data.frame with 3-4 columns: `species`, `year`, `index`, `se` (standard error). The `se` column is optional 
#' NB: Index values are assumed to be on the unbounded (logarithmic scale)
#' @param plot Logical, should a trace plot be plotted to diagnose the model output?
#' @param model The type of model to be used. See details.
#' @param parallel if \code{TRUE} the model chains will be run in parallel using one fewer cores than
#' are available on the computer as default. NOTE: this will typically not work for parallel use on cluster PCs.
#' @param n.cores if running the code in parallel this option specifies the number of cores to use. 
#' @param incl.model if \code{TRUE} the model is added as an attribute of the object returned
#' @param n.iter The number of iterations of the model to run. Defaults to 10,000 to avoid long run times
#' though much longer runs are usually required to reach convergence
#' @param n.thin Thinning rate for the Markov chains. Defaults to 5.
#' @param m.scale The measurement scale of the data. The scale of the data is assumed to be logarithmic.
#' Here you specify which log scale the data is on ('loge', 'log10', or 'logit'). Defaults to 'loge'.
#' @param num.knots If using either of the smooth models this specifies the number of knots.
#' @param rescaleYr Integer. To which year should the indicator use as a reference value (i.e. baseline). Values greater than the number of years in the dataset will be set to the final year. Defaults to 1 (the first year)
#' @param baseline Integer. What is the value of the indicator in the baseline year (defaults to 100)
#' @param errorY1 Logical. Should the indicator be presented with (`TRUE`) or without (`FALSE`) uncertainty in the baseline year. Defaults to `FALSE`.
#' @param seFromData Logical. Should the standard errors be read in from data (`TRUE`) or estimated (`FALSE`)?  Defaults to `FALSE` 
#' @param Y1perfect Logical. Should the first year of a species' index be assumed known without error (`TRUE`)? Defaults to `TRUE` 
#' @param rescale_indices Integer. A value for standardising each species time-series to start at a common value (e.g. 0). Defaults to NULL (i.e. no standardisation)
#' @param incl.2deriv Logical. Option to include estimation of second derivatives on the indicator (`TRUE`)? Defaults to `FALSE` 
#' @param save.sppars Logical. Should the species-specific parameters be monitored? Defaults to TRUE 
#' @param CI defines the credible intervals of the posterior distribution to report. Defaults the 95th percentile
#' @param seed Option to set a custom seed to initialize JAGS chains, for reproducibility. Should be an integer. This argument will be deprecated in the next version, but you can always set the outside the function yourself.
#' @details There are a number of model to choose from:
#' \itemize{
#'  \item{\code{"smooth"}}{ The default option. Indicator defined by Growth rates, with Ruppert smoother, allowing for species to join late. Error on the first year of each species' time-series is assumed to be zero. The indicator is the expected value of the geometric mean across species (with missing species imputed). 
#'  Includes three options: `seFromData` `Y1perfect` and `incl.2deriv`. See bayesian_meta_analysis for mode details. Using the default values `seFromData = FALSE` and `Y1perfect = TRUE` are the options used in Freeman  \emph{et al.} (2020).}
#'  \item{\code{"smooth_det2"}}{ Equivalent to smooth with `seFromData = TRUE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  \item{\code{"smooth_det_sigtheta"}}{ Equivalent to smooth with `seFromData = FALSE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  }
#' @return Returns a dataframe with 7 columns: Year, Index.Mprime, lowerCI.Mprime, upperCI.Mprime, Index.M, lowerCI.M and, upperCI.M. 
#' Columns headed `M` and `Mprime` are means of the M and M' parameters  as defined in Freeman et al (2020). The 'upper' and 'lower' columns are the credible intervals, the width of which is defined by the `CI` argument.
#' Note that M and M' are alternate ways of calculating the multispecies indicator: their means are nearly always virtually identical, but the uncertainty in M is usually much wider than in M'. See Freeman et al (2020) for more details.
#' @import reshape2
#' @importFrom boot inv.logit
#' @importFrom coda mcmc.list as.mcmc
#' @references Freeman, S.N., Isaac, N.J.B., Besbeas, P.T., Dennis, E.B. & Morgan, B.J.T. (2020) 
#'             A generic method for estimating and smoothing multispecies biodiversity indices, robust to intermittent data. 
#'             \emph{Journal of Agricultural Biological and Environmental Statistics}, in revision.
#' @export
#' @examples 
#' # Create some example data in the format required
#' data <- data.frame(species = rep(letters, each = 50),
#'                    year = rep(1:50, length(letters)),
#'                    index = rnorm(n = 50 * length(letters), mean = 0, sd = 1),
#'                    se = runif(n = 50 * length(letters), min = 0.01, max = .1))
#' 
#' # Run the Bayesian meta-analysis
#' bma_indicator <- bma(data, model="smooth", m.scale="logit", n.iter=100)
#' 
#' # Plot the resulting indicator
#' plot_indicator(indicator = bma_indicator[,'Index.Mprime'],
#'                CIs = bma_indicator[,c(3,4)])

bma <- function (data,
              plot = TRUE,
              model = 'smooth',
              parallel = FALSE,
              n.cores = parallel::detectCores()-1,
              incl.model = TRUE,
              n.iter = 1e4,
              n.thin = 5,
              m.scale = 'loge',
              num.knots = 12,
              seFromData = FALSE,
              Y1perfect = TRUE,
              rescale_indices = NULL,
              rescaleYr = 1,
              baseline = 100,
              errorY1 = FALSE,
              save.sppars = TRUE,
              incl.2deriv = FALSE,
              CI = 95,
              seed = NULL){

set.seed(seed = seed)

# Check if jagsUI is installed
if (!requireNamespace("jagsUI", quietly = TRUE)) {
stop("Package 'jagsUI' is needed for the 'bma' function to work. Please install this from CRAN. You will also be required to install JAGS, which you can download from https://sourceforge.net/projects/mcmc-jags/files/JAGS/",
      call. = FALSE)
}

# Check to see that column names are present
if (!all(c("species", "year", "index") %in% colnames(data))){

stop('data column names should include "species", "year" and "index"')

}

if (!("se" %in% colnames(data))){

  if(seFromData) {
  stop("Error: Standard errors have not been provided")
}

data$se <- NA
}

# reconfigure function parameters based on model choice
# Print statements when parameters are reconfigured
switch(gsub(" ", "_", tolower(model)),
      smooth = {}, # the default

      smooth_jabes = {# this is the version SF ran for the paper
      cat("'smooth_jabes' model chosen. The following parameters have been set:", "model = 'smooth'", "seFromData = FALSE", "Y1perfect = TRUE", sep = "\n")
        model = "smooth"
        seFromData = FALSE
        Y1perfect = TRUE},

      smooth_det2 = {# this is the version NI tested for ISEC 
      cat("'smooth_det2' model chosen. The following parameters have been set:", "model = 'smooth'", "seFromData = TRUE", "Y1perfect = FALSE", sep = "\n")
        model = "smooth"
        seFromData = TRUE
        Y1perfect = FALSE},

      smooth_det_sigtheta = {# this is the version NI ran for Scottish indicators
      cat("'smooth_det_sigtheta' model chosen. The following parameters have been set:", "model = 'smooth'", "seFromData = FALSE", "Y1perfect = FALSE", sep = "\n")
      model = "smooth"
      seFromData=FALSE
      Y1perfect = FALSE},
    
      smooth_det = stop("smooth_det model has been deprecated"),
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
if(min(data$index, na.rm = T) >= 0){
  print("Warning: No negative index values detected. Are you sure you transformed the data?")}

# check whether the data contain any infinite values
if(any(is.infinite(data$index))){
  stop('Dataset contains Infinite values. Fix this before proceeding')}

# This is not my preferrred behaviour
if(!m.scale %in% c('loge', 'log10', 'logit')){
   stop("m.scale must be 'loge', 'log10', or 'logit'")}

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
index <- acast(data, species ~ year, value.var = "index")

if(!is.null(rescale_indices)){
# the user has specified that each species' time series should be scaled to start at a common value.
# since the data are on the unbounded (log or logit) scale, we standardise them by addition/subtraction, rather than multiplication (as in rescale_species())
# recall that rescale_indices has to be an integer
# without mising data this would be easy, but we have to identify the first year in each species' timeseries
index <- t(apply(index, 1, function(x) rescale_indices + x - x[min(which(!is.na(x)))]))
}

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
#FY <- sapply(row.names(index), FUN = function(x){
#  min(data$year[!is.na(data$index) & data$species == x])
#})
#bugs_data[['FY']] <- FY - min(FY) + 1 # set lowest value to 1
bugs_data[['FY']] <- apply(index, 1, function(x) min(which(!is.na(x)))) # simpler alternative
}

# Setup parameters to monitor. NB Most of the model options have been deprecated, so much of this code is redundant
params = c("tau.spi")
if(model == 'smooth'){ params <- c(params, "Mprime")}
if(model %in% c('smooth', 'smooth_stoch', 'smooth_det', 'FNgr', 'smooth_stoch2', 'FNgr2')){
   params <- c(params, "logLambda", "spgrowth", "M")}
if(model %in% c('smooth_stoch', 'smooth_det', 'FNgr')){params <- c(params, "tau.sg")}
if(model %in% c('smooth', 'smooth_stoch', 'smooth_det','smooth_stoch2')){params <- c(params, "beta", "taub")}
if(incl.2deriv){ params <- c(params, "t2dash")}
if(!seFromData){ params <- c(params, "theta")}
if(save.sppars) {
params <- c(params, "spindex")
} else {
params <- params[!params %in% c("spgrowth", "sigma.obs")]
}
if(incl.2deriv){ params <- c(params, "t2dash")}

set.seed(seed = seed)
model.out <- jagsUI::jags(data = bugs_data,
                        inits = NULL,
                        param = params,
                        parallel = parallel,
                        n.cores = n.cores,
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
  year_ests <- colnames(array_sim[[x]])[grepl("^Mprime\\[",  # changed from I
                                              colnames(array_sim[[x]]))]
  ar_temp <- array_sim[[x]][,c(head(year_ests, 1), tail(year_ests,1))]
  colnames(ar_temp) <- c("First year", "Last year")
  as.mcmc(ar_temp)
}, array_sim = array_sim))

plot(comb.samples)
}

MSI1 <- MSI2 <- NULL
if("Mprime" %in% params) {
MSI1 <- scale_indicator(logI = model.out$sims.list$Mprime, FirstYr = min(data$year),
                        CI=CI, m.scale=m.scale,
                        rescaleYr=rescaleYr, errorY1=errorY1, baseline=baseline)
}
if("M" %in% params) {
MSI2 <- scale_indicator(logI = model.out$sims.list$M, FirstYr = min(data$year),
                        CI=CI, m.scale=m.scale,
                        rescaleYr=rescaleYr, errorY1=errorY1, baseline=baseline)
}

MSI <- merge(MSI1, MSI2, by="Year")  # this will fail for smooth_det

# rename the output columns
names(MSI) <- gsub(names(MSI), pattern="x", replacement = "Mprime")
names(MSI) <- gsub(names(MSI), pattern="y", replacement = "M")
names(MSI) <- gsub(names(MSI), pattern="mean", replacement = "Index")

if(incl.model) attr(MSI, 'model') <- model.out

return(MSI)
}