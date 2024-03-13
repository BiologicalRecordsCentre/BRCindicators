# Fake data taken from vignette
bayesian_meta_analysis_fake_data <- data.frame(species = rep(letters, each = 50),
                   year = rep(1:50, length(letters)), 
                   index = runif(n = 50 * length(letters), min = 0, max = 1), 
                   se = runif(n = 50 * length(letters), min = 0.01, max = .1))


bma_objects <- function (data,
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

# To capture objects
output = list()

set.seed(seed = seed)

# Check if jagsUI is installed
if (!requireNamespace("jagsUI", quietly = TRUE)) {
stop("Package 'jagsUI' is needed for the 'bma' function to work. Please insatll this from CRAN. You will also be required to install JAGS, which you can download from https://sourceforge.net/projects/mcmc-jags/files/JAGS/",
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

output$comb.samples = comb.samples

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

output$MSI = MSI

return(output)
}