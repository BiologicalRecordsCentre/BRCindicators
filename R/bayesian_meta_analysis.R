#' Bayesian Meta-analysis
#' 
#' @description Use a bayesian meta-analysis to create an indicator from species index values and standard error.
#' 
#' @param data a data.frame with 4 columns in this order: species, year, index, se (standard error) 
#' @param randwalk If \code{TRUE} a random walk feature is added to the model (defualt), else years are treated as independent.
#' @param plot Logical, should a trace plot be plotted?
#' @return Returns a dataframe with 4 columns: Year, Index, lower2.5, upper97.5. The last two columns are the credible intervals
#' @import R2jags
#' @import reshape2
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


bma <- function(data, randwalk = TRUE, plot = TRUE) {

  if(!identical(colnames(data), c('species', 'year', 'index', 'se'))){
    stop('data column names should be: "species", "year", "index", "se"')
  }
  
  # Create 2 matricies, for index and se
  index <- log(acast(na.omit(data), species ~ year, value.var='index'))
  se <- acast(na.omit(data), species ~ year, value.var='se')

  # convert to a format for bugs
  bugs_data <- list(Nsp = nrow(index),
                    Nyr = ncol(index),
                    estimate = index,
                    se =  se)

  #setup bugs
  params = c('tau.b0', 'tau.eta', 'tau.I', 'I', 'alpha')
  
  # run model demending on user specification
  if(randwalk){
    model <- jags(data = bugs_data,
                  inits = NULL,
                  param = params, 
                  model.file = bma_indicator(),
                  n.chains = 3,
                  n.thin = 2,
                  n.iter = 10000)
  } else {
    model <- jags(data = bugs_data,
                  inits = NULL,
                  param = params, 
                  model.file = bma_indicator_randwalk(),
                  n.chains = 3,
                  n.thin = 2,
                  n.iter = 10000)
  }

  
  if(plot){
    array_sim <- model$BUGSoutput$sims.array
    comb.samples <- coda::mcmc.list(
      lapply(1:3,
             FUN = function(x, array_sim){
               year_ests <- colnames(array_sim[,x,])[grepl('^I\\[', colnames(array_sim[,x,]))]
               ar_temp <- array_sim[,x,c(head(year_ests, 1), tail(year_ests, 1))]
               colnames(ar_temp) <- c('First year', 'Last year')
               coda::as.mcmc(ar_temp)
               },
             array_sim = array_sim)
    )
    plot(comb.samples)
  }
  
  
  pd <- model$BUGSoutput$summary[,c(1,3,7)]
  
  # take the raw output from model and convert it into an indicator,
  # with options for centering [not implemented]
  rescale_bayesian_indicator <- function(x, centering = 'firstyr') {
    if(centering == 'firstyr'){
      x <- x/x[1,1]
    } 
    100 * x
  }

  pd <- rescale_bayesian_indicator(pd, centering='firstyr')

  pd <- dcast(melt(pd[grepl('^I',dimnames(pd)[[1]]),]), Var1~Var2)
  names(pd) <- c('Year', 'Index', 'lower2.5', 'upper97.5')
  pd$Year  <- as.numeric(pd$Year)

  return(pd)

}
