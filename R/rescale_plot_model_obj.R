#' @export
#' 
rescale_plot_model_obj <- function(){

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