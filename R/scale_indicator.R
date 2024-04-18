scale_indicator <- function(logI, 
                        FirstYr,
                        CI = 95,
                        m.scale = "loge",
                        rescaleYr = 1,
                        errorY1 = FALSE,
                        baseline = 100){ # SCALING FUNCTION

# first check the value of rescale year
# if its a large value fix to the final year.
if(rescaleYr > nrow(logI)) rescaleYr <- ncol(logI)
if(rescaleYr < 1) stop('rescaleYr must be an integer of 1 or higher')  

# which is the baseline year? rescale year = 0 indicates no scaling
# calculate the difference of each index value from the baseline year (default is year 1)
# this makes the baseline equal to zero, so it will be exactly one after back transformation
if(errorY1){ # we retain error in the first year
logI_rescaled <- t(apply(logI, 1, function(x) x - mean(logI[,rescaleYr])))
} else { # scale every iteration of the posterior to have the same value in the reference year
logI_rescaled <- t(apply(logI, 1, function(x) x - x[rescaleYr])) 
}

# convert the CI into quantiles
# first check that CI is sensible
if((CI > 100) | (CI <= 0)) stop("Credible intervals must be between 0 and 100")
CI2q <- function(CI) {
q <- (1 - CI/100)/2
return(c(q, 1-q))
}
q <- CI2q(CI)

# summarise quantiles of the posterior distribution
pd <- data.frame(mean = apply(logI_rescaled, 2, mean),
                lowerCI = apply(logI_rescaled, 2, quantile, probs = q[1]),
                upperCI = apply(logI_rescaled, 2, quantile, probs = q[2]),
                row.names = paste0('Mprime', 1:ncol(logI_rescaled)))
if(q[1] == 0.025 & q[2] == 0.975) names(pd)[2:3] <- c("q2.5", "q97.5")

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
#  pd <- unb2b(pd[grepl("^Mprime[[:digit:]]+",dimnames(pd)[[1]]),], m.scale) # PROBLEM HERE
pd <- unb2b(pd, m.scale)
# the indicator metric is now back on the measurement scale

# rescale to baseline value in the first year
pd <- pd * baseline

pd$Year <- as.numeric(gsub("[A-z]", repl="", dimnames(pd)[[1]]))
pd$Year <- as.numeric(pd$Year) + FirstYr - 1
return(pd)
}