#' @export

# Multispecies indicators from Bayesian Meta-Analysis 
# Each function contains BUGS code for a different model
# Code by Nick Isaac, Tom August & Steve Freeman

################################################################################

# 15/10/19: We're streamlining the options for the BMA. 
      # functions now return the model text, rather than writing a temp file
      # I'm including the deprecated options for now.
# 15/11/17: There are six models, all with the prefix "bma_model"
  # the ranwalk, uniform & uniform_noeta were written by Nick Isaac
  # FNgr, SmoothStoch & SmoothDet were written by Stephen Freeman
  # Isaac and Freeman used different syntax: I have harmonized some terms
    # Terms in Freeman's code that have been changed to match syntax in Isaac's code:
      # logI is used for the multispecies indicator on the log scale (in place of 'logI2' or 'sindicator')
      # logLambda is used for the multispecies log growth rate (in place of 'growth')
      # estimate is used for the estimated log abundance (in place of 'species')
      # tau.sg used instead of taus
    # Terms in Isaac's code that have been changed to match syntax in Freeman's code:
      # spindex used for the true unknown state variable (in place of "muN")
      # s used for species loop (in place of "i")
      # nyr and nsp used in place of Nyr & Nsp
    # Terms changed from both
      # sigma.obs used for the "known" measurement error (formerly sigma in Freeman, se in Isaac)
      # tau.obs used for 1/sigma.obs^2 (formerly tau in Freeman, tau.psd in Isaac)
    # Other changes
      # sigma.obs is allowed to be data, but can be omitted and imputed from prior (as in Isaac's code)
      # variances (sigma2) removed from Freeman code

    # outstanding issues: 
    # Isaac code uses half Cauchy priors for precisions, Freeman uses uniform on stdvs
    # Isaac code does not calculate stdvs from precisions
# TA additions
  # tau.obs changed to tau.obs[s,t] in new 3 models [errors otherwise]
  # undocumented change tau.b0 -> tau.spi

################################################################################

bma_model_Smooth <- function(incl.2deriv = FALSE){

  priors<- '
    ###################  Define variance priors ###########################
  
    # process errors
    sigma.spi ~ dunif(0,1000)
    sigma2.spi <- pow(sigma.spi,2)
    tau.spi <- pow(sigma.spi,-2)
    
    # observation errors
    theta ~ dunif(0,10) # observation error is constant

  '

  smoothing <- '
    ######################### Smoothing done here   #######################
 
    beta[1] ~ dnorm(0, 0.000001)
    beta[2] ~ dnorm(0, 0.000001)
    taub ~ dgamma(0.000001, 0.000001)
    for(k in 1:num.knots){b[k]~dnorm(0,taub)}
  
    for (t in 1:(nyears - 1)){
      logLambda[t] <- m[t]
      m[t] <- mfe[t] + mre[t]
      mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
      for (k in 1:num.knots){
        temp[t,k] <- b[k] * Z[t,k]
      }
      mre[t]<-sum(temp[t,1:num.knots])
    }  
    
  '
  
  likelihood <- '
  #######################  Define likelihood  ###########################
  
  logI2[1]<-0 
  
  for (t in 2:nyears){
    logI2[t]<-logI2[t-1] + logLambda[t-1]
  }
  
  for (s in 1:nsp){
    for (t in 1:(nyears-1)){
      spgrowth[s,t] ~ dnorm(logLambda[t], tau.spi)
    }}
  
  for (s in 1:nsp){
    for (t in 1:(FY[s]-1){
      spindex[s,t] <- spindex[s,t+1] - spgrowth[s,t]
    }}

  #spindex[s,FY[s]] ~ dnorm(0, 0.0001) + estimate[s,FY[s]]
  spindex[s,FY[s]] <- estimate[s,FY[s]] # we assume the first year is known without error

    for (t in (FY[s]+1):(nyears)){
      spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
      tau.obs[s,t] <- pow(theta, -2)
  }}
  
  ####################  geomean of expected values ######################
  
  for (t in 1:nyears){
    logI[t] <- sum(spindex[,t])/nsp
  }

  '
  
  derivatives <- ifelse(incl.2deriv, "", {'
  
  #########################  second derivatives #######################
  
  I <- logI2
  t2dash[2]<-(I[2+1] - 2*I[2] + I[2-1])/1
  t2dash[nyears-1] <- (I[nyears] - 2*I[nyears-1] + I[nyears-2])/1
  t2dash[3]<-(-I[5]+16*I[4]-30*I[3]+16*I[2]-I[1] )/12
  t2dash[nyears-2] <- (-I[nyears]+16*I[nyears-1]-30*I[nyears-2]+16*I[nyears-3]-I[nyears-4])/12
  for (t in 4:(nyears-3)){
    t2dash[t]<-(2*I[t+3]-27*I[t+2]+270*I[t+1]-490*I[t]+270*I[t-1]-27*I[t-2]+2*I[t-3])/180
  }
  
  #####################################################################
  
  '})
  
  model <- paste(c("model {",
                 priors, 
                 smoothing, 
                 likelihood,
                 derivatives,
                 "}"), collapse = "\n")

  return(model)
}

################################################################################
# BEGIN DEPRECATED OPTIONS
################################################################################

bma_model_ranwalk <- function(){
  # Also known as BMA3
  
  model <- '
model {
  ###################  Define priors
  
  # one value per species
  for (s in 1:nsp){b0[s] ~ dnorm(0, tau.spi)}
  
  # one value per year
  logI[1] ~ dnorm(0, 0.0001)
  for (t in 2:nyears){logI[t] ~ dnorm(logI[t-1], tau.I)}
  
  # one value per site-species
  for (s in 1:nsp){   
    for (t in 1:nyears){
      eta[s,t] ~ dnorm(0, tau.eta) # process error
      sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  }}
  
  # Hyperpriors
  tau.spi ~ dt(0,1,1)T(0,)
  tau.I ~ dt(0,1,1)T(0,)
  tau.eta ~ dt(0,1,1)T(0,)

  ###################  Define likelihood  #######################

  # Each year-species combos is estimated with error
  for (t in 1:nyears){
    for (s in 1:nsp){ 
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t]) # Estimated effect size
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)    # Known measurement error
  
  #spindex is the true unknown species index this year (on the log scale)
  # its a simple linear function of the year and species effects, with "process error"
  spindex[s,t] <- b0[s] + logI[t] + eta[s,t]

  }}

  #########################  end likelihood ###########################

}'
  return(model)
}

################################################################################

bma_model_uniform <- function(){
  # Also known as BMA2
  
  model <- '
  model {
  ###################  Define priors

  # one value per species
  for (s in 1:nsp){b0[s] ~ dnorm(0, tau.spi)}
  
  # one value per year
  for (t in 1:nyears){logI[t] ~ dunif(-10, 10)}
  
  # one value per site-species
  for (s in 1:nsp){   
    for (t in 1:nyears){
      eta[s,t] ~ dnorm(0, tau.eta) # process error
      sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  }}
  
  # Hyperpriors
  tau.spi ~ dt(0,1,1)T(0,) # redundant if species have already been standardised
  tau.eta ~ dt(0,1,1)T(0,)

  ###################  Define likelihood  #######################
  
  # Each year-species combos is estimated with error
  for (t in 1:nyears){
    for (s in 1:nsp){ 
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t]) # Estimated effect size
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)    # Known measurement error
  
  #spindex is the true unknown species index this year (on the log scale)
  # its a simple linear function of the year and species effects, with "process error"
  spindex[s,t] <- b0[s] + logI[t] + eta[s,t]
  }}

  #########################  end likelihood ###########################

  }'
  return(model)
}

################################################################################

bma_model_uniform_noeta <- function(){
  # Also known as BMA1
  
  model <- '
  model {
  ###################  Define priors

  # one value per species
  for (s in 1:nsp){b0[s] ~ dnorm(0, tau.spi)}
  
  # one value per year
  for (t in 1:nyears){logI[t] ~ dunif(-10, 10)}
  
  # one value per site-species
  for (s in 1:nsp){   
    for (t in 1:nyears){
      eta[s,t] ~ dnorm(0, tau.eta) # process error
      sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  }}
  
  # Hyperpriors
  tau.spi ~ dt(0,1,1)T(0,) # redundant if species have already been standardised
  tau.eta ~ dt(0,1,1)T(0,)

  ###################  Define likelihood  #######################

  # Each year-species combos is estimated with error
  for (t in 1:nyears){
    for (s in 1:nsp){ 
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t]) # Estimated effect size
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)    # Known measurement error
  
  #spindex is the true unknown species index this year (on the log scale)
  # its a simple linear function of the year and species effects, with "process error"
  spindex[s,t] <- b0[s] + logI[t]
  }}

  #########################  end likelihood ###########################

  }'
  return(model)
}

################################################################################

bma_model_FNgr <- function(){
  # Indicator defined by Growth rates, as in Freeman & Newson
  # Defined by equations 1-5 in Steve Freeman's document of 2/11/17
  # Also known as "FN_indicator" in Steve's email of 2/11/17
  
  model <- '
  model {

  ###################  Define priors
  # process errors
  tau.spi <- pow(sigma.spi,-2)
  sigma.spi ~ dunif(0,30)
  tau.sg <- pow(sigma.sg,-2)
  sigma.sg ~ dunif(0,1000)

  # observation errors
  # one value per site-species
  for (s in 1:nsp){   
    for (t in 1:nyears){
      sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
    }}

  for (s in 1:nsp){
    spindex[s,1] ~ dnorm(logI[1], tau.spi)
  }
  
  for (t in 1:(nyears)){
    logI[t] ~ dnorm(0,0.000001)
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 1:(nyears-1)){logLambda[t] <- logI[t+1] - logI[t]}
  
  for (s in 1:nsp){
    for (t in 1:(nyears-1)){
      spgrowth[s,t] ~ dnorm(logLambda[t], tau.sg)
  }}
  
  for (s in 1:nsp){
    for (t in 2:(nyears)){
      spindex[s,t] <- spindex[s,1] + sum(spgrowth[s,1:(t-1)])
  }}
  
  for (s in 1:nsp){
    for (t in 1:(nyears)){
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
  }}
  

  
  #########################  end likelihood ###########################

}'
  return(model)
}

################################################################################

bma_model_smooth_stoch <- function(){
  # Indicator defined by Growth rates, with Ruppert smoother (stochastic version)
  # Defined by equation 6 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_2" in Steve's email of 2/11/17
  
  model <- '
  model {

  ###################  Define priors
  # process errors
  tau.spi <- pow(sigma.spi,-2)
  sigma.spi ~ dunif(0,30)
  tau.sg <- pow(sigma.sg,-2)
  sigma.sg ~ dunif(0,1000)

  # observation errors
  # one value per site-species
  for (s in 1:nsp){   
    for (t in 1:nyears){
      sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  }}
  
  for (s in 1:nsp){
    spindex[s,1] ~ dnorm(logI.raw[1], tau.spi)
  }
  
  ########### Smoothing done here   #############
  
  beta[1]~dnorm(0,0.000001)
  beta[2]~dnorm(0,0.000001)
  taueps ~ dgamma(0.000001,0.000001)
  taub ~ dgamma(0.000001,0.000001)
  for(k in 1:num.knots){b[k] ~ dnorm(0,taub)}
  
  for (t in 1:(nyears)){
    logI.raw[t] ~ dnorm(m[t], taueps)
    m[t] <- mfe[t] + mre[t]
    mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
    for (k in 1:num.knots){
      temp[t,k]<-b[k]*Z[t,k]
    }
    mre[t] <- sum(temp[t,1:num.knots])
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 1:(nyears-1)){
  logLambda[t] <- logI.raw[t+1] - logI.raw[t]
  }
  
  for (s in 1:nsp){
    for (t in 1:(nyears-1)){
      spgrowth[s,t] ~ dnorm(logLambda[t], tau.sg)
  }}
  
  for (s in 1:nsp){
    for (t in 2:(nyears)){
      spindex[s,t] <- spindex[s,1] + sum(spgrowth[s,1:(t-1)])
  }}
  
  for (s in 1:nsp){
    for (t in 1:(nyears)){
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
  }}
  
  logI <- m
  
  #########################  end likelihood ###########################
  
  
  }'
  return(model)
}

################################################################################

bma_model_smooth_det <- function(){
  # Indicator defined by Growth rates, with Ruppert smoother (deterministic version)
  # Defined by equation 7 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_1" in Steve's email of 2/11/17
  # this version takes standard errors on year 1 estimates. 
  # If actually zero then the model invents them (this is a fudge)

  model <- '
  model {

  ###################  Define priors
  # process errors
  tau.spi<-pow(sigma.spi,-2)
  sigma.spi~dunif(0,30)
  tau.sg<-pow(sigma.spi,-2)
  sigma.spi~dunif(0,1000)

  # observation errors
  # one value per site-species
  for (s in 1:nsp){   
    for (t in 1:nyears){
      sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
    }}
  
  for (s in 1:nsp){
    spindex[s,1]~dnorm(logI.raw[1],tau.spi)
  }
  
  ########### Smoothing done here   #############
  
  beta[1]~dnorm(0,0.000001)
  beta[2]~dnorm(0,0.000001)
  taub~dgamma(0.000001,0.000001)
  for(k in 1:num.knots){b[k]~dnorm(0,taub)}
  
  for (t in 1:(nyears)){
    logI.raw[t] <- m[t]
    m[t] <- mfe[t]+mre[t]
    mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
    for (k in 1:num.knots){
      temp[t,k] <- b[k]*Z[t,k]
      }
    mre[t] <- sum(temp[t,1:num.knots])
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 1:(nyears-1)){
    logLambda[t] <- logI.raw[t+1] - logI.raw[t]
  }
  
  for (s in 1:nsp){
    for (t in 1:(nyears-1)){
      spgrowth[s,t]~dnorm(logLambda[t],tau.sg)
    }
  }
  
  for (s in 1:nsp){
    for (t in 2:(nyears)){
      spindex[s,t]<-spindex[s,1] + sum(spgrowth[s,1:(t-1)])
    }
  }
  
  for (s in 1:nsp){
    for (t in 1:(nyears)){
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
    }
  }
  
  logI <- m
  
  #########################  end likelihood ###########################

}
  '
  return(model)
}

################################################################################

bma_model_FNgr2 <- function(){
  # Indicator defined by Growth rates, as in Freeman & Newson
  # Defined by equations 1-5 in Steve Freeman's document of 2/11/17
  # Also known as "FN_indicator" in Steve's email of 2/11/17
  # 29/11/17 Now takes an extra vector, FY, indexing the first year for which a species has data
  # this makes some other code redundant
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species are allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in logI2[1] doesn't measure the same thing as uncertainty in other years 
  # tau.spi is on the growth rates, not the index
  # logI is now estimated without uncertainty due to interspecific variation
  
  model <- '
  model {
  
  ###################  Define priors
  # process errors
  tau.spi <- pow(sigma.spi,-2)
  sigma.spi ~ dunif(0,30)
  
  # observation errors
  # one value per site-species
  for (s in 1:nsp){
   for (t in 1:nyears){
     sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  }} 

  logI2[1] <- 0
  for (t in 2:(nyears)){
    logI2[t] ~ dnorm(0,0.000001)
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 1:(nyears-1)){
    logLambda[t] <- logI2[t+1] - logI2[t]
  }
  
  for (s in 1:nsp){
  for (t in 1:(nyears-1)){
    spgrowth[s,t] ~ dnorm(logLambda[t], tau.spi)
  }}
  
  for (s in 1:nsp){
   for (t in 1:FY[s]){
    spindex[s,t] <- spindex[s,t+1] - spgrowth[s,t]
    }
   for (t in (FY[s]+1):(nyears)){
    spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
    estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
    tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
  }}
  
  #mean of the species indices  
  for (t in 1:nyears) {
    logI[t] <- sum(spindex[,t])/nsp
  }

  #########################  end likelihood ###########################
  }'
  return(model)
}

################################################################################

bma_model_smooth_stoch2 <- function(){
  # Indicator defined by Growth rates, with Ruppert smoother (stochastic version)
  # Defined by equation 6 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_2" in Steve's email of 2/11/17
  # 29/11/17 Now takes an extra vector, FY, indexing the first year for which a species has data
  # this makes some other code redundant
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species area allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in logI2[1] doesn't measure the same thing as uncertainty in other years 
  # tau.spi is on the growth rates, not the index
  # logI is now estimated without uncertainty due to interspecific variation
  
  model <- '
  model {
  
  ###################  Define priors
  # process errors
  tau.spi <- pow(sigma.spi,-2)
  sigma.spi ~ dunif(0,30)
 
  # observation errors
  # one value per site-species
  for (s in 1:nsp){
   for (t in 1:nyears){
     sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  }}
  
  logI2[1] <- 0

  ########### Smoothing done here   #############
  
  beta[1]~dnorm(0,0.000001)
  beta[2]~dnorm(0,0.000001)
  taueps ~ dgamma(0.000001,0.000001)
  taub ~ dgamma(0.000001,0.000001)
  for(k in 1:num.knots){b[k] ~ dnorm(0,taub)}
  
  for (t in 1:(nyears)){
  logLambda[t] ~ dnorm(m[t], taueps)
  m[t] <- mfe[t] + mre[t]
  mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
  for (k in 1:num.knots){
  temp[t,k]<-b[k]*Z[t,k]
  }
  mre[t] <- sum(temp[t,1:num.knots])
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 2:nyears){
    logI2[t] <- logI2[t-1] + logLambda[t-1]
  }
  
  for (s in 1:nsp){
  for (t in 1:(nyears-1)){
    spgrowth[s,t] ~ dnorm(logLambda[t], tau.spi)
  }}
  
   for (s in 1:nsp){
     for (t in 1:FY[s]){
    spindex[s,t] <- spindex[s,t+1] - spgrowth[s,t]
  }
    for (t in (FY[s]+1):(nyears)){
      spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
  }}

  #indicator is mean of the estimated species indices  
  for (t in 1:nyears) {
    logI[t] <- sum(spindex[,t])/nsp
  }


  #########################  end likelihood ###########################
  
  }'
  return(model)
}

################################################################################

bma_model_smooth_det2 <- function(){
  # Indicator defined by Growth rates, with Ruppert smoother (deterministic version)
  # Defined by equation 7 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_1" in Steve's email of 2/11/17
  # 29/11/17 Now takes an extra vector, FY, indexing the first year for which a species has data
  # this makes some other code redundant
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species area allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in logI2[1] doesn't measure the same thing as uncertainty in other years 
  # tau.spi is on the growth rates, not the index
  # logI is now estimated without uncertainty due to interspecific variation
  
  model <- '
  model {
  
  ###################  Define priors
  # process errors
  tau.spi <- pow(sigma.spi,-2)
  sigma.spi ~ dunif(0,30)

  # observation errors
  # one value per site-species
  for (s in 1:nsp){
   for (t in 1:nyears){
    sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  }}
  
  logI2[1] <- 0

  ########### Smoothing done here   #############
  
  beta[1]~dnorm(0,0.000001)
  beta[2]~dnorm(0,0.000001)
  taub~dgamma(0.000001,0.000001)
  for(k in 1:num.knots){b[k]~dnorm(0,taub)}
  
  for (t in 1:(nyears)){
    logLambda[t] <- m[t]
    m[t] <- mfe[t]+mre[t]
    mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
    for (k in 1:num.knots){
      temp[t,k] <- b[k]*Z[t,k]
    }
    mre[t] <- sum(temp[t,1:num.knots])
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 2:nyears){
    logI2[t] <- logI2[t-1] + logLambda[t-1]
  }
  
  for (s in 1:nsp){
  for (t in 1:(nyears-1)){
    spgrowth[s,t] ~ dnorm(logLambda[t],tau.spi)
  }}
  
   for (s in 1:nsp){
     for (t in 1:FY[s]){
    spindex[s,t] <- spindex[s,t+1] - spgrowth[s,t]
  }
    for (t in (FY[s]+1):(nyears)){
      spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
      tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
    }}

  # mean of the species indices  
  for (t in 1:nyears) {
    logI[t] <- sum(spindex[,t])/nsp
  }

  #########################  end likelihood ###########################
  
  }
  '
  return(model)
}

################################################################################

bma_model_smooth_det_sigtheta <- function(){
  # Indicator defined by Growth rates, with Ruppert smoother (deterministic version)
  # Defined by equation 7 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_1" in Steve's email of 2/11/17
  # 29/11/17 Now takes an extra vector, FY, indexing the first year for which a species has data
  # this makes some other code redundant
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species area allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in logI2[1] doesn't measure the same thing as uncertainty in other years 
  # tau.spi is on the growth rates, not the index
  # logI is now estimated without uncertainty due to interspecific variation
  # in this version, standard errors are not read from the data but rather estimated
  # (but they still need to be included in the dataset)
  
  model <- '
  model {
  
  ###################  Define priors
  # process errors
  tau.spi <- pow(sigma.spi,-2)
  sigma.spi ~ dunif(0,30)
  theta ~ dunif(0,30) # observation error is constant
  
  logI2[1] <- 0

  ########### Smoothing done here   #############
  
  beta[1]~dnorm(0,0.000001)
  beta[2]~dnorm(0,0.000001)
  taub~dgamma(0.000001,0.000001)
  for(k in 1:num.knots){b[k]~dnorm(0,taub)}
  
  for (t in 1:(nyears)){
    logLambda[t] <- m[t]
    m[t] <- mfe[t]+mre[t]
    mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
    for (k in 1:num.knots){
      temp[t,k] <- b[k]*Z[t,k]
    }
    mre[t] <- sum(temp[t,1:num.knots])
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 2:nyears){
    logI2[t] <- logI2[t-1] + logLambda[t-1]
  }
  
  for (s in 1:nsp){
  for (t in 1:(nyears-1)){
    spgrowth[s,t] ~ dnorm(logLambda[t],tau.spi)
  }}
  
   for (s in 1:nsp){
     for (t in 1:FY[s]){
    spindex[s,t] <- spindex[s,t+1] - spgrowth[s,t]
  }
    for (t in (FY[s]+1):(nyears)){
      spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
      tau.obs[s,t] <- pow(theta, -2)
    }}

  # mean of the species indices  
  for (t in 1:nyears) {
    logI[t] <- sum(spindex[,t])/nsp
  }

  #########################  end likelihood ###########################
  
  }
  '
  return(model)
}

################################################################################
