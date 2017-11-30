
# Multispecies indicators from Bayesian Meta-Analysis 
# Each function contains BUGS code for a different model
# Code by Nick Isaac, Tom August & Steve Freeman

################################################################################

# 15/11/17: There are six models, all with the prefix "bma_model"
  # the ranwalk, uniform & uniform_noeta were written by Nick Isaac
  # FNgr, SmoothStoch & SmoothDet were written by Stephen Freeman
  # Isaac and Freeman used different syntax: I have harmonized some terms
    # Terms in Freeman's code that have been changed to match syntax in Isaac's code:
      # logI is used for the multispecies indicator on the log scale (in place of 'tindicator' or 'sindicator')
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

bma_model_ranwalk <- function(temp_file = tempfile()){
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
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_uniform <- function(temp_file = tempfile()){
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
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_uniform_noeta <- function(temp_file = tempfile()){
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
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_FNgr <- function(temp_file = tempfile()){
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
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_smooth_stoch <- function(temp_file = tempfile()){
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
  
  
}

  '
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_smooth_det <- function(temp_file = tempfile()){
  # Indicator defined by Growth rates, with Ruppert smoother (deterministic version)
  # Defined by equation 7 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_1" in Steve's email of 2/11/17

  model <- '
  model {

  ###################  Define priors
  # process errors
  tau.spi<-pow(sigma.y,-2)
  sigma.y~dunif(0,30)
  tau.sg<-pow(sigma.s,-2)
  sigma.s~dunif(0,1000)

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
  #taueps~dgamma(0.000001,0.000001)
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
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_FNgr2 <- function(temp_file = tempfile()){
  # Indicator defined by Growth rates, as in Freeman & Newson
  # Defined by equations 1-5 in Steve Freeman's document of 2/11/17
  # Also known as "FN_indicator" in Steve's email of 2/11/17
  # 29/11/17 Now takes an extra vector, FY, indexing the first year for which a species has data
  # this makes some other code redundant
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species are allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in logI[1] doesn't measure the same thing as uncertainty in other years 
  
  
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
  #for (s in 1:nsp){   
  #  for (t in 1:nyears){
  #    sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  #  }}
  
  #for (s in 1:nsp){
  #  spindex[s,1] ~ dnorm(logI[1], tau.spi)
  #}
  
  for (t in 1:(nyears)){
  logI[t] ~ dnorm(0,0.000001)
  }
  
  ###################  Define likelihood  #######################
  
  for (t in 1:(nyears-1)){
  logLambda[t] <- logI[t+1] - logI[t]
  }
  
  for (s in 1:nsp){
  for (t in 1:(nyears-1)){
  spgrowth[s,t] ~ dnorm(logLambda[t], tau.sg)
  }}
  
  for (s in 1:nsp){
  for (t in (FY[s]+1):(nyears)){
  spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
  estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
  tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
  }}
  
  #########################  end likelihood ###########################
  
  }'
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_smooth_stoch2 <- function(temp_file = tempfile()){
  # Indicator defined by Growth rates, with Ruppert smoother (stochastic version)
  # Defined by equation 6 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_2" in Steve's email of 2/11/17
  # 29/11/17 Now takes an extra vector, FY, indexing the first year for which a species has data
  # this makes some other code redundant
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species area allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in logI[1] doesn't measure the same thing as uncertainty in other years 
  
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
  #for (s in 1:nsp){   
  #  for (t in 1:nyears){
  #    sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  #}}
  
  #for (s in 1:nsp){
  #  spindex[s,1] ~ dnorm(logI.raw[1], tau.spi)
  #}
  logI.raw[1]~dnorm(0,0.00001) # I think this is redundant
  
  ########### Smoothing done here   #############
  
  beta[1]~dnorm(0,0.000001)
  beta[2]~dnorm(0,0.000001)
  taueps ~ dgamma(0.000001,0.000001)
  taub ~ dgamma(0.000001,0.000001)
  for(k in 1:num.knots){b[k] ~ dnorm(0,taub)}
  
  for (t in 1:(nyears)){
  #logI.raw[t] ~ dnorm(m[t], taueps)
  logLambda[t] ~ dnorm(m[t], taueps)
  m[t] <- mfe[t] + mre[t]
  mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
  for (k in 1:num.knots){
  temp[t,k]<-b[k]*Z[t,k]
  }
  mre[t] <- sum(temp[t,1:num.knots])
  }
  
  ###################  Define likelihood  #######################
  
  #for (t in 1:(nyears-1)){
  #  logLambda[t] <- logI.raw[t+1] - logI.raw[t]
  #}
  for (t in 2:nyears){
  logI.raw[t] <- logI.raw[t-1] + logLambda[t-1]
  }
  
  for (s in 1:nsp){
  for (t in 1:(nyears-1)){
  spgrowth[s,t] ~ dnorm(logLambda[t], tau.sg)
  }}
  
  for (s in 1:nsp){
  for (t in (FY[s]+1):(nyears)){
  spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
  estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
  tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
  }}
  
  logI <- m
  
  #########################  end likelihood ###########################
  
  }'
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################

bma_model_smooth_det2 <- function(temp_file = tempfile()){
  # Indicator defined by Growth rates, with Ruppert smoother (deterministic version)
  # Defined by equation 7 in Steve Freeman's document of 2/11/17
  # Also known as "smooth_indicator_1" in Steve's email of 2/11/17
  # 29/11/17 Now takes an extra vector, FY, indexing the first year for which a species has data
  # this makes some other code redundant
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species area allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in logI[1] doesn't measure the same thing as uncertainty in other years 
  
  model <- '
  model {
  
  ###################  Define priors
  # process errors
  tau.spi<-pow(sigma.y,-2)
  sigma.y~dunif(0,30)
  tau.sg<-pow(sigma.s,-2)
  sigma.s~dunif(0,1000)
  
  # observation errors
  # one value per site-species
  #for (s in 1:nsp){   
  #  for (t in 1:nyears){
  #    sigma.obs[s,t] ~ dunif(0, max_se) # for the missing values
  #  }}
  
  #for (s in 1:nsp){
  #  spindex[s,1] ~ dnorm(logI.raw[1],tau.spi)
  #}
  logI.raw[1]~dnorm(0,0.00001) # I think this is redundant
  
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
  logI.raw[t] <- logI.raw[t-1] + logLambda[t-1]
  }
  
  for (s in 1:nsp){
  for (t in 1:(nyears-1)){
  spgrowth[s,t] ~ dnorm(logLambda[t],tau.sg)
  }}
  
  for (s in 1:nsp){
  for (t in (FY[s]+1):(nyears)){
  spindex[s,t]<- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)])
  estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
  tau.obs[s,t] <- pow(sigma.obs[s,t], -2)
  }}
  
  logI <- m
  
  #########################  end likelihood ###########################
  
  }
  '
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

################################################################################
