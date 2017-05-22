bma_model_ranwalk <- function(temp_file = tempfile()){
  
  model <- '
model {
# Priors 
  
  #alpha ~ dnorm(0, 0.001)
  #alpha.b0 ~ dnorm(0, 0.001)
  
  # one value per species
  for (i in 1:Nsp){b0[i] ~ dnorm(0, tau.b0)}
  #for (i in 1:Nsp){b0[i] ~ dnorm(alpha.b0, tau.b0)}

  # one value per year
  logI[1] ~ dnorm(0, 0.0001)
  for (t in 2:Nyr){logI[t] ~ dnorm(logI[t-1], tau.I)}
  
  # one value per site-species
  for (i in 1:Nsp){   
  for (t in 1:Nyr){
  eta[i,t] ~ dnorm(0, tau.eta) # process error
  se[i,t] ~ dunif(0,max_se) # for the missing values
  }}
  
  # Hyperpriors
  tau.b0 ~ dt(0,1,1)T(0,)
  tau.I ~ dt(0,1,1)T(0,)
  tau.eta ~ dt(0,1,1)T(0,)
  
  # Each year-species combos is estimated with error
  for (t in 1:Nyr){
  for (i in 1:Nsp){ 
  estimate[i,t] ~ dnorm(muN[i,t], tau.psd[i,t]) # Estimated effect size
  tau.psd[i,t] <- pow(se[i,t], -2)    # Known measurement error
  
  #muN is the true unknown species index this year (on the log scale)
  # its a simple linear function of the year and species effects, with "process error"
  muN[i,t] <- b0[i] + logI[t] + eta[i,t]
  # muN[i,t] <- alpha + b0[i] + logI[t] + eta[i,t]
  }}
  
  # derived parameter: indicator on the measurement scale
  #for (t in 1:Nyr){log(I[t]) <- logI[t]}
  
}'
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}


bma_model_uniform <- function(temp_file = tempfile()){
  
  model <- '
  model {
  # Priors 
  #alpha ~ dnorm(0, 0.001)
  
  # one value per species
  for (i in 1:Nsp){b0[i] ~ dnorm(0, tau.b0)}
  
  # one value per year
  for (t in 1:Nyr){logI[t] ~ dunif(-10, 10)}
  
  # one value per site-species
  for (i in 1:Nsp){   
  for (t in 1:Nyr){
  eta[i,t] ~ dnorm(0, tau.eta) # process error
  se[i,t] ~ dunif(0,max_se) # for the missing values
  }}
  
  # Hyperpriors
  tau.b0 ~ dt(0,1,1)T(0,) # redundant if species have already been standardised
  #tau.I ~ dt(0,1,1)T(0,) # obsolete in this model, but retained for compatability
  tau.eta ~ dt(0,1,1)T(0,)
  
  # Each year-species combos is estimated with error
  for (t in 1:Nyr){
  for (i in 1:Nsp){ 
  estimate[i,t] ~ dnorm(muN[i,t], tau.psd[i,t]) # Estimated effect size
  tau.psd[i,t] <- pow(se[i,t], -2)    # Known measurement error
  
  #muN is the true unknown species index this year (on the log scale)
  # its a simple linear function of the year and species effects, with "process error"
  muN[i,t] <- b0[i] + logI[t] + eta[i,t]
  #muN[i,t] <- b0[i] + logI[t]
  #muN[i,t] <- alpha + b0[i] + logI[t] + eta[i,t]
  }}

  # derived parameter: indicator on the measurement scale
  #for (t in 1:Nyr){log(I[t]) <- logI[t]}
  
  }'
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}


bma_model_uniform_noeta <- function(temp_file = tempfile()){
  
  model <- '
  model {
  # Priors 
  alpha ~ dnorm(0, 0.001)
  
  # one value per species
  for (i in 1:Nsp){b0[i] ~ dnorm(0, tau.b0)}
  
  # one value per year
  for (t in 1:Nyr){logI[t] ~ dunif(-10, 10)}
  
  # one value per site-species
  for (i in 1:Nsp){   
  for (t in 1:Nyr){
  eta[i,t] ~ dnorm(0, tau.eta) # process error
  se[i,t] ~ dunif(0,max_se) # for the missing values
  }}
  
  # Hyperpriors
  tau.b0 ~ dt(0,1,1)T(0,) # redundant if species have already been standardised
  tau.I ~ dt(0,1,1)T(0,) # obsolete in this model, but retained for compatability
  tau.eta ~ dt(0,1,1)T(0,)
  
  # Each year-species combos is estimated with error
  for (t in 1:Nyr){
  for (i in 1:Nsp){ 
  estimate[i,t] ~ dnorm(muN[i,t], tau.psd[i,t]) # Estimated effect size
  tau.psd[i,t] <- pow(se[i,t], -2)    # Known measurement error
  
  #muN is the true unknown species index this year (on the log scale)
  # its a simple linear function of the year and species effects, with "process error"
  #muN[i,t] <- b0[i] + logI[t] + eta[i,t]
  muN[i,t] <- b0[i] + logI[t]
  #muN[i,t] <- alpha + b0[i] + logI[t] + eta[i,t]
  }}

  # derived parameter: indicator on the measurement scale
  #for (t in 1:Nyr){log(I[t]) <- logI[t]}

  }'
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}


bma_model_test <- function(temp_file = tempfile()){
  
  model <- '
  model {
  # Priors 
  
  #alpha ~ dnorm(0, 0.001)
  #alpha.b0 ~ dnorm(0, 0.001)
  
  # one value per species
  #for (i in 1:Nsp){b0[i] ~ dnorm(0, tau.b0)}
  for (i in 1:Nsp){b0[i] ~ dnorm(alpha.b0, tau.b0)}
  b0[i] ~ 0.5
  
  # one value per year
  logI[1] ~ dunif(-10, 10)
  for (t in 2:Nyr){logI[t] ~ dnorm(logI[t-1], tau.I)}
  
  # one value per site-species
  for (i in 1:Nsp){   
  for (t in 1:Nyr){
  #eta[i,t] ~ dnorm(0, tau.eta) # process error
  se[i,t] ~ dunif(0,max_se) # for the missing values
  }}
  
  # Hyperpriors
  tau.b0 ~ dt(0,1,1)T(0,)
  tau.I ~ dt(0,1,1)T(0,)
  #tau.eta ~ dt(0,1,1)T(0,)
  
  # Each year-species combos is estimated with error
  for (t in 1:Nyr){
    for (i in 1:Nsp){ 
      estimate[i,t] ~ dnorm(muN[i,t], tau.psd[i,t]) # Estimated effect size
      tau.psd[i,t] <- pow(se[i,t], -2)    # Known measurement error
  
      #muN is the true unknown species index this year (on the log scale)
      # its a simple linear function of the year and species effects, with "process error"
      muN[i,t] <- b0[i] + logI[t]
      #muN[i,t] <- b0[i] + logI[t] + eta[i,t]
      #muN[i,t] <- alpha + b0[i] + logI[t] + eta[i,t]
    }
  }
  
  # derived parameter: indicator on the measurement scale
  #for (t in 1:Nyr){
  #  combi_logI[t] <- alpha + logI[t]
  #}
  
  }'
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}

