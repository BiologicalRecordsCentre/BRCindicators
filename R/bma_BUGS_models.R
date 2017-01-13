bma_indicator <- function(temp_file = tempfile()){
  
  model <- "
model {
  # Priors 
  alpha ~ dnorm(0, 0.0001)
  
  # one value per species
  for (i in 1:Nsp){b0[i] ~ dnorm(0, tau.b0)}
  
  # one value per year
  for (t in 1:Nyr){logI[t] ~ dnorm(0, tau.I)}
  
  # one value per site-species
  for (i in 1:Nsp){   
    for (t in 1:Nyr){
      eta[i,t] ~ dnorm(0, tau.eta) # process error
      se[i,t] ~ dunif(0,1) # for the missing values
    }}
  
  # Hyperpriors
  tau.b0 ~ dt(0,1,1)T(0,)
  tau.I ~ dt(0,1,1)T(0,)
  tau.eta ~ dt(0,1,1)T(0,)
  
  # Each year-species combos is estimated with error
  for (t in 1:Nyr){
    for (i in 1:Nsp){ 
      estimate[i,t] ~ dnorm(muN[i,t], tau.psd[i,t]) # Estimated effect size
      tau.psd[i,t] <- pow(se[i,t], -2)    # 'Known' measurement error
      
      #muN is the true unknown species index this year (on the log scale)
      # it's a simple linear function of the year and site effects with process error
      muN[i,t] <- alpha + b0[i] + logI[t] + eta[i,t]
    }}
  
  # derived parameter: indicator on the measurement scale
  for (t in 1:Nyr){log(I[t]) <- logI[t]}
}"
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}
  
bma_indicator_ranwalk <- function(temp_file = tempfile()){
  
  model <- "
model {
  
  # Priors 
  alpha ~ dnorm(0, 0.0001)
  
  # one value per species
  for (i in 1:Nsp){b0[i] ~ dnorm(0, tau.b0)}
  
  # one value per year
  logI[1] ~ dnorm(0, 0.0001)
  for (t in 2:Nyr){logI[t] ~ dnorm(logI[t-1], tau.I)}
  
  # one value per site-species
  for (i in 1:Nsp){   
  for (t in 1:Nyr){
  eta[i,t] ~ dnorm(0, tau.eta) # process error
  se[i,t] ~ dunif(0,1) # for the missing values
  }}
  
  # Hyperpriors
  tau.b0 ~ dt(0,1,1)T(0,)
  tau.I ~ dt(0,1,1)T(0,)
  tau.eta ~ dt(0,1,1)T(0,)
  
  # Each year-species combos is estimated with error
  for (t in 1:Nyr){
  for (i in 1:Nsp){ 
  estimate[i,t] ~ dnorm(muN[i,t], tau.psd[i,t]) # Estimated effect size
  tau.psd[i,t] <- pow(se[i,t], -2)    # 'Known' measurement error
  
  #muN is the true unknown species index this year (on the log scale)
  # it's a simple linear function of the year and site effects, with 'process error'
  muN[i,t] <- alpha + b0[i] + logI[t] + eta[i,t]
  }}
  
  # derived parameter: indicator on the measurement scale
  for (t in 1:Nyr){log(I[t]) <- logI[t]}
  
}"
  
  writeLines(text = model, con = temp_file)
  return(temp_file)
}