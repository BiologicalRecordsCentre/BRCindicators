
##############################################################
# Generate an 'indicator'. We will try and retrieve this
# from some observed species indices
##############################################################

#' @import reshape2

# tau.spi is the precision in growth rates around the ANNUAL mean
# that modelled by sigma in this simulation
# tau.spi for wider countryside butterflies was 30.715 (95% CI: 27.7- 33.84)
# converting this to a stdv: mean = 0.180 (95% CI 0.19, 0.172)
# We are going to assume that the simulation is on the log10 scale

# sdg is the stdv on Overall growth rates among years. This is not estimated in the model

simulate_indicator <- function(nyears=30, nsp=8, mg= -0.02, sdg=0.2, sigma=0.1, SE=0.05){

#nyears<-30; nsp<-8
Tgrowth<-rnorm(nyears,mg, sdg)     # 'growth' in the indicator

#--------------------------------------------------------
# Species TRUE relative growth rates are given by average growth + random error term
#

Sgrowth<-eta<-matrix(0,nsp,nyears)
for (s in 1:nsp){
  for (t in 1:(nyears-1)){
    eta[s,t]<-rnorm(1,0,sigma)
    Sgrowth[s,t]<-Tgrowth[t] + eta[s,t]
  }}

Tindicator<-rep(0,nyears)
for ( t in 2:nyears){
  Tindicator[t]<-sum(Tgrowth[1:(t-1)])
}

Sindicator<-matrix(0,nsp,nyears)
for ( s in 1:nsp){
  Sindicator[s,1]<-rnorm(1,0,0.5)
  for ( t in 2:nyears){
    Sindicator[s,t]<-Sindicator[s,1]+sum(Sgrowth[s,1:(t-1)])
  }}

#-------------------------------------------------------------
# observed indices ofcourse have additional error

species<-matrix(0,nsp,nyears)

for (s in 1:nsp){
  for (t in 1:nyears){
    species[s,t] <- Sindicator[s,t] + rnorm(1,0,SE)
  }}

out <- melt(species)
names(out) <- c('species', 'year', 'index')
attr(out, "data") <- list(indicator=Tindicator,
                          growth = Tgrowth,
                          spgrowth = Sgrowth,
                          spabund = Sindicator)

return(out)

}
