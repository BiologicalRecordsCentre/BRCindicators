###########################################################
####  You're going to need to create a matrix 'Z' from ####
#### the series of 'years' (something like (1,2,3.....20) #
#### and teh number of 'knots' ("degrees of freedom") #####
#### an integer probably around 1/3rd of the length   #####
#### of the series. This function does that ###############
###########################################################


makeZX <- function(num.knots, covariate){
  
  n <- length(covariate)
  X <- cbind(rep(1, n), covariate)
  
  knots <- quantile(unique(covariate),
                    seq(0, 1, length = (num.knots + 2))[-c(1, (num.knots + 2))])
  
  Z_K <- (abs(outer(covariate, knots, "-")))^3
  OMEGA_all <- (abs(outer(knots, knots, "-")))^3
  
  svd.OMEGA_all <- svd(OMEGA_all)
  sqrt.OMEGA_all <- t(svd.OMEGA_all$v %*%
                      (t(svd.OMEGA_all$u) * sqrt(svd.OMEGA_all$d)))
  Z <- t(solve(sqrt.OMEGA_all, t(Z_K)))
  
  return(list(Z = Z, X = X))
  
}