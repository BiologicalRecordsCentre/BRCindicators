# number of species
nsp = 50

# number of years
nyr = 40

#number of iterations
iter = 500

# Build a random set of data
myArray <- array(data = rnorm(n = nsp*nyr*iter,
                              mean = 0.5,
                              sd = 0.1),
                 dim = c(nsp, nyr, iter),
                 dimnames = list(paste0('SP',1:nsp),
                                 1:nyr,
                                 1:iter))

# Ensure values are bounded by 0 and 1
myArray[myArray > 1] <- 1
myArray[myArray < 0] <- 0