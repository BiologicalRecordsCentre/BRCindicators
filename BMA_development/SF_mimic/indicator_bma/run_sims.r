

#########################################  download and install the package

library(devtools)
install_github(repo = 'drnickisaac/BRCindicators')
library(BRCindicators)

# load the help file
vignette("BRCindicators")

######################################## simulate some data

data_in <- simulate_indicator(nsp=25, nyears=40, sigma=0.18, mg= -0.01, sdg=0.005)
# NB this is the dataset I generated in January. 
# It has more species that Steve's original.
# Also, the growth rate varies more among species than Steve's
# both parameters were chosen to resemble the UKBMS

# plot the simulated data
library(ggplot2)
qplot(data=data_in, y=index, x=year) + facet_wrap(~species)

# Now it is in the right format
head(data_in)

# .. but we need an SE column, which the simulations do not contain
# This does not work with NAs, so set them all to the same alue
#data_in$se <- 0.125 # if all have the same value then stoch_smooth2 does not smooth
data_in$se <- runif(n=nrow(data_in), min=0.01, max=2)


######################################## simple growth-rate model

bma_indicator_FNgr <- bma(data_in, model = 'FNgr2', 
                          m.scale = "log10", parallel=TRUE,
                          n.iter = 5000)

plot_indicator(indicator = bma_indicator_FNgr[,'Index'],
               CIs = bma_indicator_FNgr[,c(3,4)], 
               main = 'Simulation - FNgr2')


######################################## deterministic smoother

bma_indicator_smooth_det <- bma(data_in, model = 'smooth_det2', 
                                m.scale = "log10", parallel=TRUE,
                                n.iter = 5000)

plot_indicator(indicator = bma_indicator_smooth_det[,'Index'],
               CIs = bma_indicator_smooth_det[,c(3,4)], 
               main = 'Simulation - smooth_det2')


######################################## digging deeper

# look at the R code that does the work 
bma

# look at the model object
str(bma_indicator_smooth_det)

# get the bugs model output
bugs_model <- attr(bma_indicator_smooth_det, 'model')

# find out where the bugs code is stored
bugs_model$modfile

# print out the bugs code to the screen
read.delim(bugs_model$modfile)
