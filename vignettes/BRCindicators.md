---
title: "BRCindicators"
author: "Tom August"
date: "04 March, 2019"
output: 
  html_document:
    keep_md: yes
    toc: yes
vignette: >
  %\VignetteIndexEntry{BRCindicators}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---



```r
require(snowfall)
```

```
## Loading required package: snowfall
```

```
## Loading required package: snow
```

# Introduction

This document shows how to use the indicator pipeline to create biodiversity indicators such as those for DEFRA's [Biodiversity Indicators in Your Pocket](http://jncc.defra.gov.uk/page-1824). The pipeline is shared in the form of an R package called 'BRCindicators' making it easy to share and maintain.

The functions in BRCindicators work with yearly estimates of species abundance or occurrence and aggregate them into an scaled indicator value with bootstrapped confidence intervals 

This package has the ability to read in the output of occupancy models created in the R package sparta, a package for estimating species trends from occurrence data. This package can be installed [from Github](https://github.com/BiologicalRecordsCentre/sparta) and details of how to use the package are given in [the package vignette](https://github.com/BiologicalRecordsCentre/sparta/raw/master/vignette/sparta_vignette.pdf). There is no need to use sparta to create your yearly species estimates as BRCindicators can also work with other data. 

To create an indicator we first need to have species trends, let's create some using the sparta R package.

# Creating yearly estimates of occurrence in sparta

If you already have yearly estimates of abundance or occurrence for your species you can skip this stage. Here we show how you can create these estimates from raw species observation data using sparta.


```r
# Install the package from CRAN
# THIS WILL WORK ONLY AFTER THE PACKAGE IS PUBLISHED
install.packages('sparta')

# Or install the development version from GitHub
library(devtools)
install_github('biologicalrecordscentre/sparta')
```

Let's assume you have some raw data already, we can under take occupancy modelling like this


```r
# Create data
n <- 8000 # size of dataset
nyr <- 50 # number of years in data
nSamples <- 200 # set number of dates
nSites <- 100 # set number of sites
set.seed(125) # set a random seed

# Create somes dates
first <- as.Date(strptime("1950/01/01", "%Y/%m/%d")) 
last <- as.Date(strptime(paste(1950+(nyr-1),"/12/31", sep=''), "%Y/%m/%d")) 
dt <- last-first 
rDates <- first + (runif(nSamples)*dt)

# taxa are set semi-randomly
taxa_probabilities <- seq(from = 0.1, to = 0.7, length.out = 26)
taxa <- sample(letters, size = n, TRUE, prob = taxa_probabilities)

# sites are visited semi-randomly
site_probabilities <- seq(from = 0.1, to = 0.7, length.out = nSites)
site <- sample(paste('A', 1:nSites, sep=''), size = n, TRUE, prob = site_probabilities)

# the date of visit is selected semi-randomly from those created earlier
time_probabilities <- seq(from = 0.1, to = 0.7, length.out = nSamples)
time_period <- sample(rDates, size = n, TRUE, prob = time_probabilities)

myData <- data.frame(taxa, site, time_period)
```


```r
# Load the sparta package
library(sparta)
```

```
## Loading required package: lme4
```

```
## Loading required package: Matrix
```

For demonstration purposes I have a faked dataset of 8000 species observations. In my dataset the species are named after the letters in the alphabet. Below I show how I can use the Bayesian occupancy models in sparta to create yearly estimates of occurrence. For more information please see the [vignette for sparta](https://github.com/BiologicalRecordsCentre/sparta/raw/master/vignette/sparta_vignette.pdf)


```r
# Preview of my data
head(myData)
```

```
##   taxa site time_period
## 1    r  A51  1970-01-14
## 2    v  A87  1980-09-29
## 3    e  A56  1996-04-14
## 4    z  A28  1959-01-16
## 5    r  A77  1970-09-21
## 6    x  A48  1990-02-25
```

```r
# First format our data
formattedOccData <- formatOccData(taxa = myData$taxa,
                                  site = myData$site,
                                  survey = myData$time_period)
```

```
## Warning in errorChecks(taxa = taxa, site = site, survey = survey,
## closure_period = closure_period): 94 out of 8000 observations will be
## removed as duplicates
```

```r
# Here we are going to use the package snowfall to parallelise
library(snowfall)

# I have 4 cpus on my PC so I set cpus to 4
# when I initialise the cluster
sfInit(parallel = TRUE, cpus = 4)
```

```
## Warning in searchCommandline(parallel, cpus = cpus, type = type,
## socketHosts = socketHosts, : Unknown option on commandline:
## rmarkdown::render('W:/PYWELL_SHARED/Pywell Projects/BRC/Tom August/R
## Packages/BRCindicators/vignettes/BRCindicators.Rmd', encoding
```

```
## R Version:  R version 3.5.2 (2018-12-20)
```

```
## snowfall 1.84-6.1 initialized (using snow 0.4-3): parallel execution on 4 CPUs.
```

```r
# Export my data to the cluster
sfExport('formattedOccData')

# I create a function that takes a species name and runs my model
occ_mod_function <- function(taxa_name){
  
  library(sparta)
  
  # Note that this will write you results to your computer
  # the location is set to your user folder
  occ_out <- occDetFunc(taxa_name = as.character(taxa_name),
                        n_iterations = 200,
                        burnin = 15, 
                        occDetdata = formattedOccData$occDetdata,
                        spp_vis = formattedOccData$spp_vis,
                        write_results = TRUE,
                        output_dir = '~/Testing_indicator_pipe',
                        seed = 123)  
} 

# I then run this in parallel
system.time({
para_out <- sfClusterApplyLB(unique(myData$taxa), occ_mod_function)
})
```

```
##    user  system elapsed 
##    0.04    0.05  139.66
```

```r
# Stop the cluster
sfStop()
```

```
## 
## Stopping cluster
```

```r
# We can see all the files this has created
list.files('~/Testing_indicator_pipe')
```

```
##  [1] "a.rdata" "b.rdata" "c.rdata" "d.rdata" "e.rdata" "f.rdata" "g.rdata"
##  [8] "h.rdata" "i.rdata" "j.rdata" "k.rdata" "l.rdata" "m.rdata" "n.rdata"
## [15] "o.rdata" "p.rdata" "q.rdata" "r.rdata" "s.rdata" "t.rdata" "u.rdata"
## [22] "v.rdata" "w.rdata" "x.rdata" "y.rdata" "z.rdata"
```

# Installing BRCindicators

Installing the package is easy and can be done in a couple of lines


```r
library(devtools)
install_github(repo = 'biologicalrecordscentre/BRCindicators')
```

# Summarising sparta output for an indicator

Now that we have some species trends data to work with (no doubt you already have your own) we can use the first function in BRCindicators. This function reads in all the output files from sparta (which are quite large and complex) and returns a simple summary table that we can use for calculating the indicator. If you have done your analysis without using sparta you can skip to the next step.


```r
library(BRCindicators)

# All we have to supply is the directory where out data is saved
# You will note this is the 'output_dir' passed to sparta above.
trends_summary <- summarise_occDet(input_dir = '~/Testing_indicator_pipe')
```

```
## Loading data...done
```

```r
# Lets see the summary
head(trends_summary[,1:5])
```

```
##      year         a         b         c         d
## [1,] 1950 0.6936066 0.6412022 0.5689617 0.5595082
## [2,] 1951 0.7638251 0.4861749 0.4021311 0.5510929
## [3,] 1952 0.6248634 0.6434973 0.2830601 0.5293443
## [4,] 1953 0.7024590 0.7724590 0.4861202 0.3753552
## [5,] 1954 0.6174317 0.5121858 0.7556831 0.3475956
## [6,] 1955 0.4302732 0.3584699 0.6068852 0.5328962
```

Returned from this function is a summary of the data as a matrix. In each row we have the year, specified in the first column, and each subsequent column is a species. The values in the table are the mean of the posterior for the predicted proportion of sites occupied, a measure of occurrence.

# Calculating indicator values

Once we have species-year indicies we are in a position to proceed to calculating an indictor. To do this there are a number of mehods available, some of which are presented here in 'BRCindicators'

## Geometric mean

The geometric mean method is often used with data that do not have errors associated with them.

The first step is to re-scale the data so that the value for all species in the first year is the same. Once this is done we calculate the geometric mean across species for each year creating the indicator value. This function also accounts for species that have no data at the beginning of the dataset by entering them at the geometric mean for that year, this stops them dramatically changing the indicator value in the year they join the dataset. It also accounts for species that leave the dataset before the end by holding them at their last value. Finally limits to species values can be given, preventing extremely high or low values biasing the indicator.

### Rescaling and calculating geometric mean

The data I have generated in 'trends_summary' is very easy to work with but to show off what this function can do I'm going to mess it up a bit.


```r
trends_summary[1:3, 'a'] <- NA
trends_summary[1:5, 'b'] <- NA
trends_summary[2:4, 'c'] <- 1000
trends_summary[45:50, 'd'] <- NA

# Let's have a look at these changes
head(trends_summary[,1:5])
```

```
##      year         a         b            c         d
## [1,] 1950        NA        NA    0.5689617 0.5595082
## [2,] 1951        NA        NA 1000.0000000 0.5510929
## [3,] 1952        NA        NA 1000.0000000 0.5293443
## [4,] 1953 0.7024590        NA 1000.0000000 0.3753552
## [5,] 1954 0.6174317        NA    0.7556831 0.3475956
## [6,] 1955 0.4302732 0.3584699    0.6068852 0.5328962
```

```r
tail(trends_summary[,1:5])
```

```
##       year         a         b         c  d
## [45,] 1994 0.1308743 0.6617486 0.3050273 NA
## [46,] 1995 0.5324044 0.2455191 0.7130601 NA
## [47,] 1996 0.7407650 0.4489617 0.5519672 NA
## [48,] 1997 0.3051366 0.7181421 0.7871585 NA
## [49,] 1998 0.7282514 0.2540984 0.5586339 NA
## [50,] 1999 0.5635519 0.4161202 0.4713661 NA
```

Now that I have 'messed up' the data a bit we have two species with data missing at the beginning and one species with data missing at the end. We also have one species with some very high values.

Now lets run this through the re-scaling function.


```r
# Let's run this data through our scaling function (all defaults used)
rescaled_trends <- rescale_species(Data = trends_summary)

# Here's the result
head(rescaled_trends[,c('year', 'indicator', 'a', 'b', 'c', 'd')])
```

```
##      year indicator         a        b          c         d
## [1,] 1950 100.00000        NA       NA   100.0000 100.00000
## [2,] 1951 119.03406        NA       NA 10000.0000  98.49595
## [3,] 1952 126.49450        NA       NA 10000.0000  94.60885
## [4,] 1953 118.53382 118.53382       NA 10000.0000  67.08663
## [5,] 1954  94.25122 104.18620       NA   132.8179  62.12521
## [6,] 1955  99.15398  72.60485 99.15398   106.6654  95.24368
```

```r
tail(rescaled_trends[,c('year', 'indicator', 'a', 'b', 'c', 'd')])
```

```
##       year indicator         a         b         c        d
## [45,] 1994  92.77695  22.08390 183.04187  53.61122 140.7169
## [46,] 1995 102.48140  89.83859  67.91140 125.32655 140.7169
## [47,] 1996 105.83425 124.99763 124.18431  97.01306 140.7169
## [48,] 1997 106.16325  51.48914 198.64048 138.34998 140.7169
## [49,] 1998  94.26861 122.88605  70.28445  98.18479 140.7169
## [50,] 1999  95.96255  95.09446 115.10023  82.84672 140.7169
```

You can see that species 'a' and 'b' enter the dataset at the geometric mean (the indicator value), all species are indexed at 100 in the first year and the very high values in 'c' are capped at 10000 at the end 'd' has been held at it's end value.

The 'indicator' column that is returned here is our indicator, calculated as the geometric mean of all the species in the data set. 

### Confidence intervals

We can get confidence intervals for this indicator by bootstrapping across species. We have a function for that too!


```r
# This function takes just the species columns
scaled_species <- rescaled_trends[,!colnames(rescaled_trends) %in% c('year', 'indicator')]
indicator_CIs <- bootstrap_indicator(Data = scaled_species)
```

```
## Running bootstrapping for 10000 iterations...done
```

```r
# Returned are the CIs for our indicator
head(indicator_CIs)
```

```
##      quant_025 quant_975
## [1,] 100.00000  100.0000
## [2,]  87.33579  190.2349
## [3,]  95.60814  198.0757
## [4,]  88.79716  183.7160
## [5,]  80.60150  110.3239
## [6,]  86.15162  115.1020
```

### Smoothing

It is sometimes desirable to create a smoothed indicator value from the raw values. This can be achieved by fitting a GAM (general additive model) to the indicator using a spline. This spline is a smoothed curve that goes through the raw values for the indicator and is fitted using the function 'gam' in the 'mgcv' R package.


```r
# The smoothing function takes the indicator values
smoothed_indicator <- GAM_smoothing(rescaled_trends[,'indicator'])

# In this example there is little support for a non-linear trend and 
# so the line almost linear
plot(x = rescaled_trends[,'year'], y = rescaled_trends[,'indicator'])
lines(x = rescaled_trends[,'year'], y = smoothed_indicator, col = 'red')
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/smoothing_indicator-1.png)<!-- -->

```r
# But if our indicator did support a non-linear trend it might look 
# like this
eg_indicator <- jitter(sort(rnorm(50)), amount = 0.5)
eg_smoothed <- GAM_smoothing(eg_indicator)
plot(x = 1:50, y = eg_indicator)
lines(x = 1:50, y = eg_smoothed, col = 'red')
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/smoothing_indicator-2.png)<!-- -->

Where there is little support for a non-linear trend a GAM smoothed line will tend towards linear. Where there is good support for a non-linear trend the smoothed line will become more 'bendy'.

### Plotting

We now have our indicator and the confidence intervals around it. The next step is to plot it. We have included a function that creates a simple plot using ggplot2, however you could easily create your own plots in R using the data.


```r
# Plot our indicator.
plot_indicator(indicator = rescaled_trends[,'indicator'],
               smoothed_line = smoothed_indicator,
               CIs = indicator_CIs)
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/plot_indicator-1.png)<!-- -->

In this plot you can see the high upper confidence interval in years 2-4, this is due to the artificially high values we gave to species 'c'.

## Bayesian Meta-Analysis (BMA)

The Bayesian Meta-Analysis method, or BMA, is suited to data with standard errors associated with them. As with other methods we require data from more than one species, across a number of years, with an error for each species-year estimate.


```r
# Here is an example dataset for the BMA method
data <- data.frame(species = rep(letters, each = 50),
                   year = rep(1:50, length(letters)), 
                   index = runif(n = 50 * length(letters), min = 0, max = 1), 
                   se = runif(n = 50 * length(letters), min = 0.01, max = .1))
head(data)
```

```
##   species year     index         se
## 1       a    1 0.5377357 0.01769631
## 2       a    2 0.9644076 0.07268163
## 3       a    3 0.3653818 0.08335768
## 4       a    4 0.7266121 0.04494251
## 5       a    5 0.8660005 0.09993223
## 6       a    6 0.4728771 0.09887020
```

It is important that your data is in the same format and that your columns are in the same order and have the same names. Remember you can use the function `read.csv()` to read in the data from a .csv on your computer.

BMA is run using the function `bma`, here we will use the default settings and then see what we can change.


```r
bma_indicator <- bma(data)
```

```
## 
## Processing function input....... 
## 
## Done. 
##  
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 2600
##    Unobserved stochastic nodes: 1379
##    Total graph size: 6587
## 
## Initializing model
## 
## Adaptive phase..... 
## Adaptive phase complete 
##  
## 
##  Burn-in phase, 5000 iterations x 3 chains 
##  
## 
## Sampling from joint posterior, 5000 iterations x 3 chains 
##  
## 
## Calculating statistics.......
```

```
## Warning in doTryCatch(return(expr), name, parentenv, handler): At least one
## Rhat value could not be calculated.
```

```
## 
## Done.
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/runBMA-1.png)<!-- -->

The function returns a plot to your screen which is a diagnostic plot of the model. When the model has converged (i.e. reached a point where the three chains agree on the answer) the lines on the plots on the left will sit on top of one another and the plots on the right will have a nice bell shape. You can turn off this plot by setting `plot` to `FALSE`. By default the method runs the chains in series. Running them in parallel makes the models run faster (about half the time) but will slow down your computer more. We can change this with the parameter `parallel`. The number of iterations the model runs is controlled by `n.iter` and defaults to 10000. If you can it is better to run it for more iterations, though this will take longer. `m.scale` gives the scale your data is on. It is very important that this is correct, choose from 'loge' (natural log, sometimes simply called 'log'), 'log10' (log to the base 10), or 'logit' (output from models of proportions or probabilities).

Let's implement a few of these changes


```r
bma_indicator2 <- bma(data,
                     parallel = TRUE,
                     n.iter = 500,
                     m.scale = 'log10')
```

```
## 
## Processing function input....... 
## 
## Done. 
##  
## Beginning parallel processing using 3 cores. Console output will be suppressed.
## 
## Parallel processing completed.
## 
## Calculating statistics.......
```

```
## Warning in doTryCatch(return(expr), name, parentenv, handler): At least one
## Rhat value could not be calculated.
```

```
## 
## Done.
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/runBMAparameters-1.png)<!-- -->

Because we have reduced the number of interations the model no longer has a good convergence. The lines on the graphs on the left do not overlap and the graphs on the right are no longer a smooth bell shape.

The object that is returned is a data.frame with years as rows and columns giving the year value, index value and confidence intervals. You can write this to a csv using the function `write.csv`.


```r
head(bma_indicator)
```

```
##   Year    Index  lower2.5 upper97.5
## 1    1 100.0000 100.00000  100.0000
## 2    2 100.1559  98.00631  102.3808
## 3    3 100.3311  97.39855  103.5038
## 4    4 100.5416  96.95652  104.3762
## 5    5 100.5611  96.63643  104.8945
## 6    6 100.4664  95.99929  105.0179
```

We can use the plotting function in BRCindicators to plot the results of this analysis, which in this case are not all that interesting!


```r
plot_indicator(indicator = bma_indicator[,'Index'],
               CIs = bma_indicator[,c(3,4)])
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/BMAplot-1.png)<!-- -->


## Multi-species Indicator

The multi-species indicator method was developed by Statistics Netherlands and the code is made available on [their website](https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim--/msi-tool). To find out more about the inner working of this method please read the [detailed documentation](https://www.cbs.nl/-/media/_pdf/2017/22/msi_manual.pdf) on the authors website. Here is a simple example of how this method runs in `BRCindicators`.


```r
# Create some example data in the format required
nyr = 20
species = rep(letters, each = nyr)
year = rev(rep(1:nyr, length(letters)))

# Create an index value that increases with time
index = rep(seq(50, 100, length.out = nyr), length(letters))
# Add randomness to species
index = index * runif(n = length(index), 0.7, 1.3)
# Add correlated randomness across species, to years
index = index * rep(runif(0.8, 1.2, n = nyr), length(letters))

se = runif(n = nyr * length(letters), min = 10, max = 20)

data <- data.frame(species, year, index, se)

# Our species are decreasing
plot(data$year, data$index)
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi1-1.png)<!-- -->

```r
# Species index values need to be 100 in the base year. Here I use
# the first year as my base year and rescale to 100. The standard error
# in the base year should be 0.
min_year <- min(data$year)

for(sp in unique(data$species)){

  subset_data <- data[data$species == sp, ]
  multi_factor <- 100 / subset_data$index[subset_data$year == min_year]
  data$index[data$species == sp] <- data$index[data$species == sp] * multi_factor
  data$se[data$species == sp] <- data$se[data$species == sp] * multi_factor
  data$se[data$species == sp][1] <- 0

}

# Our first year is now indexed at 100
plot(data$year, data$index)
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi1-2.png)<!-- -->

```r
# Alternativly I could read in data from a csv
# data <- read.csv('path/to/my/data.csv')

# Run the MSI function
msi_out <- msi(data)
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi1-3.png)<!-- -->![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi1-4.png)<!-- -->

```r
head(msi_out$CV)
```

```
##   species   mean_CV
## 1       a 0.1927516
## 2       b 0.1918315
## 3       c 0.2048356
## 4       d 0.2195197
## 5       e 0.2048833
## 6       f 0.2096112
```

```r
# I can capture the output figures too
# pdf('test.pdf')
#   msi_out <- msi(data)
# dev.off()
```

The code returns two plots to the console, the first plot shows the coefficient of variation (CV) for each of the species. Species with high values of CV may adversly effect the relaibility of the trend estimation. Use this graph to identify the CV values of the species and use the `maxCV` parameter to set a threshold above which species will be excluded. The results of excluding species in this way can be tested by comparing trend plots. The CV values are hard to assign to species from this plot as the species are coded to numbers. To see the raw values look at the CV component of msi_out (i.e. `msi_out$CV`). The second plot shows the smoothed trend and the MSI values. These two figures can be captured in the usual way in R by using `pdf()` for example. In the example I create a dataset from random numbers but usually you would use `read.csv()` to read in data from a local file.

Here is a second example which sets some additional parameters. The parameters for `msi` get passed to `msi_tool` so to see a list of all the parameters you can change look at the help documentation in `msi_tool` usign `?msi_tool` at the R console. I cover most of hte important ones here.


```r
msi_out <- msi(data,
               nsim = 500, # The number of Mote Carlo simulations
               SEbaseyear = 10, # The year to index on
               plotbaseyear = 15, # The year to set as 100 in plots
               index_smoot = 'INDEX', # plotbaseyear uses MSI not trend
               span = 0.7, # 'wigglyness' of line, between 0 and 1
               lastyears = 5, # last X years of time series for short-term trends
               maxCV = 10, # maximum allowed Coefficient of Variation 
               changepoint = 10, # compare trends before and after this year
               truncfac = 8, # max year-to-year index ratio
               TRUNC = 5, #set all indices below TRUNC to this
               plot = TRUE # should the plots be returned?)
               )
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi2-1.png)<!-- -->![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi2-2.png)<!-- -->

This set of parameters is unrealistic but shows the options available. Note that in the second graph the year 10 point now has a se = 0, year 15 MSI is set to 100, and the short term trend is reported for the last 5 years.

The analysis also returns data which provide more insights into the analysis and let you create your own plots if required.


```r
# The returned object has 2 elements
head(msi_out$results)
```

```
##   year    MSI sd_MSI lower_CL_MSI upper_CL_MSI  Trend lower_CL_trend
## 1    1 164.57   7.50       150.50       179.96 162.23         151.14
## 2    2 164.60   7.11       151.24       179.13 156.51         149.63
## 3    3 124.16   6.20       112.59       136.94 151.27         145.21
## 4    4 154.53   6.93       141.53       168.72 146.51         140.92
## 5    5 170.74   7.36       156.90       185.78 142.41         136.76
## 6    6 143.41   6.76       130.77       157.30 139.41         134.05
##   upper_CL_trend      trend_class
## 1         173.85 moderate_decline
## 2         163.73 moderate_decline
## 3         157.31 moderate_decline
## 4         152.66 moderate_decline
## 5         148.15 moderate_decline
## 6         145.21 moderate_decline
```

The first of the two elements (`results`) returned gives all the data, and a little more, that is presented in the second figure.


```r
# The returned object has 2 elements
msi_out$trends
```

```
##                             Measure    value     significance
## 1                     overall trend   0.9691 moderate decline
## 2                  SE overall trend   0.0019                 
## 3                trend last 5 years   0.9711        uncertain
## 4             SE trend last 5 years   0.0154                 
## 5                  changepoint (10)  10.0000                 
## 6     trend before changepoint (10)   0.9789 moderate decline
## 7  SE trend before changepoint (10)   0.0046                 
## 8      trend after changepoint (10)   0.9487 moderate decline
## 9   SE trend after changepoint (10)   0.0040                 
## 10                         % change -46.1560           p<0.01
## 11                      SE % change   2.5360                 
## 12            % change last 5 years  -9.4120           p<0.01
## 13         SE % change last 5 years   2.3900                 
## 14                      changepoint       NA           p<0.01
```

```r
# I could write this as a csv too
# write.csv(msi_out$trends, file = 'path/to/my/output.csv')
```

The second element (`trends`) returned give a summary of various trend assessments across the time series.

We have also added a plot method for the MSI output which provides a plot similar to that of the second figure we have seen already. Lets use this plot method to explore the effect of changing the span value in the analysis


```r
for(i in c(0.3, 0.5, 0.7)){ # use a range of values for span
  
  msi_out <- msi(data, span = i, # span is set to i
                 nsim = 200, plot = FALSE)
  print( # print makes the plot visible in the for loop
    plot(msi_out, title = paste('MSI - span =', i)) # plot
  )
  
}
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi_span-1.png)<!-- -->![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi_span-2.png)<!-- -->![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/msi_span-3.png)<!-- -->

As the value of span gets closer to 1 the trend line gets smoother.

## Lambda Indicator

The lambda indicator calculates an indicator using growth rates from one year to the next. Formulating the indicator in terms of growth rates has two distinct advantages over the conventional approach to constructing indicators.  First, it means that the categorisation of species as ‘increasing’ or ‘decreasing’ can be made from the same set of data (the growth rates) as the construction of the indicator.  Second, it provides an elegant solution to the problem of species that join the indicator after the first year (i.e. where the first year is unreliable): other indicators typically adopt a complicated rescaling approach to ensure that species entering the indicator after the first year do not bias the overall assessment.  It also makes a simple and robust, though untestable, assumption about species that drop out of the indicator prior to the final year: specifically it assumes that their fluctuations are the same, in aggregate, as those of the species that remain in the indicator. For more details see http://webarchive.nationalarchives.gov.uk/20170302170037/http://jncc.defra.gov.uk/Docs/UKBI2015_TechBG_C4b-D1c_Bayesian_Final.docx 

Very few species’ models produced reliable occupancy estimates for every year, so a majority of the time series contain missing values.  This presents a problem for estimating growth rates for each species-year combination.  Missing values of growth rate that would be equivalent to linear interpolation of the log odds between adjacent years with reliable estimates were therefore calculated. This indicator can therefore work with missing values.

Input data is on the occupancy scale, and is therefore bounded between 0 and 1.


```r
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

str(myArray)
```

```
##  num [1:50, 1:40, 1:500] 0.555 0.578 0.489 0.562 0.55 ...
##  - attr(*, "dimnames")=List of 3
##   ..$ : chr [1:50] "SP1" "SP2" "SP3" "SP4" ...
##   ..$ : chr [1:40] "1" "2" "3" "4" ...
##   ..$ : chr [1:500] "1" "2" "3" "4" ...
```

`lambda_indicator` takes in an array of data, a three-dimensional matrix. The dimensions of this array represent species, years, and iterations.



```r
# Run the lambda_interpolation method on this data
myIndicator <- lambda_indicator(myArray)

# Plot the indicator
plot_indicator(myIndicator$summary[,'indicator'],
               myIndicator$summary[,c('lower' ,'upper')])
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/lambda_2-1.png)<!-- -->

There are a number of options available in the `lambda_indicator` function


```r
myIndicator <- lambda_indicator(myArray,
                                index = 1, # Set the index value to 1 not 100
                                year_range = c(30,40), # Set year range
                                threshold_yrs = 5) # set a threshold
plot_indicator(myIndicator$summary[,'indicator'],
               myIndicator$summary[,c('lower' ,'upper')])
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/lambda_3-1.png)<!-- -->

Note that there are a range of threshold functions that allow you to adjust which data points are used in the indicator. There are options to remove species year estimates based on their standard deviaction, rHat value and based on the number of years a species is present in the dataset. Note that the rhat threshold can only be used if you are using a directory path as you input rather than an array.

# Creating a custom pipeline function

We have demonstrated how you might run the indicator functions one at a time, however in a 'pipeline' we want data to flow through seamlessly. Additionally there are a number of parameters in the functions that we have not shown you that you might find useful. Here is an example of how you can create your own pipeline function. Our function will wrap around the functions described above, setting the parameters to meet our needs. Once we have done this it will allow use to execute our pipeline in one line.


```r
# I call my function 'run_pipeline' and the only arguement it
# takes is the directory of sparta's output
run_pipeline <- function(input_dir){

  require(sparta)
  require(BRCindicators)
  
  # Create the trends summary
  trends_summary <- summarise_occDet(input_dir = input_dir)

  # Rescale the values and get the indicator values
  # Here I set the index to 1 and change the value limits
  rescaled_trends <- rescale_species(Data = trends_summary,
                                     index = 1,
                                     max = 100,
                                     min = 0.001)
  
  # Bootstrap the indicator to get CIs
  scaled_species <- rescaled_trends[,!colnames(rescaled_trends) %in% c('year', 'indicator')]
  # This time I set the iterations to twice the default and 
  # use custom confidence intervals
  indicator_CIs <- bootstrap_indicator(Data = scaled_species,
                                       CI_limits = c(0.25, 0.75),
                                       iterations = 20000)
  
  # Get the smoothed indicator line
  smoothed_indicator <- GAM_smoothing(rescaled_trends[,'indicator'])
  
  # This time I specify the years and index value
  plot_indicator(indicator = rescaled_trends[,'indicator'],
                 year = rescaled_trends[,'year'],
                 index = 1,
                 CIs = indicator_CIs,
                 smoothed_line = smoothed_indicator)
  
  ## I'll return all my data  
  return(cbind(smoothed_indicator, indicator_CIs, as.data.frame(trends_summary)))
 }
```

Once we have created this function we can run this pipeline on a directory in one line, or put it in a loop to run across many directories. 


```r
# Now we can run the pipeline in one line, like a boss
indicator_data <- run_pipeline(input_dir = '~/Testing_indicator_pipe')
```

```
## Loading data...done
## Running bootstrapping for 20000 iterations...done
```

![](C:/Users/tomaug/AppData/Local/Temp/RtmpwDxzMK/preview-28d07a093333.dir/BRCindicators_files/figure-html/in_one-1.png)<!-- -->

```r
head(indicator_data)
```

```
##   smoothed_indicator  quant_25 quant_75 year         a         b         c
## 1          0.9410110 1.0000000 1.000000 1950 0.6936066 0.6412022 0.5689617
## 2          0.9410241 0.9161441 1.013355 1951 0.7638251 0.4861749 0.4021311
## 3          0.9410373 0.9654977 1.054758 1952 0.6248634 0.6434973 0.2830601
## 4          0.9410504 0.9418938 1.025378 1953 0.7024590 0.7724590 0.4861202
## 5          0.9410636 0.8829743 0.980575 1954 0.6174317 0.5121858 0.7556831
## 6          0.9410767 0.9152991 1.016191 1955 0.4302732 0.3584699 0.6068852
##           d         e         f         g         h         i         j
## 1 0.5595082 0.6493443 0.6055738 0.6360656 0.6808197 0.5023497 0.7240984
## 2 0.5510929 0.6752459 0.7466667 0.6198361 0.3294536 0.7777049 0.5445902
## 3 0.5293443 0.7142623 0.4499454 0.8585246 0.6218033 0.7155738 0.4503825
## 4 0.3753552 0.4351366 0.7007650 0.6812568 0.2843716 0.7617486 0.7204918
## 5 0.3475956 0.7890710 0.9208743 0.4224590 0.3489617 0.6709290 0.6630055
## 6 0.5328962 0.6496175 0.5071038 0.8983607 0.8183607 0.1828415 0.5910383
##           k         l         m         n         o         p         q
## 1 0.6745355 0.7474317 0.4280874 0.8637705 0.4754098 0.3063934 0.8561749
## 2 0.5612568 0.8439891 0.7219672 0.7681421 0.8899454 0.8829508 0.6603279
## 3 0.7230055 0.6057377 0.8093443 0.6449727 0.8291803 0.5590164 0.9317486
## 4 0.6781421 0.3494536 0.6693443 0.8687432 0.7666667 0.3861202 0.8895628
## 5 0.4267760 0.7104918 0.5197268 0.6080328 0.9333880 0.6397814 0.8000000
## 6 0.4730601 0.7218579 0.7813115 0.8701093 0.7691803 0.9534973 0.7195628
##           r         s         t         u         v         w         x
## 1 0.4616940 0.8857377 0.9212568 0.9278142 0.7791803 0.8990164 0.8789617
## 2 0.4956284 0.5505464 0.6045902 0.9059016 0.5865027 0.9578689 0.5522404
## 3 0.8715847 0.8665574 0.7143169 0.9012022 0.7806557 0.8487978 0.6573224
## 4 0.7289617 0.9055191 0.9448634 0.7619126 0.8996721 0.7424044 0.8715847
## 5 0.7063934 0.4081967 0.8507650 0.6043716 0.8098361 0.8346448 0.4581421
## 6 0.5015847 0.6086339 0.7546448 0.9357923 0.9027322 0.8691257 0.8631148
##           y         z
## 1 0.8712022 0.8991257
## 2 0.8714208 0.6720765
## 3 0.8427869 0.7418033
## 4 0.8560656 0.7687432
## 5 0.9773770 0.6566667
## 6 0.8149180 0.6977049
```
