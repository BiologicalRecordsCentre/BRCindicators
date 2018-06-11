# Nick Isaac, 8/6/18



# this version needs to have the data in a folder underneath the working directory
# It downloads the latest version of BRCindicators from my repository
# It reads the butterfly data from the W drive

# check you have the latest version
install_package <- FALSE
if(install_package){
  library(devtools)
  install_github(repo = 'drnickisaac/BRCindicators')
}

pdf('butterfly_results.pdf', onefile=TRUE, paper='a4')

# modified from my R Markdown document "testing_BMAs_butterflies_v3.Rmd


#  ```{r setup, include = FALSE}
library(BRCindicators)
library(reshape2)
library(RColorBrewer)
library(ggplot2)

test <- TRUE
iter <- 5e3
scaling_factor = 5 # inflate the standard errors
#nk <- 12
#```

## Introduction

#We are testing the new BMA indicators developed by Steve and Nick and implemented in BRCindicators by Tom.


## Load Butterfly data

#Use the butterfly data to compare the methods outputs. We used this butterfly data to compare methods including MSI previously so this is a good dataset to test on. 

# The 2016 butterfly data comes from here:
#dir <- 'W:/PYWELL_SHARED/Pywell Projects/BRC/Tom August/R Packages/BRCindicators/Other_files'
#load(file.path(dir, 'butterfly_data.rdata'))

# the 2017 data is stored locally
new_raw <- read.csv('data/all_spp_ci.csv')
all_species <- read.csv(file='data/all_species.csv')

head(new_raw)
head(all_species)


## Wider countryside

#Let's first run models with the wider countryside species.

#```{r wider, results = 'hide'}
## Wider
in_data <- subset(new_raw, SPECIES %in% all_species[all_species$type=='wider', 'BMScode'])

in_data <- in_data[,c('SPECIES','YEAR','TR0OBS','TR0SE')] 

names(in_data) <- c('species','year','index','se')
nrow(in_data) # 1007

# rescale for BMA
firstyr <- min(new_raw$YEAR)

length(unique(in_data$species))

# start the years at 1
in_data$year <- in_data$year - firstyr + 1
#```

#Now we can run the models.

#```{r runButterflies_wider, echo = FALSE, message = FALSE, fig.width = 6, fig.height = 4, warning = FALSE, results = 'hide'}
group <- 'wider'

in_data$se[in_data$se == 0] <- NA 
# SEs are on the natural log scale so change them to the log10 scale
#in_data$se <- in_data$se/log(10) NO - they are acutally on the correct scale
in_data$se <- in_data$se * scaling_factor

for(nk in c(4, 12)){
  
bma_indicator <- bma(in_data, model = 'smooth_det2',
                     n.iter = ifelse(test, 100, iter),
                     parallel = TRUE,
                     m.scale = 'log10',
                     num.knots = nk,
                     plot = FALSE)

plot_indicator(indicator = bma_indicator[,'Index'],
               CIs = bma_indicator[,c(3,4)], 
               main = paste0(group, ' Butterflies - smoothdet Indicator', nk,' knots'))

model <- attr(bma_indicator, 'model')
plot_growthrates(logLambda = model$mean$logLambda,
                 CIs = with(model, 
                            data.frame(lower=q2.5$logLambda, upper=q97.5$logLambda)),
                 main = paste0(group, ' Butterflies - smoothdet Growth rates', nk,' knots'))

# plot the modelled vs observed spindex
match_data <- in_data
spmatch <- data.frame(bms=unique(match_data$species), model=order(unique(match_data$species)))
match_data$species <- spmatch$model[match(match_data$species, spmatch$bms)]

plot_species(spindex = model$mean$spindex,
             lower = model$q2.5$spindex,
             upper = model$q97.5$spindex,
             data_in = match_data, main = paste(group, nk, 'knots'))

}

## Specialists

#Lets now run this with the specialist butterflies.

#```{r specialists, results = 'hide'}
## Specialists
in_data <- subset(new_raw, SPECIES %in% all_species[all_species$type=='specialist', 'BMScode'])
in_data <- in_data[,c('SPECIES','YEAR','TR0OBS','TR0SE')] 


names(in_data) <- c('species','year','index','se')
nrow(in_data) # 943

# rescale for BMA
firstyr <- min(new_raw$YEAR)

length(unique(in_data$species)) # 11 or 26, depending on whether the first year is to be included

# start the years at 1
in_data$year <- in_data$year - firstyr + 1
#```

#Now we can run the models for specialist butterflies.

#```{r runButterflies_specialists, echo = FALSE, message = FALSE, fig.width = 6, fig.height = 4, warning = FALSE, results = 'hide'}
group <- 'specialists'

in_data$se[in_data$se == 0] <- NA 
# SEs are on the natural log scale so change them to the log10 scale
#in_data$se <- in_data$se/log(10) NO - they are on the raw scale
in_data$se <- in_data$se * scaling_factor

for(nk in c(4, 12)){

bma_indicator <- bma(in_data, model = 'smooth_det2',
                     n.iter = ifelse(test, 100, iter),
                     parallel = TRUE,
                     m.scale = 'log10',
                     num.knots = nk,
                     plot = FALSE)

x<-plot_indicator(indicator = bma_indicator[,'Index'],
               CIs = bma_indicator[,c(3,4)],
               main = paste0(group, ' Butterflies - smoothdet ', nk,' knots'))

model <- attr(bma_indicator, 'model')
plot_growthrates(logLambda = model$mean$logLambda,
                 CIs = with(model, 
                            data.frame(lower=q2.5$logLambda, upper=q97.5$logLambda)),
                 main = paste0(group, ' Butterflies - smoothdet Growth rates', nk,' knots'))

# plot the modelled vs observed spindex
match_data <- in_data
spmatch <- data.frame(bms=unique(match_data$species), model=order(unique(match_data$species)))
match_data$species <- spmatch$model[match(match_data$species, spmatch$bms)]

plot_species(spindex = model$mean$spindex,
             lower = model$q2.5$spindex,
             upper = model$q97.5$spindex,
             data_in = match_data, main = paste(group, nk, 'knots'))
}


dev.off()
