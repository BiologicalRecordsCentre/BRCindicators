# BRCindicators
<!-- badges: start -->
[![R-CMD-check](https://github.com/BiologicalRecordsCentre/BRCindicators/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BiologicalRecordsCentre/BRCindicators/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/biologicalrecordscentre/BRCindicators/branch/master/graph/badge.svg)](https://codecov.io/gh/biologicalrecordscentre/BRCindicators?branch=master)
<!-- badges: end -->


The functions in BRCindicators work with yearly estimates of species abundance or occurrence and aggregate them into an scaled indicator value with bootstrapped confidence intervals 

Installing the package is easy and can be done in a couple of lines in R

    library(devtools)
    install_github(repo = 'biologicalrecordscentre/BRCindicators')

For more info read the vignette

    vignette('BRCindicators')
