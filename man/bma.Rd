% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_meta_analysis.R
\name{bma}
\alias{bma}
\title{Bayesian Meta-analysis}
\usage{
bma(
  data,
  plot = TRUE,
  model = "smooth",
  parallel = FALSE,
  n.cores = parallel::detectCores() - 1,
  incl.model = TRUE,
  n.iter = 10000,
  n.thin = 5,
  m.scale = "loge",
  num.knots = 12,
  seFromData = FALSE,
  Y1perfect = TRUE,
  rescale_indices = NULL,
  rescaleYr = 1,
  baseline = 100,
  errorY1 = FALSE,
  save.sppars = TRUE,
  incl.2deriv = FALSE,
  CI = 95,
  seed = NULL
)
}
\arguments{
\item{data}{a data.frame with 3-4 columns: `species`, `year`, `index`, `se` (standard error). The `se` column is optional 
NB: Index values are assumed to be on the unbounded (logarithmic scale)}

\item{plot}{Logical, should a trace plot be plotted to diagnose the model output?}

\item{model}{The type of model to be used. See details.}

\item{parallel}{if \code{TRUE} the model chains will be run in parallel using one fewer cores than
are available on the computer as default. NOTE: this will typically not work for parallel use on cluster PCs.}

\item{n.cores}{if running the code in parallel this option specifies the number of cores to use.}

\item{incl.model}{if \code{TRUE} the model is added as an attribute of the object returned}

\item{n.iter}{The number of iterations of the model to run. Defaults to 10,000 to avoid long run times
though much longer runs are usually required to reach convergence}

\item{n.thin}{Thinning rate for the Markov chains. Defaults to 5.}

\item{m.scale}{The measurement scale of the data. The scale of the data is assumed to be logarithmic.
Here you specify which log scale the data is on ('loge', 'log10', or 'logit'). Defaults to 'loge'.}

\item{num.knots}{If using either of the smooth models this specifies the number of knots.}

\item{seFromData}{Logical. Should the standard errors be read in from data (`TRUE`) or estimated (`FALSE`)?  Defaults to `FALSE`}

\item{Y1perfect}{Logical. Should the first year of a species' index be assumed known without error (`TRUE`)? Defaults to `TRUE`}

\item{rescale_indices}{Integer. A value for standardising each species time-series to start at a common value (e.g. 0). Defaults to NULL (i.e. no standardisation)}

\item{rescaleYr}{Integer. To which year should the indicator use as a reference value (i.e. baseline). Values greater than the number of years in the dataset will be set to the final year. Defaults to 1 (the first year)}

\item{baseline}{Integer. What is the value of the indicator in the baseline year (defaults to 100)}

\item{errorY1}{Logical. Should the indicator be presented with (`TRUE`) or without (`FALSE`) uncertainty in the baseline year. Defaults to `FALSE`.}

\item{save.sppars}{Logical. Should the species-specific parameters be monitored? Defaults to TRUE}

\item{incl.2deriv}{Logical. Option to include estimation of second derivatives on the indicator (`TRUE`)? Defaults to `FALSE`}

\item{CI}{defines the credible intervals of the posterior distribution to report. Defaults the 95th percentile}

\item{seed}{Option to set a custom seed to initialize JAGS chains, for reproducibility. Should be an integer. This argument will be deprecated in the next version, but you can always set the outside the function yourself.}
}
\value{
Returns a dataframe with 7 columns: Year, Index.Mprime, lowerCI.Mprime, upperCI.Mprime, Index.M, lowerCI.M and, upperCI.M. 
Columns headed `M` and `Mprime` are means of the M and M' parameters  as defined in Freeman et al (2020). The 'upper' and 'lower' columns are the credible intervals, the width of which is defined by the `CI` argument.
Note that M and M' are alternate ways of calculating the multispecies indicator: their means are nearly always virtually identical, but the uncertainty in M is usually much wider than in M'. See Freeman et al (2020) for more details.
}
\description{
Use a Bayesian meta-analysis to create an indicator from species index values, optionally incorporating standard error.
}
\details{
There are a number of model to choose from:
\itemize{
 \item{\code{"smooth"}}{ The default option. Indicator defined by Growth rates, with Ruppert smoother, allowing for species to join late. Error on the first year of each species' time-series is assumed to be zero. The indicator is the expected value of the geometric mean across species (with missing species imputed). 
 Includes three options: `seFromData` `Y1perfect` and `incl.2deriv`. See bayesian_meta_analysis for mode details. Using the default values `seFromData = FALSE` and `Y1perfect = TRUE` are the options used in Freeman  \emph{et al.} (2020).}
 \item{\code{"smooth_det2"}}{ Equivalent to smooth with `seFromData = TRUE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
 \item{\code{"smooth_det_sigtheta"}}{ Equivalent to smooth with `seFromData = FALSE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
 }
}
\examples{

# Only run if there is a JAGS installation
if(suppressWarnings(runjags::testjags(silent = TRUE)$JAGS.found)){

# Create some example data in the format required
data <- data.frame(species = rep(letters, each = 50),
                   year = rep(1:50, length(letters)),
                   index = rnorm(n = 50 * length(letters), mean = 0, sd = 1),
                   se = runif(n = 50 * length(letters), min = 0.01, max = .1))

# Run the Bayesian meta-analysis
bma_indicator <- bma(data, model="smooth", m.scale="logit", n.iter=100)

# Plot the resulting indicator
plot_indicator(indicator = bma_indicator[,'Index.Mprime'],
               CIs = bma_indicator[,c(3,4)])
   
   }
}
\references{
Freeman, S.N., Isaac, N.J.B., Besbeas, P.T., Dennis, E.B. & Morgan, B.J.T. (2020) 
            A generic method for estimating and smoothing multispecies biodiversity indices, robust to intermittent data. 
            \emph{Journal of Agricultural Biological and Environmental Statistics}, in revision.
}
