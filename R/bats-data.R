#' @name bats
#' @title Data - National Bat Monitoring Programme, UK
#' @description This dataset from the National Bat Monitoring Programme, has national (UK) abundance indices for eight bat species from 1998-2014. For more information see Barlow et al (2015).
#' @docType data
#' @usage data(bats)
#' @format 
#' There are four columns of data:
#' \itemize{
#'  \item{\code{"species"}} {The species' name (an abbreviation of the common name)}
#'  \item{\code{"year"}} {year to which the index value refers (integer)}
#'  \item{\code{"collated_index"}} {A national index of abundance. The time-series for each species is fixed to have a value of 100 in the second year of data.}
#'  \item{\code{"index"}} {A national index of abundance on the natural log scale.}
#'  }
#' @references Barlow, K. E., Briggs, P. A., Haysom, K. A., Hutson, A. M., Lechiara, N. L., Racey, P. A., … Langton, S. D. (2015). Citizen science reveals trends in bat populations: The National Bat Monitoring Programme in Great Britain. Biological Conservation, 182, 14–26. \url{https://doi.org/10.1016/j.biocon.2014.11.022}
#' @source \url{https://doi.org/10.1016/j.biocon.2014.11.022}
#' @examples 
#' data(bats)
#' head(bats)
"bats"