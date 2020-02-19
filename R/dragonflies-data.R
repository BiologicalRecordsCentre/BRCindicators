#' @name dragonflies
#' @title Data - British Dragonfly Society, UK
#' @description This dataset contains annual occupancy estimates for 45 species of dragonflies and damselflies in the UK, 1970-2015.
#' There are 2039 rows: 44 species have an estimate for every year, but the Small red-eyed damselfly has no estimates prior to 2001 (the first year in which it was recorded in UK).
#' The models are based on biological records data from the British Dragonfly Society.
#' The data are derived from Bayesian occupancy detection models, described in Outhwaite et al (2019).
#' @usage data(dragonflies)
#' @format 
#' There are six columns of data:
#' \itemize{
#'  \item{\code{"species"}} {The species' name (Latin binomial)}
#'  \item{\code{"year"}} {year to which the index value refers (integer)}
#'  \item{\code{"occupancy"}} {A national index of occupancy, defined as the proportion of occupied 1 km grid cells.}
#'  \item{\code{"index"}} {The national index of occupancy on the logit (log odds) scale.}
#'  \item{\code{"se"}} {Standard error of the posterior distribution of logit occupancy}
#'  \item{\code{"inAtlas"}} {Logical. Identifies the 38 the species was included in the 2014 Atlas of Dragonflies in Britain and Ireland.} 
#'  }
#' @references Charlotte L. Outhwaite, Gary D. Powney, Tom A. August, Richard E. Chandler, Stephanie Rorke, Oliver L. Pescott, Martin Harvey, Helen E. Roy, Richard Fox, David B. Roy, Keith Alexander, Stuart Ball, Tristan Bantock, Tony Barber, Björn C. Beckmann, Tony Cook, Jim Flanagan, Adrian Fowles, Peter Hammond, Peter Harvey, David Hepper, Dave Hubble, John Kramer, Paul Lee, Craig MacAdam, Roger Morris, Adrian Norris, Stephen Palmer, Colin W. Plant, Janet Simkin, Alan Stubbs, Peter Sutton, Mark Telfer, Ian Wallace & Nick J. B. Isaac. (2019). Annual estimates of occupancy for bryophytes, lichens and invertebrates in the UK, 1970-2015. Scientific Data, 6(1), 259. \url{https://doi.org/10.1038/s41597-019-0269-1}
#' @source \url{https://doi.org/10.1038/s41597-019-0269-1}
#' @examples 
#' data(dragonflies)
#' head(dragonflies)
"dragonflies"