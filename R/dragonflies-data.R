#' @name dragonflies
#' @title Data - British Dragonfly Society, UK
#' @description This dataset from the British Dragonfly Society, UK, has data for
#' UK dragonfly species
#' The data contains: species, latin names for drangonfly species; year, integer
#' from 1970-2015; index, the index of occupancy (the proportion of occupied
#' sites). This value is the mean of the posterior distribution; index_logit,
#' the index of occupancy on the logit scale. This value is the mean of the
#' posterior distribution; sd standard deviation of the posterior distribution
#' of logit occupancy; inAtlas, logical, Was a trend estimate for
#' this species included in the 2014 Atlas of dragonflies?
#' @docType data
#' @usage data(dragonflies)
#' @references Charlotte L. Outhwaite, Gary D. Powney, Tom A. August, Richard E. Chandler, Stephanie Rorke, Oliver L. Pescott, Martin Harvey, Helen E. Roy, Richard Fox, David B. Roy, Keith Alexander, Stuart Ball, Tristan Bantock, Tony Barber, Björn C. Beckmann, Tony Cook, Jim Flanagan, Adrian Fowles, Peter Hammond, Peter Harvey, David Hepper, Dave Hubble, John Kramer, Paul Lee, Craig MacAdam, Roger Morris, Adrian Norris, Stephen Palmer, Colin W. Plant, Janet Simkin, Alan Stubbs, Peter Sutton, Mark Telfer, Ian Wallace & Nick J. B. Isaac. (2019). Annual estimates of occupancy for bryophytes, lichens and invertebrates in the UK, 1970-2015. Scientific Data, 6(1), 259. \url{https://doi.org/10.1038/s41597-019-0269-1}
#' @source \url{https://doi.org/10.1038/s41597-019-0269-1}
#' @examples 
#' data(dragonflies)
#' head(dragonflies)
"dragonflies"