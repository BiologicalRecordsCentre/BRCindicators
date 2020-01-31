#' @name butterflies_wc
#' @title Data - UK Butterfly Monitoring Scheme - Wider Countryside Species
#' @description This dataset from the UK Butterfly monitoring scheme has data for 24 species regarded as wider countryside species.
#' The data contains: species, Latin name of the butterfly species; year, integer
#' from 1976-2017; index, the index of abundance on the log10 scale. The
#' time-series for each species is fixed to have a mean of exactly 2.
#' Since 10^2 = 100, the time-series are centered on a value of 100; se, standard
#' error in the index value (on the log scale)
#' @docType data
#' @usage data(butterflies_wc)
#' @references Botham, M.S.; Brereton, T.; Harris, S.; Harrower, C.; Middlebrook, I.; Randle, Z.; Roy, D.B. (2019). United Kingdom Butterfly Monitoring Scheme: collated indices 2017. NERC Environmental Information Data Centre. \url{https://catalogue.ceh.ac.uk/documents/ace3c3ef-df89-40b9-ba8b-106997fd6d9c}
#' @source \url{https://catalogue.ceh.ac.uk/documents/ace3c3ef-df89-40b9-ba8b-106997fd6d9c}
#' @examples 
#' data(butterflies_wc)
#' head(butterflies_wc)
"butterflies_wc"