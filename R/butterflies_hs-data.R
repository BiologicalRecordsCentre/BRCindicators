#' @name butterflies_hs
#' @title Data - UK Butterfly Monitoring Scheme - Habitat Specialists
#' @description This dataset from the UK Butterfly Monitoring Scheme has national abundance indices for 26 species regarded as habitat specialists from 1976-2017.
#' There are 1000 rows of data. Ten species have a complete time-series (42 years of data); 15 species join late (9 in the 1970s, 3 in the 1980s and 3 in the 1990s); the Swallowtail (Papilio machaon) spans the full range of years but has no data for 1978.
#' @docType data
#' @format 
#' There are four columns of data:
#' \itemize{
#'  \item{\code{"species"}} {The species' name (Latin binomial)}
#'  \item{\code{"year"}} {year to which the index value refers (integer)}
#'  \item{\code{"index"}} {A "collated" index of abundance on the log10 scale. The time-series for each species is fixed to have a mean of exactly 2. Since 10^2 = 100, the time-series are centered on a value of 100.}
#'  \item{\code{"se"}} {standard errors on the index value (on the log10 scale), as estimated by the underlying statistical model.}
#'  }
#' @usage data(butterflies_hs)
#' @references Botham, M.S.; Brereton, T.; Harris, S.; Harrower, C.; Middlebrook, I.; Randle, Z.; Roy, D.B. (2019). United Kingdom Butterfly Monitoring Scheme: collated indices 2017. NERC Environmental Information Data Centre. \url{https://catalogue.ceh.ac.uk/documents/ace3c3ef-df89-40b9-ba8b-106997fd6d9c}
#' @source \url{https://catalogue.ceh.ac.uk/documents/ace3c3ef-df89-40b9-ba8b-106997fd6d9c}
#' @examples 
#' data(butterflies_hs)
#' head(butterflies_hs)
"butterflies_hs"