#' Extract BMA BUGS code
#' 
#' @description Gets a copy of the BUGS code and writes it to the screen and (optionally) to a text file
#' @param option text string defining which particular variant is desired 
#' @param print.screen Logical, should the code be printed to the console?
#' @param save.local Logical, should text file of the BUGS code be produced?
#' @param incl.2deriv Logical, Should there be a calculation of the indicator's second derivative?
#' @param Y1perfect Logical, should the first year for each species be assumed to be known without error (`TRUE`)
#' @param seFromData Logical, Should the standard errors be read in from data (`TRUE`) or estimated (`FALSE`)
#' @details There are a number of model to choose from:
#' \itemize{
#'  \item{\code{"smooth"}}{ The default option. Indicator defined by Growth rates, with Ruppert smoother, allowing for species to join late. Error on the first year of each species' time-series is assumed to be zero. The indicator is the expected value of the geometric mean across species (with missing species imputed). 
#'  Includes three options: `seFromData` `Y1perfect` and `incl.2deriv`. See bayesian_meta_analysis for mode details.}
#'  \item{\code{"smooth_JABES"}}{ Equivalent to smooth with `seFromData = TRUE` and `Y1perfect = TRUE`. This is the version implemented in the JABES paper. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  \item{\code{"smooth_det2"}}{ Equivalent to smooth with `seFromData = TRUE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  \item{\code{"smooth_det_sigtheta"}}{ Equivalent to smooth with `seFromData = FALSE` and `Y1perfect = FALSE`. Retained for backwards compatability. Choosing this option will overwrite user-entered options for `seFromData` and `Y1perfect`.}
#'  \item{\code{"smooth_det"}}{ Specific variant of smooth_det2 - under review. Likely to be deprecated}
#'  }
#' @export
#' @examples
#' get_bmaBUGScode(option="smooth", print.screen=TRUE, save.local=FALSE)
#' get_bmaBUGScode(option="smooth_JABES", print.screen=TRUE, save.local=FALSE)
#' get_bmaBUGScode(option="smooth_det", print.screen=TRUE, save.local=FALSE)
#' get_bmaBUGScode(option="smooth_det2", print.screen=TRUE, save.local=FALSE)
#' get_bmaBUGScode(option="smooth_det_sigtheta", print.screen=TRUE, save.local=FALSE)

get_bmaBUGScode <- function(option="smooth", 
                            print.screen = FALSE,
                            save.local = FALSE,
                            incl.2deriv = FALSE,
                            Y1perfect = FALSE,
                            seFromData = FALSE){
  
  switch(tolower(option),
         smooth = {model_code <- bma_model_Smooth(incl.2deriv, 
                                                  seFromData,
                                                  Y1perfect)},
         smooth_jabes = {model_code <- bma_model_Smooth(# this is the version SF ran for the paper
                                                  incl.2deriv, 
                                                  seFromData = FALSE,
                                                  Y1perfect = TRUE)},
         smooth_det2 = {model_code <- bma_model_Smooth(# this is the version NI tested for ISEC 
                                                  incl.2deriv, 
                                                  seFromData = TRUE,
                                                  Y1perfect = FALSE)},
         smooth_det_sigtheta = {model_code <- bma_model_Smooth(# this is the version NI ran for Scottish indicators
                                                  incl.2deriv, 
                                                  seFromData=FALSE,
                                                  Y1perfect = FALSE)},
         smooth_det = {model_code <- bma_model_smooth_det()},
         # Deprecated options can still be called directly from this function
         random_walk = model_code <- bma_model_ranwalk(),
         uniform = model_code <- bma_model_uniform(),
         uniform_noeta = model_code <- bma_model_uniform_noeta(),
         fngr = model_code <- bma_model_FNgr(),
         smooth_stoch = model_code <- bma_model_smooth_stoch(),
         fngr2 = model_code <- bma_model_FNgr2(),
         smooth_stoch2 = model_code <- bma_model_smooth_stoch2(),
         {stop(paste("Model type not known. Check the help file for details"))})

  if(print.screen) cat(model_code, sep = '\n')
  if(save.local) writeLines(text = model_code, con = paste0("bma_BUGS_code_",option,".txt"))
  return(model_code)
}


