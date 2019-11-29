#' Extract BMA BUGS code
#' 
#' @description Gets a copy of the BUGS code and writes it 
#' @param option text string defining which particular variant is desired 
#' @param print.screen Logical, should the code be printed to the console?
#' @param save.local Logical, should text file of the BUGS code be produced?
#' @export
#' @examples
#' get_bmaBUGScode(option="smooth", print.screen=TRUE, save.local=FALSE)

get_bmaBUGScode <- function(option="smooth", 
                                print.screen = FALSE,
                                save.local = FALSE) {
  
  switch(tolower(option),
         smooth = {model_code <- bma_model_Smooth()},
         smooth_det2 = {model_code <- bma_model_smooth_det2()},
         smooth_det_sigtheta = {model_code <- bma_model_smooth_det_sigtheta()},
         smooth_det = {model_code <- bma_model_smooth_det()},
         random_walk = stop("Random walk model has been deprecated"),
         uniform = stop("Uniform model has been deprecated"),
         uniform_noeta = stop("Uniform model has been deprecated"),
         fngr = stop("This model option has been deprecated"),
         smooth_stoch = stop("This model option has been deprecated"),
         fngr2 = stop("This model option has been deprecated"),
         smooth_stoch2 = stop("This model option has been deprecated"),
         {stop(paste("Model type not known. Must be one of 'smooth',",
                     "'smooth_det', 'smooth_det2', 'smooth_det_sigtheta'"))})

  if(print.screen) cat(model_code, sep = '\n')
  if(save.local) writeLines(text = model_code, con = "bma_BUGS_code.txt")
  return(model_code)
}


