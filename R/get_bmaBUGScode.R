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
         random_walk = {model_code <- bma_model_ranwalk()},
         uniform = {model_code <- bma_model_uniform()},
         uniform_noeta = {model_code <- bma_model_uniform_noeta()},
         fngr = {model_code <- bma_model_FNgr()},
         smooth_stoch = {model_code <- bma_model_smooth_stoch()},
         smooth_det = {model_code <- bma_model_smooth_det()},
         fngr2 = {model_code <- bma_model_FNgr2()},
         smooth_stoch2 = {model_code <- bma_model_smooth_stoch2()},
         smooth_det2 = {model_code <- bma_model_smooth_det2()},
         smooth_det_sigtheta = {model_code <- bma_model_smooth_det_sigtheta()},
         {stop(paste("model type not know. Must be one of 'smooth', 'random_walk',",
                     "'uniform', 'uniform_noeta', 'FNgr', 'smooth_stoch',",
                     "'smooth_det', 'FNgr2', 'smooth_stoch2', 'smooth_det2', 'smooth_det_sigtheta'"))})

  if(print.screen) cat(model_code, sep = '\n')
  if(save.local) writeLines(text = model_code, con = "bma_BUGS_code.txt")
  return(model_code)
}


