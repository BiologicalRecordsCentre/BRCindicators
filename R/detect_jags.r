#' Detect JAGS installation
#' 
#' @description small function that detects whether JAGS is installed.
#' @return TRUE or FALSE, indicating whether the JAGS installation has been detected 
#' @importFrom runjags testjags
#' @export
#' @examples
#' 
#' detect_jags()

#' @importFrom runjags testjags

detect_jags <- function(){
    return(suppressWarnings(runjags::testjags(silent = TRUE)$JAGS.found))
}