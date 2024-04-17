#' @importFrom runjags testjags

detect_jags <- function(){
    return(runjags::testjags(silent = TRUE)$JAGS.available)
}