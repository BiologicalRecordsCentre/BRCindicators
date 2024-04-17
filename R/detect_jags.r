#' @import runjags

detect_jags <- function(){
    return(runjags::testjags(silent = TRUE)$JAGS.available)
}