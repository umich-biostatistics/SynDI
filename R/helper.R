
####################################
########## Function 1: Expit function
###################################

#' Expit function
#' 
#' @param x vector to expit
#' 
#' @return numeric vector with the value of the expit function 
#' y = expit(x) = exp(x)/(1+exp(x)).
#'
#' Expit helper function.

expit <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}
