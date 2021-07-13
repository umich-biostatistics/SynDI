
####################################
########## Function 1: Expit function
###################################

#' Expit function
#'
#' Expit helper function.

expit <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}
