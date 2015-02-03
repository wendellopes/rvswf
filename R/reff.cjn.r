#' Cylindrical Bessel function
#' 
#' @details Only a parser to the built-in \code{R} function.
#' @param n The order of the function
#' @param x The argument of the funcion
#' @return The Cylindrical Bessel function
#' @export
reff.cjn<-function(n,x){
   return(besselJ(x,n))
}
