#' Refference function Derivative of the Spherical Bessel function.
#' 
#' @details Uses the definition and built-in \code{R} function to give the
#' @param x The argument of the function
#' @param n The order of the function
#' @return The value of the derivative
#' @export
reff.cdj<-function(n,x){
   return(.5*(besselJ(x,n-1)-besselJ(x,n+1)))
}
