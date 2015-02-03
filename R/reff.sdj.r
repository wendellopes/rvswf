#' Reference function Derivative of Spherical Bessel function .
#' 
#' @details The derivative of the Spherical Bessel function \eqn{j_n(x)}.
#' @param x The argument of the function
#' @param n The order of the function
#' @return The derivative \eqn{j_n'(x)}.
#' @import reff.sjn
#' @export
reff.sdj<-function(x,n){
   a<-n/(2*n+1)
   b<-(n+1)/(2*n+1)
	return(a*reff.sjn(x,n-1)-b*reff.sjn(x,n+1))
}
