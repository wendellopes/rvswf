#' Ricatti-Bessel function
#' 
#' @details Only a parser to the built-in \code{R} function.
#' @param n The order of the function
#' @param x The argument of the funcion
#' @return The Ricatti-Bessel function \eqn{\psi_n(x)=xj_n(x)}
#' @export
reff.rjn<-function(n,x){
	return(sqrt(pi/2)*besselJ(x,n+.5)*sqrt(x))
}
