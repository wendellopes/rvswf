#' Reference function Derivative of Ricatti-Bessel function .
#' 
#' @details The derivative of the Ricatti-Bessel function \eqn{\psi_n(x)=xj_n(x)}.
#' @param x The argument of the function
#' @param n The order of the function
#' @return The derivative \eqn{\psi_n'(x)=j_n(x)+xj_n'(x)}.
#' @import reff.rjn
#' @export
reff.rdj<-function(n,x){
   a<-n/(2*n+1)
   b<-(n+1)/(2*n+1)
	return(b*reff.rjn(x,n-1)-a*reff.rjn(x,n+1))
}
