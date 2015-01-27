#' Auxiliary function \eqn{S_n(x)=n/x}. Used on calculations of Bessel functions.
#' 
#' @param n The order of the Bessel function.
#' @param x The argument of the Bessel function.
#' @return The value \code{n/x}.
lcfe.afs<-function(n,x){
   return(n/x)
}
