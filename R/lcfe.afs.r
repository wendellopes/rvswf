#' Auxiliary function \eqn{S_n(x)=n/x}. Used on calculations of Bessel functions.
#' 
#' @param n The order of the Bessel function.
#' @param x The argument of the Bessel function.
#' @return The value \code{n/x}.
#' @export
#' @examples
#' u<-0:100
#' n<-sample(u,10)
#' x<-sample(u,10)
#' print(data.frame(n,x,"n/x"=n/x,"lcfe.afs"=lcfe.afs(n,x)))
lcfe.afs<-function(n,x){
   return(n/x)
}
