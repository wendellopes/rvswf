#' Compare results for Cylindrical Bessel Functions.
#' 
#' @details Compare results using vswf, built in \code{R} and \code{gsl} algorithms.
#' @param x The argument of \eqn{J_n(x)}.
#' @param n The order of the Cylindrical Bessel function.
#' @return Table comparing built-in \code{R} functions, \code{gsl} and \code{rvswf}.
#' @import bess.cyl
#' @export
#' @examples
#' x<-5
#' nmax<-10
#' print(comp.cyl(5,3))
comp.cyl<-function(nmax,x){
   u<-besselJ(x,0:nmax)
   v<-bess.cyl(nmax,x)$Jn
   w<-bessel_Jn_array(0,nmax,x)
   return(data.frame(R=u,VSWF=v,GSL=w))
}