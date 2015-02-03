#' Compare results for Spherical Bessel Functions.
#' 
#' @details Compare results using vswf, built in \code{R} and \code{gsl} algorithms.
#' @param x The argument of \eqn{j_n(x)}.
#' @param n The order of the Spherical Bessel function.
#' @return Table comparing built-in \code{R} functions, \code{gsl} and \code{rvswf}.
#' @import reff.sjn, bess.sph
#' @export
#' @examples
#' x<-5
#' nmax<-10
#' print(comp.sph(5,3))
comp.sph<-function(nmax,x){
   u<-reff.sjn(x,0:nmax)
   v<-bess.sph(nmax,x)$jn
   w<-bessel_jl_array(nmax,x)
   t<-bessel_jl_steed_array(nmax,x)
   return(data.frame(R=u,VSWF=v,GSL=w,GSL.STEED=t))
}

