#' Compare results for Ricatti-Bessel Functions.
#' 
#' @details Compare results using vswf, built in \code{R} and \code{gsl} algorithms.
#' @param x The argument of \eqn{\psi_n(x)=xj_n(x)}.
#' @param n The order of the Ricatti-Bessel function.
#' @return Table comparing built-in \code{R} functions, \code{gsl} and \code{rvswf}.
#' @include reff.sjn.r bess.ric.r
#' @export
#' @examples
#' x<-5
#' nmax<-10
#' print(comp.ric(5,3))
comp.ric<-function(nmax,x){
   u<-x*reff.sjn(x,0:nmax)
   v<-bess.ric(nmax,x)$Rn
   w<-x*bessel_jl_array(nmax,x)
   t<-x*bessel_jl_steed_array(nmax,x)
   return(data.frame(R=u,VSWF=v,GSL=w,GSL.STEED=t))
}
