#' General expression for wave guides.
#' 
#' @details This expression does not take in account the value of the Fourier 
#' transform of azimuthal component. It does not depend on the position of the
#' expansion, only on the wave vector.
#' @param gamma The transversal component of the wave vector \eqn{\gamma}.
#' @param kz The longitudinal component of the wave vector \eqn{k_z}.
#' @param lmax The maximum value of \eqn{l} of the expansion.
#' @return The values \eqn{A} and \eqn{B}.
#' @include vswf.qlm.r
#' @export
#' @seealso \code{\link{vswf.cwg}}, \code{\link{vswf.rwg}}.
#' @examples 
#' th<-pi/3
#' u<-vswf.gwg(sin(th),cos(th),10)
#' plot(Re(u$A),type='b')
#' points(Im(u$A),type='b',pch=4,col='red')
#' plot(Re(u$B),type='b')
#' points(Im(u$B),type='b',pch=4,col='red')
vswf.gwg<-function(gama,kz,lmax){
   k<-sqrt(gama^2+kz^2)
   LMAX=lmax*(lmax+2)+1
   #----------------------------------------
   u<-vswf.qlm(kz/k,lmax)
   Qlm<-u$Qlm
   dQlm<-u$dQlm
   ll<-u$l
   mm<-u$m
   llp1<-1/sqrt(ll*(ll+1))
   llp1[1]<-0
   #----------------------------------------
   A<-2*(1i^ll)*((k/gama)^2)*Qlm*mm*llp1
   B<-2*(1i^(ll-1))*dQlm*llp1
   return(data.frame(A,B))
}
