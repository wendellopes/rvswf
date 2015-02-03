#' Beam Shape Coefficients for a Plane Wave, as done by Mie
#' 
#' @details We use the expression showed on Jackson's Classical Electrodynamics
#' book.
#' @param lmax The maximum value of \eqn{l}.
#' @param norm If TRUE, the Beam Shape Coefficient will be divided by \eqn{\sqrt{2}}.
#' @param s The polarity of the plane wave.
#' @return Beam Shape Coefficients for a Mie Plane Wave (z-direction).
#' @include vswf.jlm.r
#' @export
#' @seealso \code{\link{vswf.gpw}}, \code{\link{vswf.jlm}}.
#' @examples
#' lm<-5
#' a<-vswf.mpw(lm)
#' plot(Re(a$GTE),type='b')
#' points(Im(a$GTE),pch=4,col='red',type='b')
#' plot(Re(a$GTM),type='b')
#' points(Im(a$GTM),pch=4,col='red',type='b')
vswf.mpw<-function(lmax,norm=TRUE,s=1){
#-------------------------------------------------------------------------------
   if(lmax<1){lmax<-1}                         # Pelo menos 1 termo
   LMAX=lmax*(lmax+2)+1                        # Vetor para lmax
#-------------------------------------------------------------------------------
   gte<-rep(0,LMAX)
   gtm<-rep(0,LMAX)
   glm<-function(l,m,norm=TRUE){
      k<-(1i^l)*sqrt(4*pi*(2*l+1))
      k<-rep(k,2*l+1)
      k[m!=s]<-0
      if(norm){
         k<-k/sqrt(2)
      }
      return(k)
   }
   for(l in 0:lmax){
      m<--l:l
      gte[vswf.jlm(l,m)]<-glm(l,m,norm)
   }
   gtm<--1i*s*gte
   U<-data.frame(GTE=gte,GTM=gtm)
   return(U)
}
