#' Lorentz-Mie Optical Force Calculation.
#' @param m The ratio between the refractive indices.
#' @param x The form factor value.
#' @param by \code{LD} using logarithmic derivative; \code{RB} using the ratio
#' between Ricatti-Bessel functions.
#' @return The coefficients \eqn{A_n}, \eqn{B_n} and \eqn{C_n} used to calculate
#' the Optical force.
#' @import lmie.exp
#' @export
#' @examples
#' lmax<-20
#' lambda=.532  # Green 532 nm
#' m<-ridx.pol(lambda)/ridx.h2o(lambda)
#' g<-vswf.mpw(lmax)
#' 
lmie.ofc<-function(m,x,GTE,GTM,lmax=floor(x+2+4*x^(1/3))){
   u<-lmie.exp(m,x,by="LD",NMAX=lmax)   
   nu<-nrow(u)
   nm<-nu*(nu+2)
   Anm<-Bnm<-Cnm<-l<-m<-j<-1:nm
   # Calculations
   n0<-1:(nu-1)
   n1<-n0+1
   An<-u$an[n1]+Conj(u$an[n0])-2*u$an[n1]*Conj(u$an[n0])
   Bn<-u$bn[n1]+Conj(u$bn[n0])-2*u$bn[n1]*Conj(u$bn[n0])
   Cn<-u$an[n0]+Conj(u$bn[n0])-2*u$an[n1]*Conj(u$an[n0])
   for(k in 1:nu){
      s<--k:k
      n<-k*(k+1)+s
      l[n]<-k
      m[n]<-s
      j[n]<-n
      Anm[n]<-An[k]
      Bnm[n]<-Bn[k]
      Cnm[n]<-Cn[k]
   }
   return(data.frame(Anm,Bnm,Cnm))
}