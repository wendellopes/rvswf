#' Lorentz-Mie Optical Force Calculation.
#' @param m The ratio between the refractive indices.
#' @param x The form factor value.
#' @param by \code{LD} using logarithmic derivative; \code{RB} using the ratio
#' between Ricatti-Bessel functions.
#' @return The coefficients \eqn{A_n}, \eqn{B_n} and \eqn{C_n} used to calculate
#' the Optical force.
lmie.ofc<-function(m,x){
   k<-lmie.exp(m,x,by="LD")
   n0<-1:(nrow(k)-1)
   n1<-n0+1
   An<-k$an[n1]+Conj(k$an[n0])-2*k$an[n1]*Conj(k$an[n0])
   Bn<-k$bn[n1]+Conj(k$bn[n0])-2*k$bn[n1]*Conj(k$bn[n0])
   Cn<-k$an[n0]+Conj(k$bn[n0])-2*k$an[n1]*Conj(k$an[n0])
   return(data.frame(An,Bn,Cn))
}