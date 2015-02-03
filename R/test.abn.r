#' Test the calculation of Mie Coefficients 
#' 
#' @detail Must be rewritten
#' @param n The order of the function
#' @param m The relative index of refraction
#' @param x The argument of the function
#' @return Data frame with values \eqn{A_n}, \eqn{B_n} and \eqn{C_n}.
#' @import test.abc
test.abn<-function(n,m,x){
   u.1<-test.abc(n,x)
   u.m<-test.abc(n,m*x)
   Ta<-(u.m$An/m-u.1$An)/(u.m$An/m-u.1$Bn)
   Tb<-(u.m$An*m-u.1$An)/(u.m$An*m-u.1$Bn)
   an<-u.1$Cn*Ta
   bn<-u.1$Cn*Tb
   return(data.frame(An=u.1$An,Bn=u.1$Bn,Cn=u.1$Cn,Ta,Tb,an,bn))
}
