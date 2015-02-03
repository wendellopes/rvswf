#' Normalized Legendre Functions and its derivatives
#' 
#' @details Calculate the normalized Legendre Functions \eqn{Q_l^m}, in the 
#' way that \eqn{Y_{lm}(\theta,\phi)=Q_l^m(\cos\theta) e^{im\phi}}.
#' @param cth The argument of \eqn{Q_l^m(x)}, so \eqn{-1<x=\cos\theta<1}.
#' @param lmax The maximum value of \eqn{l}.
#' @return An array of \eqn{Q_l^m} and its derivative \eqn{Q_l^m{}'}.
#' @include vswf.jlm.r
#' @export
#' @examples
#' u<-vswf.qlm(.5,5)
#' print(u)
vswf.qlm<-function(cth,lmax){
   sth<-sqrt(1-cth^2)
#-------------------------------------------------------------------------------
   LMAX=lmax*(lmax+2)+1
#-------------------------------------------------------------------------------
   alfaQ<-function(l,m){
      return(sqrt(((2*l-1)*(2*l+1))/((l-m)*(l+m))))
   }
   betaQ<-function(l,m){
      return(sqrt((2*l+1)/(2*l-3))*sqrt(((l+m-1)*(l-m-1))/((l-m)*(l+m))))
   }
   gammaQ<-function(l){
      return(sqrt((2*l+1)/(2*l)))
   }
   deltaQ<-function(l){
      return(sqrt(2*l+1))
   }
#-------------------------------------------------------------------------------
   if(lmax<2){
      LIST<-4
   }else{
      LIST<-LMAX
   }
   ls  <-rep(0,LIST)
   ms  <-rep(0,LIST)
   Qlm <-rep(0,LIST)
   dQlm<-rep(0,LIST)
#-------------------------------------------------------------------------------
   ls[1:4]<-c(0, 1,1,1)
   ms[1:4]<-c(0,-1,0,1)
#-------------------------------------------------------------------------------
   Qlm[vswf.jlm(0,0)]<-1/sqrt(4*pi)                      #Q00
   Qlm[vswf.jlm(1,1)]<--gammaQ(1)*sth*Qlm[vswf.jlm(0,0)] #Q11
   Qlm[vswf.jlm(1,0)]<-sqrt(3)*cth*Qlm[1]                #Q10
   Qlm[vswf.jlm(1,-1)]<--Qlm[vswf.jlm(1,1)]              #Q11*(-1)
#-------------------------------------------------------------------------------
   dQlm[vswf.jlm(0, 0)]<-0                                              # dQ00
   dQlm[vswf.jlm(1, 1)]<--(cth/sth^2)*Qlm[vswf.jlm(1, 1)]               # dQ11
   dQlm[vswf.jlm(1, 0)]<--(cth*Qlm[vswf.jlm(1, 0)]-
      sqrt(3)*Qlm[vswf.jlm(0,0)])/(sth^2)                               # dQ10
   dQlm[vswf.jlm(1,-1)]<--(cth/sth^2)*Qlm[vswf.jlm(1,-1)]               # dQ11
#-------------------------------------------------------------------------------
   if(lmax>1){
      for(l in 2:lmax){
         m.p<-0:(l-2)
         m.m<-(-l):(-1)
         mm<-(-l):l
         ls[vswf.jlm(l,mm)]<-rep(l,2*l+1)
         ms[vswf.jlm(l,mm)]<-mm
         csc2<-1/(sth^2)
         csp<-(-1)^m.m
         Qlm[vswf.jlm(l,l  )]=-gammaQ(l)*sth*Qlm[vswf.jlm(l-1,l-1)]     #OK
         Qlm[vswf.jlm(l,l-1)]=deltaQ(l)*cth*Qlm[vswf.jlm(l-1,l-1)]      #OK
         Qlm[vswf.jlm(l,m.p)]=alfaQ(l,m.p)*cth*Qlm[vswf.jlm(l-1,m.p)]-
            betaQ(l,m.p)*Qlm[vswf.jlm(l-2,m.p)]                         #OK 
         Qlm[vswf.jlm(l,m.m)]=csp*Qlm[vswf.jlm(l,abs(m.m))]             #OK
         dQlm[vswf.jlm(l,mm)]<--l*cth*csc2*Qlm[vswf.jlm(l,mm)]+
            sqrt((2*l+1)/(2*l-1))*
            sqrt((l+mm)*(l-mm))*csc2*Qlm[vswf.jlm(l-1,mm)]
      }
   } 
   return(data.frame(l=ls,m=ms,Qlm,dQlm))
}
