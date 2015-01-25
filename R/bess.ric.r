#' Calculates Ricatti-Bessel functions from 0 to nmax.
#'
#' @details \code{bess.ric} calculates the Ricatti-Bessel
#' functions using downward recurrence, from \eqn{R_nmax(x)} to \eqn{R_0(x)}.
#' The system of equations is given by \eqn{S_n(x)=n/x}, 
#' \eqn{\rho_n=R_n(x)/R_{n+1}(x)}{r[n]=R_n/R_{n+1}} and 
#' \eqn{C_n=R_n'(x)/R_n(x)}. The system can be solved by means of
#' the recurrence relations of the Ricatti-Bessel functions
#' \deqn{    \rho_{n-1}+\frac{1}{\rho_n}=S_{2n+1}   }{r[n-1]+1/r[n]=S[2n+1]}
#' \deqn{(n+1\rho_{n-1}-\frac{n}{\rho_n}=(2n+1)C_{n}}{r[n-1]-n/r[n]=(2n+1)C[n]}
#' that can be rewriten
#' \deqn{\rho_{n}=S_{n+1}+C_{n+1}                }{  r[n]=S[n+1]+C[n+1]}
#' \deqn{\frac{1}{\rho_n}=S_{n+1}-C_n.           }{1/r[n]=S[n+1]-C[n  ].}
#' The logarithmic derivatives obeys the relation,
#' \deqn{(S_{n+1}-C_n)(S_{n+1}+C_{n+1})=1.         }{(S[n+1]-C[n])(S[n+1]+C[n+1])=1.}
#' The values can be calculated upward or downward.
#' @param nmax The maximum order of \eqn{R_n(x)}
#' @param x The argument of \eqn{R_n(x)}
#' @param code If you prefer to use native R or C language.
#' The algorithm is the same.
#' @return An array of Ricatti-Bessel functions and its derivatives.
#' from 0 to \code{nmax} at point \code{x}
#' @examples
#' x<-30
#' nmax<-50
#' a<-bess.ric(nmax,x,code="C")
#' b<-bess.ric(nmax,x,code="R")
#' d<-sqrt(x*pi/2)*besselJ(x=x,nu=.5+(0:nmax))
#' plot(a$Rn,type='b')
#' points(b$Rn,col='red',pch=4)
#' points(d,col='blue',pch=3)
bess.ric<-function(nmax,x,code="C"){
   if(abs(x)<1e-10){
      return(data.frame(Rn=rep(0,nmax+1),dRn=rep(0,nmax+1)))
   }
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("bess_ric",
            nmax=as.integer(nmax),
            x=as.double(x),
            Rn=as.double(dummy),
            Dn=as.double(dummy))
      return(data.frame(Rn=u$Rn,dRn=u$Dn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   Cn<-rep(0,nmax+1)  # Vector 
	   rn<-rep(1,nmax+1)  # Vector
	   rm<-rn
	   Cn[nmax+1]<-lcfe.rbl(nmax,x) # Last element
	   rn[nmax+1]<-lcfe.sbd(nmax,x) # Last element
	   Cm<-Cn
	   # DOWNWARD RECURRENCE
	   RN<-1
      for(n in nmax:1){
      	  # original
      	  rn[n]<-S(n,x)+Cn[n+1]
      	  Cn[n]<-S(n,x)-1/rn[n]
      	  # modified (permits one step normalization)
      	  rm[n]<-S(n,x)*rm[n+1]+Cm[n+1]
      	  Cm[n]<-S(n,x)*rm[n]-rm[n+1]
      	  # Normalization
      	  if(abs(rm[n])>1e10){
      	  	 #cat("renorming...\n")
      	  	 #print(c(rn[n-1],Cn[n-1]))
      	  	 Cm<-Cm/rm[n] # this must be done first
      	  	 rm<-rm/rm[n] # otherwise the result will be wrong.
      	  }
      	  #print(c(rn[n-1],Cn[n-1]))
      }
      # one step normalization taking care about zeros
      # Bessel function
      j0<-bess.szr(x)
      j1<-bess.sun(x)
      R0<-x*j0
      R1<-x*j1
      dR0<-cos(x)
      dR1<--j1+R0
      if(abs(R0)>1e-10){ # R_0 != 0
         jn<-(rm/rm[1])*R0 # create functions for normalizations
      }else{
         jn<-(rm/rm[2])*R1
      }
      # Its Derivative
      if(dR0>1e-10){ # R_0' != 0
      	  djn<-(Cm/Cm[1])*dR0
      }else{
      	  djn<-(Cm/Cm[2])*dR1
      }
      # Return results
      return(data.frame(Rn=jn,dRn=djn))
   }
}
