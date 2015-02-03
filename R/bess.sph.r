#' Calculates Spherical Bessel functions from 0 to nmax.
#'
#' @details \code{bess.sph} calculates the Spherical Bessel
#' functions using downward recurrence, from \eqn{j_nmax(x)} to \eqn{j_0(x)}.
#' The system of equations is given by \eqn{S_n(x)=n/x}, 
#' \eqn{\rho_n=j_n(x)/j_{n+1}(x)}{r[n]=j_n/j_{n+1}} and 
#' \eqn{c_n=j_n'(x)/j_n(x)}. The system can be solved by means of
#' the recurrence relations of the Spherical Bessel functions
#' \deqn{ \rho_{n-1}+\frac{1  }{\rho_n}=S_{2n+1}   }{ r[n-1]+    1/r[n]=S[2n+1]}
#' \deqn{n\rho_{n-1}-\frac{n+1}{\rho_n}=(2n+1)c_{n}}{nr[n-1]-(n+1)/r[n]=(2n+1)c[n]}
#' that can be rewriten
#' \deqn{\rho_n=S_{n+2}+c_{n+1}          }{  r[n]=S[n+2]+c[n+1]}
#' \deqn{\frac{1}{\rho_n}=S_n-c_n.       }{1/r[n]=S[n  ]-c[n  ].}
#' The logarithmic derivatives obeys the relation,
#' \deqn{(S_{n+2}-c_{n+1})(S_n+c_n)=1.   }{(S[n+2]-c[n])(S[n]+C[n])=1.}
#' The values can be calculated upward or downward.
#' @param nmax The maximum order of \eqn{j_n(x)}
#' @param x The argument of \eqn{j_n(x)}
#' @param code If you prefer to use native R or C language.
#' The algorithm is the same.
#' @return An array of Spherical Bessel functions and its derivatives 
#' from 0 to \code{nmax} at point \code{x}, and also the logarithmic
#' derivative \eqn{c_n=j_n'/j_n} and the ratio \eqn{\rho_n=j_n/j_{n+1}}.
#' @useDynLib rvswf
#' @import lcfe.sbl, lcfe.sbd
#' @export
#' @examples
#' x<-30
#' nmax<-50
#' a<-bess.sph(nmax,x,code="C")
#' b<-bess.sph(nmax,x,code="R")
#' d<-sqrt(pi/(2*x))*besselJ(x=x,nu=.5+(0:nmax))
#' plot(a$jn,type='b')
#' points(b$jn,col='red',pch=4)
#' points(d,col='blue',pch=3)
#-------------------------------------------------------------------------------
bess.sph<-function(nmax,x,code="C"){
   if(abs(x)<1e-10){
      return(data.frame(Rn=rep(0,nmax+1),Dn=rep(0,nmax+1)))
   }
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(1,nmax+1)
      u<-.C("bess_sph",
            nmax=as.integer(nmax),
            x=as.double(x),
            jn=as.double(dummy),
            dn=as.double(dummy),
            NMAX=as.integer(2000))
      return(data.frame(jn=u$jn,djn=u$dn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   dn<-rep(0,nmax+1)  # Vector 
	   jn<-rep(1,nmax+1)  # Vector
	   dn[nmax+1]<-lcfe.sbl(nmax,x) # Last element
	   # DOWNWARD RECURRENCE
      for(n in nmax:1){
      	jn[n]<-S(n+1,x)*jn[n+1]+dn[n+1]
      	dn[n]<-S(n-1,x)*jn[n  ]-jn[n+1]
      	if(abs(jn[n])>1e2){
      	  	dn<-dn/jn[n] # this must be done first
      	   jn<-jn/jn[n] # otherwise the result will be wrong.
      	}
      }
      # Bessel function
      b0<-bess.szr(x)
      b1<-bess.sun(x)
      if(abs(b0)>1e-10){   # j_0 != 0
         jn<-(jn/jn[1])*b0 # create functions for nojnalizations
      }else{               # j_0 == 0
         jn<-(jn/jn[2])*b1
      }
      # Its Derivative
      if(abs(b1)>1e-10){         # j_0' != 0
      	  djn<-(dn/dn[1])*(-b1) # j_0'=-j_1
      }else{                     # j_0' == 0
      	  djn<-(dn/dn[2])*b0    # j_1' = (2/x)j_0'+j_0 = j_0 # APROXIMADO!
      }
      return(data.frame(jn,djn))
   }
}
