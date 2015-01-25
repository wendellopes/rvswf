#' Calculates Cylindrical Bessel functions from 0 to nmax.
#'
#' @details \code{bess.cyl} calculates the Cylindrical Bessel
#' function using downward recurrence, from \eqn{J_nmax(x)} to \eqn{J_0(x)}.
#' The system of equations is given by \eqn{S_n(x)=n/x}, 
#' \eqn{\gamma_n=J_n(x)/J_{n+1}(x)}{g[n]=J_n/J_{n+1}} and 
#' \eqn{D_n=J_n'(x)/J_n(x)}. The system can be solved by means of
#' the recurrence relations of the Cylindrical Bessel functions
#' \deqn{\gamma_{n-1}+\frac{1}{\gamma_n}=2S_{n}}{g[n-1]+1/g[n]=2 S[n]}
#' \deqn{\gamma_{n-1}-\frac{1}{\gamma_n}=2D_{n}}{g[n-1]-1/g[n]=2 D[n]}
#' that can be rewriten
#' \deqn{\gamma_{n}=S_{n+1}+D_{n+1}            }{  g[n]=S[n+1]+D[n+1]}
#' \deqn{\frac{1}{\gamma_n}=S_n-D_n.           }{1/g[n]=S[n  ]-D[n  ].}
#' The logarithmic derivatives obeys the relation,
#' \deqn{(S_n-D_n)(S_{n+1}+D_{n+1})=1.         }{(S[n]-D[n])(S[n+1]+D[n+1])=1.}
#' The values can be calculated upward or downward.
#' @param nmax The maximum order of \eqn{J_n(x)}
#' @param x The argument of \eqn{J_n(x)}
#' @param code If you prefer to use native R or C language.
#' The algorithm is the same.
#' @return An array of Cylindrical Bessel functions and its derivatives.
#' from 0 to \code{nmax} at point \code{x}
#' @examples
#' x<-30
#' nmax<-50
#' a<-bess.cyl(nmax,x,code="C")
#' b<-bess.cyl(nmax,x,code="R")
#' d<-besselJ(x=x,nu=0:nmax)
#' plot(a$Jn,type='b')
#' points(b$Jn,col='red',pch=4)
#' points(d,col='blue',pch=3)
bess.cyl<-function(nmax,x,code="C"){ # PROBLEMAS COM ZEROS #
   if(abs(x)<1e-10){
      return(data.frame(Jn=c(1,rep(0,nmax)),dJn=rep(0,nmax+1)))
   }
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("bess_cyl",
            nmax=as.integer(nmax),
            x=as.double(x),
            Jn=as.double(dummy),
            Dn=as.double(dummy),
            NMAX=as.integer(2000))
      return(data.frame(Jn=u$Jn,dJn=u$Dn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   Dn<-rep(0,nmax+1)  # Vector 
	   gn<-rep(1,nmax+1)  # Vector
	   gm<-gn                        # gn modified: gm -> J_{nmax}=1 -> J_{nmax}'=gn
	   Dm<-Dn                        # gn modified: gm -> J_{nmax}=1 -> J_{nmax}'=gn
	   Dn[nmax+1]<-lcfe.cbl(nmax,x)  # Last element - D_n=Dn[n+1]
	   gn[nmax+1]<-lcfe.cbd(nmax,x)  # Last element - gamma_n=gn[n+1]
	   Dm<-Dn                        # Dn modified - Dm
	   # DOWNWARD RECURRENCE         # Many choices
         for(n in nmax:1){              # position (nmax+1):2 -> element nmax:1
      	   # gamma
#            gn[n]=S(2*n,x)-1/gn[n+1]           # [OK]
      	   gn[n]<-S(n,x)+Dn[n+1]
      	   Dn[n]<-S(n-1,x)-1/gn[n]              # [OK]
      	   # modified (permits one step normalization)
      	   gm[n]<-S(n,x)*gm[n+1]+Dm[n+1]
      	   Dm[n]<-S(n-1,x)*gm[n]-gm[n+1]
      	   # Normalization
      	   if(abs(gm[n])>1e2){
      	   	 #cat("renorming...\n")
      	   	 Dm<-Dm/gm[n] # this must be done first
      	   	 gm<-gm/gm[n] # otherwise the result will be wrong.
      	   }
      }
      # one step normalization taking care about zeros
      # Bessel function
      J0<-besselJ(x,0)
      J1<-besselJ(x,1)
      if(abs(J0>1e-10)){
         Jn<-(gm/gm[1])*J0 # create functions for normalizations
      }else{
         Jn<-(gm/gm[2])*J1   	
      }
      # Its Derivative
      if(abs(J1)>1e-10){
      	  dJn<-(Dm/Dm[1])*(-J1)
      }else{
      	  dJn<-(Dm/Dm[2])*.5*(Jn[1]-Jn[3])
      }
      return(data.frame(Jn,dJn))
   }
}
# no<-4
# xo<-5
# cat("Dn    =",lcfe.cbl(no-1,xo),"\n")
# cat("gamma =",lcfe.cbd(no-1,xo),"\n")
# print(bess.cyl(no,xo))
