#' Calculates the logarithmic derivative of the Cylindrical Bessel functions.
#'
#' @details \code{lcfa.cyl} calculates the logarithmic derivative of 
#' Cylindrical Bessel functions \eqn{D_n(x)=J_n'(x)/J_n(x)}{D[n]=J_n'/J_n}
#' using downward recurrence.
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
#' The values are calculated by downward recurrence, and the inicial values
#' calculated by Lentz method.
#' @param nmax The maximum order of \eqn{J_n(x)}
#' @param x The argument of \eqn{J_n(x)}
#' @param code If you prefer to use native R or C language.
#' @param NMAX Maximum number of iterations
#' @return An array of the logarithmic derivative of Cylindrical Bessel 
#' functions and the ratio between two consecutive Cylindrical Bessel functions.
#' from 0 to \code{nmax} at point \code{x}
#' @seealso \code{\link{lcfe.cbi}}, \code{\link{lcfe.cbl}}, \code{\link{cfe.cbd}}.
#' @examples
#' x<-5
#' nmax<-50
#' a<-lcfa.cyl(nmax,x,code="C")
#' b<-lcfa.cyl(nmax,x,code="R")
#' # RATIO BETWEEN CONSECUTIVE CYLINDRICAL BESSEL FUNCTIONS g_n=J_n/J_{n+1}
#' plot(a$gn)
#' points(b$gn,col='red',pch=4)
#' # LOGARITMIC DERIVATIVES \eqn{D_n=J_n'/J_n}
#' plot(a$Dn)
#' points(b$Dn,col='red',pch=4)
lcfa.cyl<-function(nmax,x,code="C",NMAX=200){ # PROBLEMAS COM ZEROS #
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      u<-.C("lcfa_cyl",
            nmax=as.integer(nmax),
            x=as.double(x),
            gn=as.double(dummy),
            Dn=as.double(dummy),
            NMAX=as.integer(NMAX))
      return(data.frame(gn=u$gn,Dn=u$Dn))
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   Dn<-rep(0,nmax+1)  # Vector 
	   gn<-rep(1,nmax+1)  # Vector
	   Dn[nmax+1]<-lcfe.cbl(nmax,x)  # Last element - D_n=Dn[n+1]
	   gn[nmax+1]<-lcfe.cbd(nmax,x)  # Last element - gamma_n=gn[n+1]
	   # DOWNWARD RECURRENCE
      for(n in nmax:1){             # position (nmax+1):2 -> element nmax:1
      	gn[n]<-S(n,x)+Dn[n+1]      # [OK]
      	Dn[n]<-S(n-1,x)-1/gn[n]    # [OK]
      }
      return(data.frame(gn,Dn))
   }
}