#' Array of logarithmic derivative and ratio of Ricatti-Bessel functions.
#' 
#' @details Calculate by downward recurrence the derivative of Ricatti-Bessel functions
#' \eqn{C_n=\psi_n'(x)/\psi_n(x)}, where \eqn{\psi_n(x)=xj_n(x)}, and simultaneously
#' the ratio \eqn{\rho_n=\psi_n(x)/\psi_{n+1}(x)=j_n(x)/j_{n+1}(x)}, from 0 to 
#' \code{nmax}. The starting values are calculated by Lentz continued fraction
#' method.
#' @param nmax Maximun value of \code{n}.
#' @param x The argument of the functions. Can be complex.
#' @param code If the \code{C} code or \code{R} built in.
#' @param NMAX Maximum number of iterations on the Lentz algorithm.
#' @seealso \code{\link{lcfe.rbl}}, \code{\link{lcfe.sbd}}, \code{\link{lcfe.sbi}}.
#' @examples
#' nmax<-10
#' x<-5
#' u.c<-lcfa.ric(nmax,x,code="C")
#' u.r<-lcfa.ric(nmax,x,code="R")
#' u<-data.frame(
#'    # Logarithmic Derivatives
#'    C.LogDev=u.c$Cn,R.LogDev=u.r$Cn,
#'    # Ratio between Ricatti-Bessel functions
#'    C.RicRat=u.c$rn,R.RicRat=u.r$rn)
#' print(u)
#-------------------------------------------------------------------------------
# RICCATI-BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
#---------------------------------------
# \rho_{n-1}+\frac{1}{\rho_n}=S_{2n+1}
# (n+1)\rho_{n-1}-n\frac{1}{rho_n}=(2n+1)C_{n}
#
# \rho_{n}=S_{n+1}+C_{n+1}
# \frac{1}{\rho_n}=S_{n+1}-C_n
#
# (S_{n+1}-C_n)(S_{n+1}+C_{n+1})=1
#---------------------------------------
lcfa.ric<-function(nmax,x,code="C",NMAX=200){ # PROBLEMAS COM ZEROS #
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      if(is.complex(x)){
         u<-.C("lcfc_ric",
               nmax=as.integer(nmax),
               x=as.complex(x),
               rn=as.complex(dummy),
               Cn=as.complex(dummy),
               NMAX=as.integer(NMAX))
         return(data.frame(rn=u$rn,Cn=u$Cn))
      }else{
         u<-.C("lcfa_ric",
               nmax=as.integer(nmax),
               x=as.double(x),
               rn=as.double(dummy),
               Cn=as.double(dummy),
               NMAX=as.integer(NMAX))
         return(data.frame(rn=u$rn,Cn=u$Cn))
      }
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   Cn<-rep(0,nmax+1)  # Vector 
	   rn<-rep(1,nmax+1)  # Vector
	   Cn[nmax+1]<-lcfe.rbl(nmax,x)  # Last element - c_n=Cn[n+1]
	   rn[nmax+1]<-lcfe.sbd(nmax,x)  # Last element - rho_n=rn[n+1]
	   # DOWNWARD RECURRENCE     
      for(n in nmax:1){             # position (nmax+1):2 -> element nmax:1
      	rn[n]<-S(n,x)+Cn[n+1]      # [OK]
      	Cn[n]<-S(n,x)-1/rn[n]      # [OK]
      }
      return(data.frame(rn,Cn))
   }
}
