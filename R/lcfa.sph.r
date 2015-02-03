#' Array of logarithmic derivative and ratio of Spherical Bessel functions.
#' 
#' @details Calculate by downward recurrence the derivative of Spherical Bessel functions
#' \eqn{c_n=j_n'(x)/j_n(x)}  and simultaneously
#' the ratio \eqn{\rho_n=j_n(x)/j_{n+1}(x)}, from 0 to 
#' \code{nmax}. The starting values are calculated by Lentz continued fraction
#' method.
#' @param nmax Maximun value of \code{n}.
#' @param x The argument of the functions. Can be complex.
#' @param code If the \code{C} code or \code{R} built in.
#' @param NMAX Maximum number of iterations on the Lentz algorithm.
#' @seealso \code{\link{lcfe.sbl}}, \code{\link{lcfe.sbd}}, \code{\link{lcfe.sbi}}.
#' @export
#' @examples
#' nmax<-10
#' x<-5
#' u.c<-lcfa.sph(nmax,x,code="C")
#' u.r<-lcfa.sph(nmax,x,code="R")
#' u<-data.frame(
#'    # Logarithmic Derivatives
#'    C.LogDev=u.c$cn,R.LogDev=u.r$cn,
#'    # Ratio between Spherical Bessel functions
#'    C.SphRat=u.c$rn,R.SphRat=u.r$rn)
#' print(u)
#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS [DONE]
#-------------------------------------------------------------------------------
#---------------------------------------
# \rho_{n-1}+\frac{1}{\rho_n}=S_{2n+1}
# n\rho_{n-1}-(n+1)\frac{rho_n}=(2n+1)c_{n}
#
# \rho_{n}=S_{n+2}+c_{n+1}
# \frac{1}{\rho_n}=S_n-c_n
#
# (S_n-c_n)(S_{n+2}+c_{n+1})=1
#---------------------------------------
lcfa.sph<-function(nmax,x,code="C",NMAX=200){ # PROBLEMAS COM ZEROS #
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      dummy<-rep(0,nmax+1)
      if(is.complex(x)){
         u<-.C("lcfc_sph",
               nmax=as.integer(nmax),
               x=as.complex(x),
               rn=as.complex(dummy),
               cn=as.complex(dummy),
               NMAX=as.integer(NMAX))
         return(data.frame(rn=u$rn,cn=u$cn))  
      }else{
         u<-.C("lcfa_sph",
               nmax=as.integer(nmax),
               x=as.double(x),
               rn=as.double(dummy),
               cn=as.double(dummy),
               NMAX=as.integer(NMAX))
         return(data.frame(rn=u$rn,cn=u$cn))
      }
   }else{
      S<-function(n,x){
         return(n/x)
      }
	   cn<-rep(0,nmax+1)  # Vector 
	   rn<-rep(1,nmax+1)  # Vector
	   cn[nmax+1]<-lcfe.sbl(nmax,x)  # Last element - c_n=cn[n+1]
	   rn[nmax+1]<-lcfe.sbd(nmax,x)  # Last element - rho_n=rn[n+1]
	   # DOWNWARD RECURRENCE     
      for(n in nmax:1){             # position (nmax+1):2 -> element nmax:1
      	rn[n]<-S(n+1,x)+cn[n+1]    # [OK]
      	cn[n]<-S(n-1,x)-1/rn[n]    # [OK]
      }
      return(data.frame(rn,cn))
   }
}
