#' Calculation of Spherical Bessel Functions.
#' 
#' @details We use Lentz method for calculating the starting value
#'  j_lmax(x)=bo+(a1/b1+)(a2/b2+)(a3/b3+)...(an/bn),
#'  and then we use downward recurrence to calculate the other values,
#'  from lmax-1 to 0.
#'  @param lmax The maximum value of \eqn{l}.
#'  @param xo The argument of \eqn{j_l}. Can be real or complex.
#'  @param compare Make the comparison with the GSL (Gnu Scientific Library).
#'  @param verbose Shows the normalizations on downward recurrence.
#'  @return An array with the Spherical Bessel functions and its derivatives.
#'  @seealso \code{\link{bess.sph}}, \code{\link{bess.cyl}}, \code{\link{bess.ric}},
#'  \code{\link{bess.czr}}, \code{\link{bess.cun}}, \code{\link{bess.szr}},
#'  \code{\link{bess.sun}}.
#'  @examples
#'  vswf.sbf(5,3,compare=TRUE)
#-------------------------------------------------------------------------------
vswf.sbf<-function(lmax,xo,compare=FALSE,verbose=FALSE){
   lmax.length<-lmax+1
   if(lmax<2*abs(xo)){
      lmax<-as.integer(2*abs(xo))
   }
   lmax<-lmax+10
   Tk<-function(x,k){return((2*k+1)/x)}
   eo<-.Machine$double.xmin
   ACC<-10^-50
#-------------------------------------------------------------------------------
# L > x
# Case x>L: Calculate j_l for L2=int(1.5x)
#-------------------------------------------------------------------------------
#   print(lmax/xo)
#-------------------------------------------------------------------------------
   fn<-lmax/xo
   if(abs(fn)==0){fn<-eo}
   Cn<-fn
   Dn<-0
   N<-1
   DN<-10
   fna<-c(fn)
   while(abs(DN-1)>ACC){
      an<--1
      bn<-Tk(xo,N)
      Cn<-bn+an/Cn
      if(abs(Cn)==0){Cn<-eo}
      Dn<-bn+an*Dn
      if(abs(Dn)==0){Dn<-eo}
      Dn<-1/Dn
      DN<-Cn*Dn
      fn<-fn*DN
      N<-N+1
      fna<-c(fn,fna)
      if(N>2000){
         break
      }
   }
#   print(N)
#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS CALCULATION
# DOWNWARD RECURRENCE 
# VALID FOR lmax>x
   f<-fn
   gn<-c(1)
   dgn<-c(f)
   gno<-1
   dgno<-f  # dgno=gno*dgno
   RN<-1
   for(n in lmax:1){
      gnom<-gno*(n+1)/xo+dgno
      dgno<-gnom*(n-1)/xo-gno
      gno<-gnom
      gn<-c(gno,gn)
      dgn<-c(dgno,dgn)
      if(abs(gno)>1e100){
         if(verbose){
            cat("RENORMALIZACAO ",RN,"N =",n,"\n")
         }
         RN<-RN+1
         gn<-gn/gno
         dgn<-dgn/gno
         dgno<-dgno/gno
         gno<-1
      }
   }
#-------------------------------------------------------------------------------
# NORMALIZATION
   xu<-sin(xo)/xo
   # Verify if xo is not a zero of j_0(x)
   if(abs(xu)>1e-5){                   # this way, normalize by j_0
      gn<-(sin(xo)/xo)*gn/gn[1]
   }else{                         # otherwise normalize by j_1
      gn<-(-cos(xo)+sin(xo)/xo)*(1/xo)*gn/gn[2]
   }
   if(abs(gn[2])>1e-5){
      dgn<-dgn*(-gn[2])/dgn[1]
   }else{
      dgn<-((1-(2/(xo^2)))*(sin(xo)/xo)+2*cos(xo)/(xo^2))*dgn/dgn[2]
   }
   if(compare & !is.complex(xo)){ #perform comparison with gnu scientif library
      library(gsl)
      gsl<-bessel_jl_steed_array(lmax=lmax,x=xo)
      u<-data.frame(jn_gsl=gsl,jn=gn,djn=dgn)
      return(u[1:lmax.length,])
   }else{
      u<-data.frame(jn=gn,djn=dgn)
      return(u[1:lmax.length,])
   }
}
