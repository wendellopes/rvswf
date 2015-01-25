#' Calculates the ratio of Cylindrical Bessel Functions.
#' 
#' @details Using Lentz method is possible to calculate the ratio 
#' \eqn{\gamma_n(x)=J_{n}(x)/J_{n+1}(x)}. By downward recurrence one can 
#' calculate \eqn{\gamma_n}.
#' @param n The order of \eqn{\gamma_n(x)}.
#' @param x The argument of \eqn{\gamma_n(x)} of type complex.
#' @param NMAX The maximum number of iterations.
#' @param code Choice between \code{C} or native {R} code.
#' @return The value of \eqn{\gamma_n} for complex arguments.
#' @seealso \code{\link{lcfa.cyl}}, \code{\link{cfe.cbd}}.
lcfe.cbd<-function(n,x,NMAX=2000,code="C"){
   n<-as.integer(n)
   nmaxo<-NMAX
   fn<-0
   if(!code%in%c("C","R")){
      stop("Code must be \"C\" or \"R\"")
   }
   if(code=="C"){
      if(is.complex(x)){
         u<-.C("lcfc_cbd",
               n=as.integer(n),
               x=as.complex(x),
               NMAX=as.integer(NMAX),
               fn=as.complex(fn)
         )
      }else{
         u<-.C("lcfe_cbd",
               n=as.integer(n),
               x=as.double(x),
               NMAX=as.integer(NMAX),
               fn=as.double(fn)
         )
      }

      if(u$NMAX>nmaxo){
         stop(paste("The accuracy was not reached in",NMAX,"iterations!"))
      }
      return(u$fn)
   }
   else{
      n<-as.integer(n)
      # Constants
      eo<-.Machine$double.xmin
      ACC<-10^-50
      # initialization of calculations
      fn<-lcfe.afs(2*(n+1),x)    # bo
      if(abs(fn)<eo){fn<-eo} # migth be zero
      Pn<-fn
      Qn<-0
      # Loop Parameters
      j<-0
      Dn<-10
      while(abs(Dn-1)>ACC){
         if(j>NMAX){
            NMAX<-j
            break
         }
         j<-j+1
         an<--1
         bn<-lcfe.afs(2*(n+j+1),x);
         Pn<-bn+an/Pn
         if(abs(Pn)<eo){Pn<-eo} # migth be zero
         Qn<-bn+an*Qn
         if(abs(Qn)<eo){Qn<-eo} # migth be zero
         Qn<-1/Qn
         Dn<-Pn*Qn
         fn<-fn*Dn
      }
      if(NMAX>nmaxo){
         stop(paste("The accuracy was not reached in",NMAX,"iterations!"))
      }
   return(fn)
   }
}
