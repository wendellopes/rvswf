#' Checks the results of Lentz method and internal R Bessel function
#' 
#' @details In order to check the results of Lentz method for calculation
#' of Cylindrical Bessel ratio of consecutive function \eqn{\rho_n}{r[n]}.
#' @param n Order of the ratio given by \eqn{\gamma_n=j_n/j_{n+1}}{r[n]=j_n/j_{n+1}}.
#' @param x Argument of Bessel functions.
#' @param code If C or native R function.
#' @return Data frame with the values calculated by the algorithms.
#' @seealso \code{\link{lcfe.sbi}}, \code{\link{lcfe.sbl}}, \code{\link{lcfe.sbd}},
#' \code{\link{lcfe.afs}}, \code{\link{besselJ}}.
#' @examples 
#' comp.sbd(5,4,code="C")
#' comp.sbd(5,4,code="R")
comp.sbd<-function(n,x,code="C"){
   a<-besselJ(x,n+ .5)/besselJ(x,n+1.5)
   b<-besselJ(x,n+2.5)/besselJ(x,n+1.5)
   s<-(2*(n+1)+1)/x
   d<-  lcfe.sbd(n  ,x,code=code)
   e<-1/lcfe.sbd(n+1,x,code=code)
   g<-1/lcfe.sbi(n  ,x,code=code)
   h<-  lcfe.sbi(n+1,x,code=code)
   #------------------------------------
   #------------------------------------
   rnames<-c("  rho_{n  }",       # a,s
             "1/rho_{n+1}",       # b,d
             "S_{2n+2}   ")       # g,e
   m<-as.data.frame(array(0,c(3,3)))
   rownames(m)<-rnames
   names(m)<-c("Internal","Direct","Inverse")
   m$Internal<-c(a,b,s)
   m$Direct  <-c(d,e,d+e)
   m$Inverse <-c(g,h,g+h)
   return(m)
   #------------------------------------ 
 
}
