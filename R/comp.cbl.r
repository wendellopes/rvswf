#' Checks the results of Lentz method and internal R Bessel function
#' 
#' @details In order to check the results of Lentz method for calculation
#' of Cylindrical Bessel logarithmic derivative \eqn{D_n}.
#' @param n Order of the logarithmic derivative given by \eqn{D_n=J_n'/J_n}.
#' @param x Argument of Bessel functions.
#' @param code If C or native R function.
#' @return Data frame with the values calculated by the algorithm.
#' @seealso \code{\link{lcfe.cbi}}, \code{\link{lcfe.cbl}},
#' \code{\link{lcfe.afs}}, \code{\link{besselJ}}.
#' @import lcfe.cbl,lcfe.afs,reff.cdj
#' @export
#' @examples 
#' comp.cbl(5,4,code="C")
#' comp.cbl(5,4,code="R")
comp.cbl<-function(n,x,code="C"){
   a<-reff.cdj(x,n  )/besselJ(x,n  )
   b<-reff.cdj(x,n+1)/besselJ(x,n+1)
   c<-lcfe.cbl(n  ,x,code=code)
   d<-lcfe.cbl(n+1,x,code=code)
   e<-lcfe.afs(n  ,x)
   f<-lcfe.afs(n+1,x)
   g<-n/x
   h<-(n+1)/x
   rnames=c("D_{n  }",             # a,c
           "D_{n+1}",             # b,d
           "S_{n  }",             # g,e
           "S_{n+1}",             # h,f
           "1")                   # 1
   m<-as.data.frame(array(0,c(5,2)))
   rownames(m)<-rnames
   names(m)<-c("Reference","Calculated")
   m$Reference <-c(a,b,g,h,(g-a)*(h+b))   # 
   m$Calculated<-c(c,d,e,f,(e-c)*(f+d))   # 
   return(m)
}
