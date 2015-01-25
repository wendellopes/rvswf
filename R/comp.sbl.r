#' Checks the results of Lentz method and internal R Bessel function
#' 
#' @details In order to check the results of Lentz method for calculation
#' of Cylindrical Bessel logarithmic derivative \eqn{c_n}.
#' @param n Order of the logarithmic derivative given by \eqn{c_n=j_n'/j_n}.
#' @param x Argument of Bessel functions.
#' @param code If C or native R function.
#' @return Data frame with the values calculated by the algorithm.
#' @seealso \code{\link{lcfe.sbi}}, \code{\link{lcfe.sbl}}, \code{\link{lcfe.sbd}},
#' \code{\link{lcfe.afs}}, \code{\link{besselJ}}.
#' @examples 
#' comp.sbl(5,4,code="C")
#' comp.sbl(5,4,code="R")
comp.sbl<-function(n,x,code="C"){
   a<-reff.sdj(x,n  )/reff.sjn(x,n  )
   b<-reff.sdj(x,n+1)/reff.sjn(x,n+1)
   s<-lcfe.sbl(n  ,x,code=code)
   d<-lcfe.sbl(n+1,x,code=code)
   e<-lcfe.afs(n  ,x)
   f<-lcfe.afs(n+2,x)
   g<-n/x
   h<-(n+2)/x
   #------------------------------------
   rnames=c("c_{n  }",             # a,s
            "c_{n+1}",             # b,d
            "S_{n  }",             # g,e
            "S_{n+1}",             # h,f
            "1")                   # 1
   m<-as.data.frame(array(0,c(5,2)))
   rownames(m)<-rnames
   names(m)<-c("Reference","Calculated")
   m$Reference <-c(a,b,g,h,(g-a)*(h+b))   # 
   m$Calculated<-c(s,d,e,f,(e-s)*(f+d))   # 
   return(m)
}
