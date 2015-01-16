#' Checks the results of Lentz method and internal R Bessel function
#' 
#' @details In order to check the results of Lentz method for calculation
#' of Cylindrical Bessel ratio of consecutive function \eqn{\gamma_n}{g[n]}.
#' @param n Order of the ratio given by \eqn{\gamma_n=J_n/J_{n+1}}{g[n]=J_n/J_{n+1}}.
#' @param x Argument of Bessel functions.
#' @return Data frame with the values calculated by the algorithms.
#' @seealso \code{\link{lcfe.cbi}}, \code{\link{lcfe.cbl}}, \code{\link{lcfe.cbd}},
#' \code{\link{lcfe.afs}}, \code{\link{besselJ}}.
#' @examples 
#' comp.cbd(5,4)
comp.cbd<-function(n,x){
   #------------------------------------
   # Cylindrical Bessel Function
   # \gamma_n+1/\gamma_{n+1}=S_{2(n+1)}
   # S_n=n/x
   # \gamma_n=J_n(x)/J_{n+1}(x)
   # S_n=lcfe.afs(n,x)
   # \gamma_n=lcfe.cbd(n,x)
   # 1/\gamma_n=lcfe.cbi(n,x)
   #------------------------------------
   a<-besselJ(x,n)/besselJ(x,n+1)
   b<-besselJ(x,n+2)/besselJ(x,n+1)
   s<-(2*(n+1))/x
   d<-lcfe.cbd(n,x)
   e<-1/lcfe.cbd(n+1,x)
   #f<-lcfe.afs(2*(n+1),x)
   g<-1/lcfe.cbi(n,x)
   h<-lcfe.cbi(n+1,x)
   #------------------------------------
   #------------------------------------
   rnames<-c("  gamma_{n  }",       # a,c
             "1/gamma_{n+1}",       # b,d
             "S_{2n+2}     ")       # g,e
   m<-as.data.frame(array(0,c(3,3)))
   rownames(m)<-rnames
   names(m)<-c("Internal","Direct","Inverse")
   m$Internal<-c(a,b,s)
   m$Direct  <-c(d,e,d+e)
   m$Inverse <-c(g,h,g+h)
   return(m)
   #------------------------------------  
}
