#' Reference function Derivative of Cylindrical Hankel function (1,2).
#' 
#' @details The derivative of the Cylindrical Hankel 
#' function given by \eqn{h_n^{(1,2)}=j_n(x)\pm y_n(x)}.
#' @param x The argument of the function
#' @param n The order of the function
#' @param type 1 or 2.
#' @include reff.chn.r
#' @export 
reff.cdh<-function(n,x,type=1){
   if(abs(type)!=1){
      stop("type must be plus or minus 1!")
   }
   return(.5*(reff.chn(x,n-1,type)-reff.chn(x,n+1,type)))
}
