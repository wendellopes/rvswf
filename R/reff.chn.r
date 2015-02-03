#' Reference function Cylindrical Hankel function (1,2).
#' 
#' @details The Cylindrical Hankel 
#' function given by \eqn{h_n^{(1,2)}=j_n(x)\pm y_n(x)}.
#' @param x The argument of the function
#' @param n The order of the function
#' @param type 1 or 2.
#' @export 
reff.chn<-function(x,n,type=1){
   if(abs(type)!=1){
      stop("type must be plus or minus 1!")
   }
   if(type==1){
      return(besselJ(x,n)+1i*besselY(x,n))
   }
   if(type==2){
      return(besselJ(x,n)-1i*besselY(x,n)) 
   }
}
