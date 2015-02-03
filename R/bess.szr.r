#' Spherical Bessel function of order one \eqn{j_0(x)}
#' 
#' @details It uses the series solution to calculate the values of 
#' the Spherical Bessel function for \eqn{x \to 0}.
#' @param x Complex argument of the Spherical Bessel function of order one.
#' @return Complex values for \eqn{j_0(x)}.
#' @export
#' @examples
#' print(bess.szr(3+2i))
#' print(bess.szr(0:10))
bess.szr<-function(x){
   Kj<-function(x,j){
      return(x/j)
   }
   if(length(x)>1){
      return(sapply(x,bess.szr))
   }else{
      if(abs(x)>1){
         return(sin(x)/x)
      }else{
         S<-1
         for(j in seq(10,1,-1)){
            S<-1-Kj(x,2*j)*Kj(x,2*j+1)*S
         }
         return(S)
      }
   }
}