#' Spherical Bessel function of order one \eqn{j_1(x)}
#' 
#' @details It uses the series solution to calculate the values of 
#' the Spherical Bessel function for \eqn{x \to 0}.
#' @param x Complex argument of the Spherical Bessel function of order one.
#' @return Complex values for \eqn{j_1(x)}.
#' @examples
#' print(bess.sun(3+2i))
#' print(bess.sun(0:10))
bess.sun<-function(x){
   Kj<-function(x,j){
      return(x/j)
   }
   if(length(x)>1){
      return(sapply(x,bess.sun))
   }else{
      if(abs(x)>1){
         return((sin(x)/x-cos(x))/x)
      }else{
         S<-0
         for(j in seq(10,1,-1)){
            S<-(2*j)/(2*j+1)-Kj(x,2*j+1)*Kj(x,2*j+2)*S
         }
         S<-x*S/2
         return(S)
      }
   }
}