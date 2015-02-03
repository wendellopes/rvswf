#' Cylindrical Bessel function of order one \eqn{J_1(x)}.
#' 
#' @details Cylindrical Bessel function calculated by series summation.
#' @export
bess.cun<-function(x){
   Kj<-function(x,j){
      return(x/j)
   }
   if(x>1){
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