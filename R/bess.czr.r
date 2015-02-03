#' Cylindrical Bessel function of order zero \eqn{J_0(x)}
#' 
#' @param x Argument of the Spherical Bessel function of order zero.
#' @export
bess.czr<-function(x){
   Kj<-function(x,j){
      return(x/j)
   }
   if(x>1){
      return(sin(x)/x)
   }else{
      S<-1
      for(j in seq(10,1,-1)){
         S<-1-Kj(x,2*j)*Kj(x,2*j+1)*S
      }
      return(S)
   }
}