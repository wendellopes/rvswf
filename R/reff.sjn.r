#' Reference function for Spherical Bessel Functions.
#' 
#' @details It uses de built in Cylindrical Bessel function to calculate the
#' Spherical Bessel function, taking care about the zero.
#' @param n The order of the function
#' @param x The argument of the function
#' @return The value of the spherical Bessel function \eqn{j_n(x)}
#' @export
reff.sjn<-function(n,x){
   if(x==0){
      if(n==0){
         u<-1
      }else{
         u<-0
      }
   }else{
	   u<-sqrt(pi/2)*besselJ(x,n+.5)/sqrt(x)
   }
   return(u)
}
