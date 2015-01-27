#' Reference function for Spherical Bessel Functions.
#' 
#' @details It uses de built in Cylindrical Bessel function to calculate the
#' Spherical Bessel function, taking care about the zero.
reff.sjn<-function(x,n){
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
