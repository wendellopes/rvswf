#-------------------------------------------------------------------------------
# Derivative of Cylindrical Hankel
#-------------------------------------------------------------------------------
reff.cdh<-function(x,n,type=1){
   if(abs(type)!=1){
      stop("type must be plus or minus 1!")
   }
   return(.5*(reff.ch1(x,n-1,type)-reff.ch1(x,n+1,type)))
}
#-------------------------------------------------------------------------------
# Spherical Bessel Functions
#-------------------------------------------------------------------------------