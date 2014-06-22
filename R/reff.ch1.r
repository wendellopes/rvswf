#-------------------------------------------------------------------------------
# Cylindrical Hankel
#-------------------------------------------------------------------------------
reff.ch1<-function(x,n,type=1){
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
