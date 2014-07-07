#-------------------------------------------------------------------------------
# Lorentz-Mie Expansion
#-------------------------------------------------------------------------------
lmie.exp<-function(m,x,by="LD",...){
   if(!by%in%c("LD","RB")){
      stop("by must be Logorithm Derivative (LD) or 
             Ratio between Bessel functions (RB)")
   }
   if(by=="LD"){
   	  return(lmie.log(m,x,...))
   }else{
      return(lmie.rho(m,x,...))
   }
}
