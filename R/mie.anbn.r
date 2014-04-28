#-------------------------------------------------------------------------------
mie.anbn<-function(m,x,by="LD",...){
   if(!by%in%c("LD","RB")){
      stop("by must be Logorithm Derivative (LD) or 
             Ratio between Bessel functions (RB)")
   }
   if(by=="LD"){
   	  return(mie.abld(m,x,...))
   }else{
      return(mie.abrb(m,x,...))
   }
}
