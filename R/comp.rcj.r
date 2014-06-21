#-------------------------------------------------------------------------------
# RHO FOR CYLINDRICAL BESSEL J
#-------------------------------------------------------------------------------
comp.rcj<-function(n,x){
   a<-besselJ(x,n)/besselJ(x,n+1)
   b<-besselJ(x,n+1)/besselJ(x,n+2)
   c<-(2*(n+1))/x
   cat("gamma[n]=J[n]/J[n+1]\n")
   cat("a=gamma[n]\n")
   cat("b=gamma[n+1]\n")
   cat("c=(2n+2)/x\n")
   cat("u=a+1/b=c\n")
   return(data.frame(a,b,c,u=a+1/b))
}
