#-------------------------------------------------------------------------------
# RHO FOR SPHERICAL AND RICCATI BESSEL H
#-------------------------------------------------------------------------------
comp.rsh<-function(n,x,k=1){
   if(abs(k)!=1){
      stop("k must be plus or minus 1!")
   }
   a<-(besselJ(x,n+ .5)+k*1i*besselY(x,n+ .5))/
      (besselJ(x,n+1.5)+k*1i*besselY(x,n+1.5))
   b<-(besselJ(x,n+1.5)+k*1i*besselY(x,n+1.5))/
      (besselJ(x,n+2.5)+k*1i*besselY(x,n+2.5))
   c<-(2*(n+1)+1)/x
   cat("rho[n]=h[n]/h[n+1]\n")
   cat("a=rho[n]\n")
   cat("b=rho[n+1]\n")
   cat("c=(2n+3)/x\n")
   cat("u=a+1/b=c\n")
   return(data.frame(a,b,c,u=a+1/b))
}
