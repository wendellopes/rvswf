#-------------------------------------------------------------------------------
# RHO FOR SPHERICAL AND RICCATI BESSEL J
#-------------------------------------------------------------------------------
comp.rsj<-function(n,x){
   a<-besselJ(x,n+ .5)/besselJ(x,n+1.5)
   b<-besselJ(x,n+1.5)/besselJ(x,n+2.5)
   c<-(2*(n+1)+1)/x
   cat("rho[n]=j[n]/j[n+1]\n")
   cat("a=rho[n]\n")
   cat("b=rho[n+1]\n")
   cat("c=(2n+3)/x\n")
   cat("u=a+1/b=c\n")
   return(data.frame(a,b,c,u=a+1/b))
}
