#-------------------------------------------------------------------------------
# RHO FOR SPHERICAL AND RICCATI BESSEL J
#-------------------------------------------------------------------------------
#------------------------------------
# Cylindrical Bessel Function
# \rho_n+1/\rho_{n+1}=S_{2(n+1)+1}
# S_n=n/x
# \rho_n=j_n(x)/j_{n+1}(x)
# S_n=lcfe.afs(n,x)
# \rho_n=lcfe.sbd(n,x)
# 1/\rho_n=lcfe.sbi(n,x)
#------------------------------------
comp.rsj<-function(n,x){
   a<-besselJ(x,n+ .5)/besselJ(x,n+1.5)
   b<-besselJ(x,n+2.5)/besselJ(x,n+1.5)
   c<-(2*(n+1)+1)/x
   d<-lcfe.sbd(n,x)
   e<-1/lcfe.sbd(n+1,x)
   #f<-lcfe.afs(2*(n+1)+1,x)
   g<-1/lcfe.sbi(n,x)
   h<-lcfe.sbi(n+1,x)
   #------------------------------------
   #------------------------------------
   names<-c("  rho_{n  }",       # a,c
            "1/rho_{n+1}",       # b,d
            "S_{2n+3}   ")       # g,e
   v.ref<-c(a,b,c)
   v.dir<-c(d,e,d+e)
   v.inv<-c(g,h,g+h)
   return(cbind(names,v.ref,v.dir,v.inv))
   #------------------------------------  
}