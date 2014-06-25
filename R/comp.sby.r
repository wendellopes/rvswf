#-------------------------------------------------------------------------------
# RHO FOR SPHERICAL AND RICCATI BESSEL Y
#-------------------------------------------------------------------------------
comp.sby<-function(n,x){
   #------------------------------------
   # Cylindrical Bessel Function y
   # \rho*_n+1/\rho*_{n+1}=S_{2(n+1)+1}
   # \rho_{n}*=-1/rho_{-n-2}
   # S_n=n/x
   # \rho*_n=y_n/y_{n+1}=-j_{-n-1}/j_{-n-2}
   # S_n=lcfe.afs(n,x)
   # \rho_n=lcfe.sbd(n,x)
   # 1/\rho_n=lcfe.sbi(n,x)
   #------------------------------------
   a<-besselY(x,n+ .5)/besselY(x,n+1.5)
   b<-besselY(x,n+2.5)/besselY(x,n+1.5)
   c<-(2*(n+1)+1)/x
   #
   d<--besselJ(x,-n- .5)/besselJ(x,-n-1.5)
   e<--besselJ(x,-n-2.5)/besselJ(x,-n-1.5)
   f<-lcfe.afs(2*(n+1)+1,x)
   #
   g<--1/lcfe.sbd(-n-2,x)
   h<--lcfe.sbd(-n-3,x)
   #
   i<--lcfe.sbi(-n-2,x)
   j<--1/lcfe.sbi(-n-3,x)
   #------------------------------------
   #------------------------------------
   names<-c("  rho*_{n  }",       # a,c
            "1/rho*_{n+1}",       # b,d
            "S_{2n+3}    ")       # g,e
   v.rfy<-c(a,b,c)
   v.rfj<-c(d,e,f)
   v.dir<-c(g,h,g+h)
   v.inv<-c(i,j,i+j)
   return(cbind(names,v.rfy,v.rfj,v.dir,v.inv))
   #------------------------------------ 
}
