#-------------------------------------------------------------------------------
# GAMMA FOR CYLINDRICAL BESSEL J
#-------------------------------------------------------------------------------
comp.rcj<-function(n,x){
   #------------------------------------
   # Cylindrical Bessel Function
   # \gamma_n+1/\gamma_{n+1}=S_{2(n+1)}
   # S_n=n/x
   # \gamma_n=J_n(x)/J_{n+1}(x)
   # S_n=lcfe.afs(n,x)
   # \gamma_n=lcfe.cbd(n,x)
   # 1/\gamma_n=lcfe.cbi(n,x)
   #------------------------------------
   a<-besselJ(x,n)/besselJ(x,n+1)
   b<-besselJ(x,n+2)/besselJ(x,n+1)
   c<-(2*(n+1))/x
   d<-lcfe.cbd(n,x)
   e<-1/lcfe.cbd(n+1,x)
   #f<-lcfe.afs(2*(n+1),x)
   g<-1/lcfe.cbi(n,x)
   h<-lcfe.cbi(n+1,x)
   #------------------------------------
   #------------------------------------
   names<-c("  gamma_{n  }",       # a,c
            "1/gamma_{n+1}",       # b,d
            "S_{2n+2}     ")       # g,e
   v.ref<-c(a,b,c)
   v.dir<-c(d,e,d+e)
   v.inv<-c(g,h,g+h)
   return(cbind(names,v.ref,v.dir,v.inv))
   #------------------------------------  
}
