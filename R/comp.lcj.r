#-------------------------------------------------------------------------------
# RHO FOR CYLINDRICAL BESSEL J
#-------------------------------------------------------------------------------
comp.lcj<-function(n,x){
   #------------------------------------
   # Cylindrical Bessel Function
   # (S_n-D_n)(S_{n+1}+D_{n+1})=1
   # S_n=n/x
   # D_n=J_n'(x)/J_n(x)
   # S_n=lcfe.afs(n,x)
   # D_n=lcfe.cbl(n,x)
   #------------------------------------
   a<-test.cdj(x,n  )/besselJ(x,n  )
   b<-test.cdj(x,n+1)/besselJ(x,n+1)
   c<-lcfe.cbl(n  ,x)
   d<-lcfe.cbl(n+1,x)
   e<-lcfe.afs(n  ,x)
   f<-lcfe.afs(n+1,x)
   g<-n/x
   h<-(n+1)/x
   #------------------------------------
   cat("a<-test.cdj(x,n  )/besselJ(x,n  )\n")
   cat("b<-test.cdj(x,n+1)/besselJ(x,n+1)\n")
   cat("c<-lcfe.cbl(n  ,x)               \n")
   cat("d<-lcfe.cbl(n+1,x)               \n")
   cat("e<-lcfe.afs(n  ,x)               \n")
   cat("f<-lcfe.afs(n+1,x)               \n")
   #------------------------------------
   names=c("D_{n  }",             # a,c
           "D_{n+1}",             # b,d
           "S_{n  }",             # g,e
           "S_{n+1}",             # h,f
           "1")                   # 1
   v.ref=c(a,b,g,h,(g-a)*(h+b))   # 
   v.cal=c(c,d,e,f,(e-c)*(f+d))   # 
   return(cbind(names,v.ref,v.cal))
}
