#-------------------------------------------------------------------------------
# RHO FOR SPHERICAL BESSEL j
#-------------------------------------------------------------------------------
comp.lsj<-function(n,x){
   #------------------------------------
   # Spherical Bessel Function
   # (S_n-c_n)(S_{n+2}+c_{n+1})=1
   # S_n=n/x
   # c_n=j_n'(x)/j_n(x)
   # S_n=lcfe.afs(n,x)
   # c_n=lcfe.cbl(n,x)
   #------------------------------------
   a<-reff.sdj(x,n  )/reff.sjn(x,n  )
   b<-reff.sdj(x,n+1)/reff.sjn(x,n+1)
   c<-lcfe.sbl(n  ,x)
   d<-lcfe.sbl(n+1,x)
   e<-lcfe.afs(n  ,x)
   f<-lcfe.afs(n+2,x)
   g<-n/x
   h<-(n+2)/x
   #------------------------------------
   cat("a<-reff.sdj(x,n  )/reff.sjn(x,n  )\n")
   cat("b<-reff.sdj(x,n+1)/reff.sjn(x,n+1)\n")
   cat("c<-lcfe.cbl(n  ,x)                \n")
   cat("d<-lcfe.cbl(n+1,x)                \n")
   cat("e<-lcfe.afs(n  ,x)                \n")
   cat("f<-lcfe.afs(n+1,x)                \n")
   #------------------------------------
   names=c("c_{n  }",             # a,c
           "c_{n+1}",             # b,d
           "S_{n+1}",             # g,e
           "S_{n+2}",             # h,f
           "1")                   # 1
   v.ref=c(a,b,g,h,(g-a)*(h+b))   # 
   v.cal=c(c,d,e,f,(e-c)*(f+d))   # 
   return(cbind(names,v.ref,v.cal))
}
