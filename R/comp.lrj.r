#-------------------------------------------------------------------------------
# RHO FOR CYLINDRICAL BESSEL J
#-------------------------------------------------------------------------------
comp.lrj<-function(n,x){
   #------------------------------------
   # Riccati Bessel Function
   # (S_{n+1}-C_n)(S_{n+1}+C_{n+1})=1
   # S_n=n/x
   # D_n=(x j_n(x))'/(x j_n(x))
   # S_n=lcfe.afs(n,x)
   # C_n=lcfe.rbl(n,x)
   #------------------------------------
   a<-reff.rdj(x,n  )/reff.rjn(x,n  )
   b<-reff.rdj(x,n+1)/reff.rjn(x,n+1)
   c<-lcfe.rbl(n  ,x)
   d<-lcfe.rbl(n+1,x)
   e<-lcfe.afs(n+1,x)
   f<-lcfe.afs(n+1,x)
   g<-(n+1)/x
   h<-(n+1)/x
   #------------------------------------
   cat("a<-reff.rdj(x,n  )/reff.rjn(x,n  )\n")
   cat("b<-reff.rdj(x,n+1)/reff.rjn(x,n+1)\n")
   cat("c<-lcfe.rbl(n  ,x)                \n")
   cat("d<-lcfe.rbl(n+1,x)                \n")
   cat("e<-lcfe.afs(n  ,x)                \n")
   cat("f<-lcfe.afs(n+1,x)                \n")
   #------------------------------------
   names=c("C_{n  }",             # a,c
           "C_{n+1}",             # b,d
           "S_{n+1}",             # g,e
           "S_{n+1}",             # h,f
           "1")                   # 1
   v.ref=c(a,b,g,h,(g-a)*(h+b))   # 
   v.cal=c(c,d,e,f,(e-c)*(f+d))   # 
   return(cbind(names,v.ref,v.cal))
}
