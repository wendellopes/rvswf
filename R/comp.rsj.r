#-------------------------------------------------------------------------------
# RHO FOR SPHERICAL AND RICCATI BESSEL J
#-------------------------------------------------------------------------------
comp.rsj<-function(n,x){
   a<-besselJ(x,n+ .5)/besselJ(x,n+1.5)
   b<-besselJ(x,n+1.5)/besselJ(x,n+2.5)
   c<-(2*(n+1)+1)/x
   #------------------------------------
   cat("a=rho[n]=j[n]/j[n+1]\n")    # 1
   cat("b=rho[n+1]\n")              # 2
   cat("c=(2n+3)/x\n")              # 3
   cat("u=a+1/b=c\n")               # 4
   cat("v=lcfe.sbd(n,x)\n")         # 5
   cat("w=1/lcfe.sbi(n,x)\n")       # 6
   cat("x=lcfe.sbd(n,x,\"R\")\n")   # 7
   cat("y=1/lcfe.sbi(n,x,\"R\")\n") # 8
   #------------------------------------
   names=c("a","b","c",             # 1,2,3
           "a+1/b",                 # 4
           "C.sbd","C.sbi",         # 5,6
           "R.sbd","R.sbi")         # 7,8
   value=c(a,b,c,                   # 1,2,3
           a+1/b,                   # 4
           lcfe.sbd(n,x),           # 5
           1/lcfe.sbi(n,x),         # 6
           lcfe.sbd(n,x,code="R"),  # 7
           1/lcfe.sbi(n,x,code="R"))# 8
   return(cbind(names,value))
}
