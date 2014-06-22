#-------------------------------------------------------------------------------
# RHO FOR CYLINDRICAL BESSEL J
#-------------------------------------------------------------------------------
comp.rcj<-function(n,x){
   a<-besselJ(x,n)/besselJ(x,n+1)
   b<-besselJ(x,n+1)/besselJ(x,n+2)
   c<-(2*(n+1))/x
   #------------------------------------
   cat("a=gamma[n]=J[n]/J[n+1]\n")  # 1
   cat("b=gamma[n+1]\n")            # 2
   cat("c=(2n+2)/x\n")              # 3
   cat("u=a+1/b=c\n")               # 4
   cat("v=lcfe.cbd(n,x)\n")         # 5
   cat("w=1/lcfe.cbi(n,x)\n")       # 6
   cat("x=lcfe.cbd(n,x,\"R\")\n")   # 7
   cat("y=1/lcfe.cbi(n,x,\"R\")\n") # 8
   #------------------------------------
   names=c("a","b","c",             # 1,2,3
           "a+1/b",                 # 4
           "C.cbd","C.cbi",         # 5,6
           "R.cbd","R.cbi")         # 7,8
   value=c(a,b,c,                   # 1,2,3
           a+1/b,                   # 4
           lcfe.cbd(n,x),           # 5
           1/lcfe.cbi(n,x),         # 6
           lcfe.cbd(n,x,code="R"),  # 7
           1/lcfe.cbi(n,x,code="R"))# 8
   return(cbind(names,value))
}
