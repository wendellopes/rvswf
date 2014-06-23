#-------------------------------------------------------------------------------
# RHO FOR SPHERICAL AND RICCATI BESSEL Y
#-------------------------------------------------------------------------------
comp.rsy<-function(n,x){
   a<-besselY(x,n+ .5)/besselY(x,n+1.5)
   b<-besselY(x,n+1.5)/besselY(x,n+2.5)
   c<-(2*(n+1)+1)/x
   d<--besselJ(x,-n- .5)/besselJ(x,-n-1.5)
   e<--besselJ(x,-n-1.5)/besselJ(x,-n-2.5)
   #------------------------------------
   cat("a=y[n]/y[n+1]\n")               # 1
   cat("b=y[n+1]/y[n+2]\n")             # 2
   cat("c=(2n+3)/x\n")                  # 3
   cat("d=j[-n-1]/j[-n-2]\n")           # 4
   cat("e=j[-n-2]/j[-n-3]\n")           # 5
   cat("u=a+1/b=c\n")                   # 6
   cat("v=-1/lcfe.sbd(-n-2,x)\n")       # 7
   cat("w=-lcfe.sbi(-n-2,x)\n")         # 8
   cat("x=-1/lcfe.sbd(-n-2,x,\"R\")\n") # 9
   cat("y=-lcfe.sbi(-n-2,x,\"R\")\n")   # 10
   #------------------------------------
   names=c("a","b","c","d","e",         # 1,2,3,4,5
           "a+1/b",                     # 6
           "C.sbd","C.sbi",             # 7,8
           "R.sbd","R.sbi")             # 9,10
   value=c(a,b,c,d,e,                   # 1,2,3,4,5
           a+1/b,                       # 6
           -1/lcfe.sbd(-n-2,x),         # 7
           -lcfe.sbi(-n-2,x),           # 8
           -1/lcfe.sbd(-n-2,x,code="R"),# 9
           -lcfe.sbi(-n-2,x,code="R"))  # 10
   return(cbind(names,value))
}