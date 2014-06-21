tst.abcd<-function(n,x){
   An<-.5/x+besselDJ(x,n+.5)/besselJ(x,n+.5)
   Bn<-.5/x+besselDH2(x,n+.5)/besselH2(x,n+.5)
   Cn<-besselJ(x,n+.5)/besselH2(x,n+.5)
   return(data.frame(An,Bn,Cn))
}

