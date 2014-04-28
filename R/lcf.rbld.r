#-------------------------------------------------------------------------------
# Logarithmic Derivative Riccati-Bessel [OK]
lcf.rbld<-function(n,x,NMAX=2000){
   # Constants
   eo<-.Machine$double.xmin
   ACC<-10^-50
   # initialization of calculations
   fn<-lcf.afsn(n+1,x)    # bo
   if(abs(fn)<eo){fn<-eo} # migth be zero
   Pn<-fn
   Qn<-0
   # Loop Parameters
   j<-0
   Dn<-10
   while(abs(Dn-1)>ACC){
      j<-j+1
      an<--1
      bn<-lcf.afsn(2*(n+j)+1,x);
      Pn<-bn+an/Pn
      if(abs(Pn)<eo){Pn<-eo} # migth be zero
      Qn<-bn+an*Qn
      if(abs(Qn)<eo){Qn<-eo} # migth be zero
      Qn<-1/Qn
      Dn<-Pn*Qn
      fn<-fn*Dn
      if(j==NMAX){
         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
         break
      }
   }
   return(fn)
}
