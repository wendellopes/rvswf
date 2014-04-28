#-------------------------------------------------------------------------------
# Calculation of Spherical Bessel Functions by Continued Fractions
# Lentz Method
# fn=bo+(a1/b1+)(a2/b2+)(a3/b3+)...(an/bn)
# DOWNWARD RECURRENCE WITH EXTRA RENORMALIZATIONS
# Works for Lmax>x. If Lmax<x, performs calculations for Lmax'>x and then
# discards terms greater than Lmax. (Not the best performance choice)
#-------------------------------------------------------------------------------
vswf.sbf<-function(NN,xo,compare=FALSE,verbose=FALSE){
   NN.length<-NN+1
   if(NN<2*xo){
      NN<-as.integer(2*xo)
   }
   NN<-NN+10
   Tk<-function(x,k){return((2*k+1)/x)}
   eo<-.Machine$double.xmin
   ACC<-10^-50
#-------------------------------------------------------------------------------
# L > x
# Case x>L: Calculate j_l for L2=int(1.5x)
#-------------------------------------------------------------------------------
#   print(NN/xo)
#-------------------------------------------------------------------------------
   fn<-NN/xo
   if(fn==0){fn<-eo}
   Cn<-fn
   Dn<-0
   N<-1
   DN<-10
   fna<-c(fn)
   while(abs(DN-1)>ACC){
      an<--1
      bn<-Tk(xo,N)
      Cn<-bn+an/Cn
      if(Cn==0){Cn<-eo}
      Dn<-bn+an*Dn
      if(Dn==0){Dn<-eo}
      Dn<-1/Dn
      DN<-Cn*Dn
      fn<-fn*DN
      N<-N+1
      fna<-c(fn,fna)
      if(N>2000){
         break
      }
   }
#   print(N)
#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS CALCULATION
# DOWNWARD RECURRENCE 
# VALID FOR lmax>x
   f<-fn
   gn<-c(1)
   dgn<-c(f)
   gno<-1
   dgno<-f  # dgno=gno*dgno
   RN<-1
   for(n in NN:1){
      gnom<-gno*(n+1)/xo+dgno
      dgno<-gnom*(n-1)/xo-gno
      gno<-gnom
      gn<-c(gno,gn)
      dgn<-c(dgno,dgn)
      if(abs(gno)>1e100){
         if(verbose){
            cat("RENORMALIZACAO ",RN,"N =",n,"\n")
         }
         RN<-RN+1
         gn<-gn/gno
         dgn<-dgn/gno
         dgno<-dgno/gno
         gno<-1
      }
   }
#-------------------------------------------------------------------------------
# NORMALIZATION
   xu<-(xo/pi)-as.integer(xo/pi)
   # Verify if xo is not a zero of j_0(x)
   if(xu>1e-3){                   # this way, normalize by j_0
      gn<-(sin(xo)/xo)*gn/gn[1]
   }else{                         # otherwise normalize by j_1
      cat("Zero de j_l(x)\n")
      gn<-(-cos(xo)+sin(xo)/xo)*(1/xo)*gn/gn[2]
   }
   # Verify if xo is not a zero of j_0'(x)=-j_1(x)
   # x=4.4934 (example)
   if(gn[2]>1e-2){
      dgn<-dgn*(-gn[2])/dgn[1]
   }else{
      cat("Zero de j_l'(x)\n")
      dgn<-((1-(2/(xo^2)))*(sin(xo)/xo)+2*cos(xo)/(xo^2))*dgn/dgn[2]
   }
   if(compare){ #perform comparison with gnu scientif library
      library(gsl)
      gsl<-bessel_jl_steed_array(lmax=NN,x=xo)
      u<-data.frame(jn_gsl=gsl,jn=gn,djn=dgn)
      return(u[1:NN.length,])
   }else{
      u<-data.frame(jn=gn,djn=dgn)
      return(u[1:NN.length,])
   }
}
