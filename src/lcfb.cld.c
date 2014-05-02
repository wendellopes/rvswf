#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <float.h>
//#-------------------------------------------------------------------------------
//# Cylindrical Bessel [OK]
//lcf.cbld<-function(n,x,NMAX=2000){
//   # Constants
//   eo<-.Machine$double.xmin
//   ACC<-10^-50
//   # initialization of calculations
//   fn<-lcf.afsn(n,x)    # bo
//   if(abs(fn)<eo){fn<-eo} # migth be zero
//   Pn<-fn
//   Qn<-0
//   # Loop Parameters
//   j<-0
//   Dn<-10
//   while(abs(Dn-1)>ACC){
//      j<-j+1
//      an<--1
//      bn<-lcf.afsn(2*(n+j),x);
//      Pn<-bn+an/Pn
//      if(abs(Pn)<eo){Pn<-eo} # migth be zero
//      Qn<-bn+an*Qn
//      if(abs(Qn)<eo){Qn<-eo} # migth be zero
//      Qn<-1/Qn
//      Dn<-Pn*Qn
//      fn<-fn*Dn
//      if(j==NMAX){
//         cat("LIMITE DE",NMAX,"ITERACOES EXCEDIDO COM ACC =",abs(Dn-1),"\n")
//         break
//      }
//   }
//   return(fn)
//}
double lcf_cbld(int n,double x,int NMAX){
   double eo = DBL_MIN;
   double ACC=10^-50;
   double fn=lcf_afsn(n,x);
   if(abs(fn)<eo){fn<-eo;}
   double Pn=fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;
   while(abs(Dn-1)>ACC){
      j=j+1;
      an=-1;
      bn=lcf_afsn(2*(n+j),x);
      Pn=bn+an/Pn;
      if(abs(Pn)<eo){Pn=eo;}// # migth be zero
      Qn=bn+an*Qn;
      if(abs(Qn)<eo){Qn=eo;}// # migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      fn=fn*Dn;
      if(j==NMAX){
         printf("LIMITE DE %d ITERACOES EXCEDIDO COM ACC = %f \n",NMAX,abs(Dn-1));
         break;
      }
   }
   return(fn);
}
