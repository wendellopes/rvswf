#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <float.h>
//#-------------------------------------------------------------------------------
//# Calculations must be performed for psi(n,x) and xi(n,x), where
//# 1<=n<=N and x=x and x=Mx. Number of points: np=2*N*nx*ny*nz
//#-------------------------------------------------------------------------------
//# FRACAO CONTINUADA PELO METODO DE LENTZ
//# f=bo+(a1/b1+)(a2/b2+)(a3/b3+)...(an/bn+...)
//# input: vetores an e bn em que an=(a0,a1,...,aN) e bn=(b0,b1,...,bN)
//# a0 may be zero or any number
//#-------------------------------------------------------------------------------
//# NAMES
//# lcf : Lentz Continued Fraction Evaluation
//# cbld : Cylindrical Bessel Log Derivative
//# cbri : Cylindrical Bessel Ratio Inverse (Barnett)
//# cbrd : Cylindrical Bessel Ratio Direct
//# sbld : Spherical Bessel Log Derivative
//# sbri : Spherical Bessel Ratio Inverse (Barnett)
//# sbrd : Spherical Bessel Ratio Direct
//# rbld : Riccati-Bessel Log Derivative
//# rbri : Riccati-Bessel Ratio Inverse (Barnett)
//# rbrd : Riccati-Bessel Ratio Direct
//# afsn : Auxiliary function Sn
//#-------------------------------------------------------------------------------
//# AUXILIARY FUNCTION
//lcf.afsn<-function(n,x){
//   return(n/x)
//}
double lcf_afsn(int n,double x){
   double s=n/x;
   return(s);
}
//#-------------------------------------------------------------------------------
//# J_{n}/J_{n+1} [OK] DIRECT
//lcf.cbrd<-function(n,x,NMAX=2000){
//   # Constants
//   eo<-.Machine$double.xmin
//   ACC<-10^-50
//   # initialization of calculations
//   fn<-lcf.afsn(2*(n+1),x)    # bo
//   if(abs(fn)<eo){fn<-eo} # migth be zero
//   Pn<-fn
//   Qn<-0
//   # Loop Parameters
//   j<-0
//   Dn<-10
//   while(abs(Dn-1)>ACC){
//      j<-j+1
//      an<--1
//      bn<-lcf.afsn(2*(n+j+1),x);
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
double lcf_cbrd(int n,double x,int NMAX){
   const double eo = DBL_MIN;
   double ACC=10^-50;
   double fn=lcf_afsn(2*(n+1),x);
   if(abs(fn)<eo){fn<-eo;} // migth be zero
   double Pn=fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10;
   double an;
   double bn;
   while(abs(Dn-1)>ACC){
      j=j+1;
      an=-1;
      bn=lcf_afsn(2*(n+j+1),x);
      Pn=bn+an/Pn;
      if(abs(Pn)<eo){Pn<-eo;} // migth be zero
      Qn=bn+an*Qn;
      if(abs(Qn)<eo){Qn<-eo;} // migth be zero
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
//#-------------------------------------------------------------------------------
//# j_{n}/j_{n+1} [OK] DIRECT
//lcf.sbrd<-function(n,x,NMAX=2000){
//   # Constants
//   eo<-.Machine$double.xmin
//   ACC<-10^-50
//   # initialization of calculations
//   fn<-lcf.afsn(2*(n+1)+1,x)    # bo
//   if(abs(fn)<eo){fn<-eo} # migth be zero
//   Pn<-fn
//   Qn<-0
//   # Loop Parameters
//   j<-0
//   Dn<-10
//   while(abs(Dn-1)>ACC){
//      j<-j+1
//      an<--1
//      bn<-lcf.afsn(2*(n+j+1)+1,x); 
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
double lcf_sbrd(int n,double x,int NMAX){
   const double eo = DBL_MIN;
   double ACC=10^-50;
   double fn=lcf_afsn(2*(n+1)+1,x);
   if(abs(fn)<eo){fn<-eo;} // migth be zero
   double Pn=fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10;
   double an;
   double bn;   
   while(abs(Dn-1)>ACC){
      j=j+1;
      an=-1;
      bn=lcf_afsn(2*(n+j+1)+1,x);
      Pn=bn+an/Pn;
      if(abs(Pn)<eo){Pn=eo;} // migth be zero
      Qn=bn+an*Qn;
      if(abs(Qn)<eo){Qn=eo;} // migth be zero
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
//#-------------------------------------------------------------------------------
//# J_{n+1}/J_{n} [OK] BARNETT
//lcf.cbri<-function(n,x,NMAX=2000){
//   # Constants
//   eo<-.Machine$double.xmin
//   ACC<-10^-50
//   # initialization of calculations
//   fn<-0    # bo
//   if(abs(fn)<eo){fn<-eo} # migth be zero
//   Pn<-fn
//   Qn<-0
//   # Loop Parameters
//   j<-0
//   Dn<-10
//   while(abs(Dn-1)>ACC){
//      an<-(-1)^sign(j)
//      j<-j+1
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
double lcf_cbri(int n,double x,int NMAX){
   const double eo = DBL_MIN;
   double ACC=10^-50;
   double fn=0.0;
   if(abs(fn)<eo){fn<-eo;} // migth be zero
   double Pn=fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10;
   double an;
   double bn;
   while(abs(Dn-1)>ACC){
      an=(-1)^(j/abs(j));
      j=j+1;
      bn=lcf_afsn(2*(n+j),x);
      Pn=bn+an/Pn;
      if(abs(Pn)<eo){Pn=eo;} //# migth be zero
      Qn=bn+an*Qn;
      if(abs(Qn)<eo){Qn=eo;} //# migth be zero
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
//#-------------------------------------------------------------------------------
//# j_{n+1}/j_{n} [OK] BARNETT
//lcf.sbri<-function(n,x,NMAX=2000){
//   # Constants
//   eo<-.Machine$double.xmin
//   ACC<-10^-50
//   # initialization of calculations
//   fn<-0    # bo
//   if(abs(fn)<eo){fn<-eo} # migth be zero
//   Pn<-fn
//   Qn<-0
//   # Loop Parameters
//   j<-0
//   Dn<-10
//   while(abs(Dn-1)>ACC){
//      an<-(-1)^sign(j)
//      j<-j+1
//      bn<-lcf.afsn(2*(n+j)+1,x);
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
double lcf_sbri(int n,double x,int NMAX){
   const double eo = DBL_MIN;
   double ACC=10^-50;
   double fn=0.0;
   if(abs(fn)<eo){fn<-eo;} // migth be zero
   double Pn=fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10;
   double an;
   double bn;
   while(abs(Dn-1)>ACC){
      an=(-1)^(j/abs(j));
      j=j+1;
      bn=lcf_afsn(2*(n+j)+1,x);
      Pn=bn+an/Pn;
      if(abs(Pn)<eo){Pn=eo;} //# migth be zero
      Qn=bn+an*Qn;
      if(abs(Qn)<eo){Qn=eo;} //# migth be zero
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

//#-------------------------------------------------------------------------------
//# Spherical Bessel [OK]
//lcf.sbld<-function(n,x,NMAX=2000){
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
//      bn<-lcf.afsn(2*(n+j)+1,x);
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
double lcf_sbld(int n,double x,int NMAX){;
   double eo=DBL_MIN;
   double ACC=10^-50;
   double fn=lcf_afsn(n,x);    // bo;
   if(abs(fn)<eo){fn=eo;} // migth be zero;
   double Pn=fn;
   double Qn=0.0;
   // Loop Parameters;
   double j=0;
   double Dn=10;
   double an;
   double bn;
   while(abs(Dn-1)>ACC){;
      j=j+1;
      an=-1;
      bn=lcf_afsn(2*(n+j)+1,x);;
      Pn=bn+an/Pn;
      if(abs(Pn)<eo){Pn=eo;} // migth be zero;
      Qn=bn+an*Qn;
      if(abs(Qn)<eo){Qn=eo;} // migth be zero;
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
//#-------------------------------------------------------------------------------
//# Riccati-Bessel [OK]
//lcf.rbld<-function(n,x,NMAX=2000){
//   # Constants
//   eo<-.Machine$double.xmin
//   ACC<-10^-50
//   # initialization of calculations
//   fn<-lcf.afsn(n+1,x)    # bo
//   if(abs(fn)<eo){fn<-eo} # migth be zero
//   Pn<-fn
//   Qn<-0
//   # Loop Parameters
//   j<-0
//   Dn<-10
//   while(abs(Dn-1)>ACC){
//      j<-j+1
//      an<--1
//      bn<-lcf.afsn(2*(n+j)+1,x);
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
double lcf_rbld(int n,double x,int NMAX){
   double eo=DBL_MIN;;
   double ACC=10^-50;
   double fn=lcf_afsn(n+1,x);    // bo;
   if(abs(fn)<eo){fn=eo;} // migth be zero;
   double Pn=fn;
   double Qn=0.0;
   double an;
   double bn;
   // Loop Parameters;
   int j=0;
   double Dn=10.0;
   while(abs(Dn-1)>ACC){;
      j=j+1;
      an=-1;
      bn=lcf_afsn(2*(n+j)+1,x);
      Pn=bn+an/Pn;
      if(abs(Pn)<eo){Pn=eo;} // migth be zero;
      Qn=bn+an*Qn;
      if(abs(Qn)<eo){Qn=eo;} // migth be zero;
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
//#-------------------------------------------------------------------------------
//# MIE CODE
//# Folowing the notation of Jianqi Shen and Xiaoshu Cai
//# Progress in Electromagnetics Research Symposium 2005, Hangzhou, China
//# August 22-26
//#-------------------------------------------------------------------------------
//# MIE CODE
//# an=A_n(x)Ta(mx)
//# bn=A_n(x)Tb(mx)
//# A_n(x)=\psi_n(x)/\xi_n(x)
//# Ta(x)=(L_n(mx)/m-L_n(x))/(L_n(mx)/m-B_n(x))
//# Tb(x)=(mL_n(mx)-L_n(x))/(mL_n(mx)-B_n(x))
//# B_n(x)=\xi_n'(x)/\xi_n(x)
//# DOWNWARD RECURRENCE
//# L_n(x), L_n(mx)
//# L_{n-1}=S_n-1/(S_n+L_n)
//# UPWARD RECURRENCE
//# A_n(x)=A_{n-1}(x)(B_n(x)+S_n(x))/(L_n(x)+S_n(x))
//# B_n(x)=-S_n(x)+1/(S_n(x)-B_{n-1}(x))
//# STARTING VALUES
//# A_1=1/(1+i(cos(x)+x sin(x))/(sin(x)-x cos(x)))
//# B_0=-i
//#-------------------------------------------------------------------------------
//#-------------------------------------------------------------------------------
//# WAVELENGTH (MICROMETERS)
//l<-.532
//# INDEX OF REFRACTION POLIESTIRENE
//np<-1.5725+.0031080/l^2+.00034779/l^4
//# INDEX OF REFRACTION POLIESTIRENE
//nw<-1.3253+.0027553/l^2+.00003779/l^4
//m<-np/nw
//#-------------------------------------------------------------------------------
//# MIE COEFFICIENTS BY MEANS OF LOGARITHMIC DERIVATIVES
//#-------------------------------------------------------------------------------
//mie.abld<-function(m,x,NMAX=floor(x+7.5*x^(1/3))+2){
//   # STOPPING CRITERIUM
//   # DOWNWARD RECURRENCE
//   An.1<-rep(-17,NMAX)
//   An.m<-rep(-17,NMAX)
//   # CALCULATIONS
//   An.1[NMAX]<-lcf.rbld(NMAX,  x)
//   An.m[NMAX]<-lcf.rbld(NMAX,m*x)
//   for(n in NMAX:2){
//      An.1[n-1]<-lcf.afsn(n,  x)-1/(n/(  x)+An.1[n])
//      An.m[n-1]<-lcf.afsn(n,m*x)-1/(n/(m*x)+An.m[n])
//   }
//   # UPWARD RECURRENCE
//   Cn<-Bn<-rep(1,NMAX)
//   Cn[1]<-1/(1+1i*(cos(x)+x*sin(x))/(sin(x)-x*cos(x)))
//   Bn[1]=-lcf.afsn(1,x)+1/(lcf.afsn(1,x)+1i)
//   for(n in 2:NMAX){
//      Bn[n]<--lcf.afsn(n,x)+1/(lcf.afsn(n,x)-Bn[n-1])
//      Cn[n]<-Cn[n-1]*(Bn[n]+n/x)/(An.1[n]+n/x)
//   }
//   # OTHER ExPRESSIONS
//   Ta<-(An.m/m-An.1)/(An.m/m-Bn)
//   Tb<-(An.m*m-An.1)/(An.m*m-Bn)
//   # MIE COEFFICIENTS
//   an<-Cn*Ta
//   bn<-Cn*Tb
//   # SCATTERING 
//   u<-data.frame(Cn,Ta,Tb,an,bn)
//   return(u)
//}
//#-------------------------------------------------------------------------------
//# MIE COEFFICIENTS BY MEANS OF RATIO BETWEEN BESSEL FUNCTIONS
//#-------------------------------------------------------------------------------
//# TODO: How to calculate \gamma_n=\zeta_{n}/\zeta_{n+1}
//#       ** solution: upward recurrence 
//# TODO: How to calculate C_n=\psi_n/\zeta_n
//#       ** 
//mie.abrb<-function(m,x,NMAX=floor(x+7.5*x^(1/3))+2,DIRECT=TRUE){
//   # STOPPING CRITERIUM
//   # We need the zeroth therm at postion [1]
//   rho.1<-rho.m<-g<-Cn<-rep(-17,NMAX+1)
//   # Calculatin \gamma_0=g[1]
//      p0<-sin(x)
//      q0<-cos(x)
//      z0<-p0+1i*q0
//      p1<-p0/x-q0
//      q1<-q0/x+p0
//      z1<-p1+1i*q1
//   # Starting series
//   g[1]<-z0/z1
//   Cn[1]<-p0/z0
//   # STARTING VALUES 
//   if(DIRECT){
//      rho.1[NMAX+1]<-lcf.sbrd(NMAX,  x)
//      rho.m[NMAX+1]<-lcf.sbrd(NMAX,m*x)
//   }else{
//      rho.1[NMAX+1]<-1/lcf.sbri(NMAX,  x)
//      rho.m[NMAX+1]<-1/lcf.sbri(NMAX,m*x)
//   }
//   # DOWNWARD RECURRENCE
//   for(n in NMAX:1){
//      rho.1[n]<-lcf.afsn(2*n+1,  x)-1/rho.1[n+1]
//      rho.m[n]<-lcf.afsn(2*n+1,m*x)-1/rho.m[n+1]
//   }
//   # UPWARD RECURRENE 
//   for(n in 1:NMAX){
//      g[n+1]<-1/(lcf.afsn(2*n+1,x)-g[n])
//      Cn[n+1]<-Cn[n]*g[n]/rho.1[n]
//   }
//   # OTHER ExPRESSIONS
//   n<-1:(NMAX+1)
//   k<-(1-1/m^2)*n/x
//   Ta<-(rho.m/m-rho.1+k)/(rho.m/m-g+k)
//   Tb<-(rho.m*m-rho.1  )/(rho.m*m-g  )
//   # MIE COEFFICIENTS
//   n<-1:NMAX
//   Cn<-Cn[n+1]
//   Ta<-Ta[n]
//   Tb<-Tb[n]
//   an<-Cn*Ta
//   bn<-Cn*Tb
//   # SCATTERING 
//   u<-data.frame(Cn,Ta,Tb,an,bn)
//   #return(u[2:(NMAX+1),])
//   return(u)
//}
//#-------------------------------------------------------------------------------
//mie.anbn<-function(m,x,...){
//   return(mie.abld(m,x,...))
//}
//#-------------------------------------------------------------------------------
//# PARTICLE DEPENDENT FORCE CONSTANTS
//#-------------------------------------------------------------------------------
//ofc.abcn<-function(m,x){
//   k<-mie.anbn(m,x)
//   n0<-1:(nrow(k)-1)
//   n1<-n0+1
//   An<-k$an[n1]+Conj(k$an[n0])-2*k$an[n1]*Conj(k$an[n0])
//   Bn<-k$bn[n1]+Conj(k$bn[n0])-2*k$bn[n1]*Conj(k$bn[n0])
//   Cn<-k$an[n0]+Conj(k$bn[n0])-2*k$an[n1]*Conj(k$an[n0])
//   return(data.frame(An,Bn,Cn))
//}
//#-------------------------------------------------------------------------------
//# G^TE, G^TM
//#-------------------------------------------------------------------------------
//G<-function(m,x){
//   u<-ofc.abcn(m,x)
//   pM<-nrow(u)
//   An<-Bn<-Cn<-TE<-TM<-l<-m<-j<-1:(pM*(pM+2))
//   for(k in 1:pM){
//      s<--k:k
//      n<-k*(k+1)+s
//      l[n]<-k
//      m[n]<-s
//      j[n]<-n
//      TE[n]<-rnorm(n)
//      TM[n]<-rnorm(n)
//      An[n]<-u$An[k]
//      Bn[n]<-u$Bn[k]
//      Cn[n]<-u$Cn[k]
//   }
//   return(data.frame(j,l,m,An,Bn,Cn,TE,TM))
//}
//#
//# TESTING FUNCTIONS
//besselDJ<-function(x,n){
//   return(.5*(besselJ(x,n-1)-besselJ(x,n+1)))
//}
//besselH1<-function(x,n){
//   return(besselJ(x,n)+1i*besselY(x,n))
//}
//besselDH1<-function(x,n){
//   return(.5*(besselH1(x,n-1)-besselH1(x,n+1)))
//}
//besselH2<-function(x,n){
//   return(besselJ(x,n)-1i*besselY(x,n))
//}
//besselDH2<-function(x,n){
//   return(.5*(besselH2(x,n-1)-besselH2(x,n+1)))
//}
//mie.abcd<-function(n,x){
//   An<-.5/x+besselDJ(x,n+.5)/besselJ(x,n+.5)
//   Bn<-.5/x+besselDH2(x,n+.5)/besselH2(x,n+.5)
//   Cn<-besselJ(x,n+.5)/besselH2(x,n+.5)
//   return(data.frame(An,Bn,Cn))
//}
//mie.coef<-function(n,m,x){
//   u.1<-mie.abcd(n,x)
//   u.m<-mie.abcd(n,m*x)
//   Ta<-(u.m$An/m-u.1$An)/(u.m$An/m-u.1$Bn)
//   Tb<-(u.m$An*m-u.1$An)/(u.m$An*m-u.1$Bn)
//   an<-u.1$Cn*Ta
//   bn<-u.1$Cn*Tb
//   return(data.frame(An=u.1$An,Bn=u.1$Bn,Cn=u.1$Cn,Ta,Tb,an,bn))
//}
//#-------------------------------------------------------------------------------
//# SCATTERING COEFFICIENTS
//#-------------------------------------------------------------------------------
//mie.qnmx<-function(m,x){
//   ab<-mie.anbn(m,x)
//   n<-nrow(ab)
//   n<-1:n
//   an<-ab$an
//   bn<-ab$bn
//   Qsca<-(2/(x^2))*sum((2*n+1)*(abs(an)^2+abs(bn)^2))
//   Qext<-(2/(x^2))*sum((2*n+1)*Re(an+bn))
//   Q<-data.frame(Qsca,Qext)
//   return(Q)
//}
//mie.scat<-function(m,x){
//   u<-data.frame(Qsca=x,Qext=x)
//   for(i in 1:length(x)){
//      u[i,]<-mie.qnmx(m,x[i])
//   }
//   return(u)
//}
//#-------------------------------------------------------------------------------
//# TESTS
//#-------------------------------------------------------------------------------
//RCn<-function(n,x){
//   a<-besselJ(x,n+.5)
//   b<-besselY(x,n+.5)
//   return(1/(1+1i*b/a))
//}
//rhonxJ<-function(n,x){
//   a<-besselJ(x,n+ .5)/besselJ(x,n+1.5)
//   b<-besselJ(x,n+1.5)/besselJ(x,n+2.5)
//   c<-(2*n+3)/x
//   return(c(a,b,a+1/b,c))
//}
//rhonxY<-function(n,x){
//   a<-besselY(x,n+ .5)/besselY(x,n+1.5)
//   b<-besselY(x,n+1.5)/besselY(x,n+2.5)
//   c<-(2*n+3)/x
//   d<--1/rhonxJ(-n-2,x)[1]
//   return(c(a,b,a+1/b,c,d))
//}
//rhonxH<-function(n,x,k=1){
//   a<-(besselJ(x,n+ .5)+k*1i*besselY(x,n+ .5))/
//      (besselJ(x,n+1.5)+k*1i*besselY(x,n+1.5))
//   b<-(besselJ(x,n+1.5)+k*1i*besselY(x,n+1.5))/
//      (besselJ(x,n+2.5)+k*1i*besselY(x,n+2.5))
//   c<-(2*n+3)/x
//   return(c(a,b,a+1/b,c))
//}
