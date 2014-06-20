#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <float.h>
#include <gsl/gsl_sf_bessel.h>
#include <omp.h>
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
/*------------------------------------------------------------------------------
 * VSWF: Vector Spherical Wave Functions dynamic libraries written in C.
 * author: Wendel Lopes Moreira <wendellopes@gmail.com>
 * date: 2013-02-18
 * version 1.1
 * depends: Gnu Scientific Library <http://www.gnu.org/software/gsl/>
------------------------------------------------------------------------------*/
//------------------------------------------------------------------------------
// Psi_m(\vec{k},\vec{r}): Basic function for cylindrical simetries
// used in Cylindrical Wave Guides and Bessel Beams.
//------------------------------------------------------------------------------
double complex psi_mkr(int *m,int *s,
      double *gamma,double *kz,
      double *x,double *y,double *z){
   double rho=sqrt(*x*(*x)+*y*(*y));
   double cph,sph;
   cph=*x/rho;
   sph=*y/rho;
   double complex eiph=cpow(cph+I*sph,*m*(*s));
   double complex eikz=cexp(I*(*kz)*(*z)); 
   double jn=gsl_sf_bessel_Jn(*m,rho*(*gamma));
   double complex u=jn*eiph*eikz;
   return(u);
}
//------------------------------------------------------------------------------
// RECTANGULAR WAVE GUIDE - TM MODE
//------------------------------------------------------------------------------
void rwg_tmm(
      double *kx, double *ky, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double gamma2=*kx*(*kx)+*ky*(*ky);
   double k=sqrt(*kx*(*kx)+*ky*(*ky)+*kz*(*kz));
   //
   double ckx=cos(*kx*(*x));
   double skx=sin(*kx*(*x));
   double cky=cos(*ky*(*y));
   double sky=sin(*ky*(*y));
   //
   double complex ETMx=I*(*kz*(*kx)/gamma2)*ckx*sky*cexp(I*(*kz)*(*z));
   double complex ETMy=I*(*kz*(*ky)/gamma2)*skx*cky*cexp(I*(*kz)*(*z));
   //
   double complex ETMm=(ETMx-I*ETMy)/sqrt(2.0);
   double complex ETMz=sin(*kx*(*x))*sin(*ky*(*y))*cexp(I*(*kz)*(*z));
   double complex ETMp=(ETMx+I*ETMy)/sqrt(2.0);
   //
   double complex HTMx=-I*(k*(*ky)/gamma2)*skx*cky*cexp(I*(*kz)*(*z));
   double complex HTMy= I*(k*(*kx)/gamma2)*ckx*sky*cexp(I*(*kz)*(*z));
   //
   double complex HTMm=(HTMx-I*HTMy)/sqrt(2.0);
   double complex HTMz=0.0+I*0.0;
   double complex HTMp=(HTMx+I*HTMy)/sqrt(2.0);
   //
   Em=&ETMm;
   Ez=&ETMz;
   Ep=&ETMp;
   //
   Hm=&HTMm;
   Hz=&HTMz;
   Hp=&HTMp;
}
//------------------------------------------------------------------------------
// RECTANGULAR WAVE GUIDE - TE MODE
//------------------------------------------------------------------------------
void rwg_tem(
      double *kx, double *ky, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double gamma2=*kx*(*kx)+*ky*(*ky);
   double k=sqrt(*kx*(*kx)+*ky*(*ky)+*kz*(*kz));
   //
   double ckx=cos(*kx*(*x));
   double skx=sin(*kx*(*x));
   double cky=cos(*ky*(*y));
   double sky=sin(*ky*(*y));
   //
   double complex ETEx=-I*(k*(*ky)/gamma2)*ckx*sky*cexp(I*(*kz)*(*z));
   double complex ETEy= I*(k*(*kx)/gamma2)*skx*cky*cexp(I*(*kz)*(*z));
   //
   *Em=(ETEx-I*ETEy)/sqrt(2.0);
   *Ez=0.0+I*0.0;
   *Ep=(ETEx+I*ETEy)/sqrt(2.0);
   //
   double complex HTEx=-I*(*kz*(*kx)/gamma2)*skx*cky*cexp(I*(*kz)*(*z));
   double complex HTEy=-I*(*kz*(*ky)/gamma2)*ckx*sky*cexp(I*(*kz)*(*z));
   //
   *Hm=(HTEx-I*HTEy)/sqrt(2.0);
   *Hz=cos(*kx*(*x))*cos(*ky*(*y))*cexp(I*(*kz)*(*z));
   *Hp=(HTEx+I*HTEy)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void wfd_rwg(
                    int *TE,
                    int *nx, int *ny, int *nz,
                    double *kx, double *ky, double *kz,
                    double *x,  double *y,  double *z,
                    double *rx, double *ry, double *rz,
                    double complex *Hm, double complex *Hz, double complex *Hp,
                    double complex *Em, double complex *Ez, double complex *Ep
                    ){ 
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           if(*TE==1){
              rwg_tem(kx,ky,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Hm[i],&Hz[i],&Hp[i],
                    &Em[i],&Ez[i],&Ep[i]);
           }
           if(*TE==0){
              rwg_tmm(kx,ky,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Hm[i],&Hz[i],&Hp[i],
                    &Em[i],&Ez[i],&Ep[i]);
           }
           i++; 
         }    
      }    
   }    
}
//------------------------------------------------------------------------------
// CYLINDRICAL WAVE GUIDE - TM MODE
//------------------------------------------------------------------------------
void cwg_all(
      int *MD,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double ctg=*kz/(*gamma);
   double csc=k/(*gamma);
   int msm=*M-*S;
   int msp=*M+*S;
   //
   *Em= psi_mkr(&msm,S,gamma,kz,x,y,z)*( I*(*S)*ctg)/sqrt(2.0);
   *Ez= psi_mkr(M,   S,gamma,kz,x,y,z);
   *Ep= psi_mkr(&msp,S,gamma,kz,x,y,z)*(-I*(*S)*ctg)/sqrt(2.0);
   // If TM, H <- -H, E <- E
   // If TE, E <-  H, H <- E
   *Hm=-(*MD)*psi_mkr(&msm,S,gamma,kz,x,y,z)*(*S*csc)/sqrt(2.0);
   *Hz=0;
   *Hp=-(*MD)*psi_mkr(&msp,S,gamma,kz,x,y,z)*(*S*csc)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void wfd_cwg(
                    int *TE, int *m, int *s,
                    int *nx, int *ny, int *nz,
                    double *gamma, double *kz,
                    double *x,  double *y,  double *z,
                    double *rx, double *ry, double *rz,
                    double complex *Hm, double complex *Hz, double complex *Hp,
                    double complex *Em, double complex *Ez, double complex *Ep
                    ){ 
   int MD;
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           if(*TE==1){
              // MODO TE
              MD=1;
              cwg_all(&MD,m,s,
                    gamma,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Em[i],&Ez[i],&Ep[i],
                    &Hm[i],&Hz[i],&Hp[i]);
           }
           if(*TE==0){
              // MODO TM
              MD=-1;
              cwg_all(&MD,m,s,
                    gamma,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Hm[i],&Hz[i],&Hp[i],
                    &Em[i],&Ez[i],&Ep[i]);
           }
           i++; 
         }    
      }    
   }    
}
//------------------------------------------------------------------------------
// BESSEL BEAMS Z ORIENTED
//------------------------------------------------------------------------------
void def_bbz(
      int *MD,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double cth=*kz/k;
   double sth=(*gamma)/k;
   int msm=*M-*S;
   int msp=*M+*S;
   // TM MODE
   *Em=(*MD)*psi_mkr(&msm,S,gamma,kz,x,y,z)*( I*(*S)*cth*sth)/sqrt(2.0);
   *Ez=(*MD)*psi_mkr(M,   S,gamma,kz,x,y,z)*sth*sth;
   *Ep=(*MD)*psi_mkr(&msp,S,gamma,kz,x,y,z)*(-I*(*S)*cth*sth)/sqrt(2.0);
   //
   *Hm=psi_mkr(&msm,S,gamma,kz,x,y,z)*(*S*sth)/sqrt(2.0);
   *Hz=0;
   *Hp=psi_mkr(&msp,S,gamma,kz,x,y,z)*(*S*sth)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void wfd_bbz(
                    int *TE, int *m, int *s,
                    int *nx, int *ny, int *nz,
                    double *gamma, double *kz,
                    double *x,  double *y,  double *z,
                    double *rx, double *ry, double *rz,
                    double complex *Hm, double complex *Hz, double complex *Hp,
                    double complex *Em, double complex *Ez, double complex *Ep
                    ){ 
//------------------------------------------------------------------------------
   int MD;
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           if(*TE==1){
              MD=1;
              def_bbz(&MD,m,s,
                    gamma,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Hm[i],&Hz[i],&Hp[i],
                    &Em[i],&Ez[i],&Ep[i]);
           }
           if(*TE==0){
              MD=-1;
              def_bbz(&MD,m,s,
                    gamma,kz,&x[ix],
                    &y[iy],&z[iz],&Em[i],&Ez[i],
                    &Ep[i],&Hm[i],&Hz[i],&Hp[i]);
           }
           i++; 
         }    
      }    
   }    
}

//------------------------------------------------------------------------------
// BESSEL BEAMS P ORIENTED
//------------------------------------------------------------------------------
void def_bbp(
      int *P,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double cth=*kz/k;
   double sth=(*gamma)/k;

   int pp1=1+*P;
   int pm1=1-*P;
   int mm1=*M-1;
   int mspm1=*M-1+*S*(*P-1);
   int mspm0=*M-1+*S*(*P);
   int mspp1=*M-1+*S*(*P+1);

   *Em=pp1*psi_mkr(&mm1,S,gamma,kz,x,y,z)/2
      -*P*(sth*sth)*psi_mkr(&mspm1,S,gamma,kz,x,y,z)/2;
   *Ez= -I*(*P)*(*S)*cth*sth*psi_mkr(&mspm0,S,gamma,kz,x,y,z)/sqrt(2.0);
   *Ep=pm1*psi_mkr(&mm1,S,gamma,kz,x,y,z)/2
      +*P*(sth*sth)*psi_mkr(&mspp1,S,gamma,kz,x,y,z)/2;
   *Hm=-I*(*P)*cth*psi_mkr(&mm1,S,gamma,kz,x,y,z)*pp1/2.0;
   *Hz=-(*S)*sth*psi_mkr(&mspm0,S,gamma,kz,x,y,z)/sqrt(2.0);
   *Hp=-I*(*P)*cth*psi_mkr(&mm1,S,gamma,kz,x,y,z)*pm1/2.0;
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void wfd_bbp(
                    int *p, int *m, int *s, 
                    int *nx, int *ny, int *nz,
                    double *gamma, double *kz,
                    double *x,  double *y,  double *z,
                    double *rx, double *ry, double *rz,
                    double complex *Hm, double complex *Hz, double complex *Hp,
                    double complex *Em, double complex *Ez, double complex *Ep
                    ){ 
//------------------------------------------------------------------------------
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           def_bbp(p,m,s,
                 gamma,kz,
                 &x[ix],&y[iy],&z[iz],
                 &Hm[i],&Hz[i],&Hp[i],
                 &Em[i],&Ez[i],&Ep[i]);
           i++; 
         }    
      }    
   }    
//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
// VECTOR SPHERICAL WAVE FUNCTIONS CALCULATIONS
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// CONSTANTS FOR Qlm (NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS)
//------------------------------------------------------------------------------
double alfaQ(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt(((2*l-1)*(2*l+1))/((l-m)*(l+m)));
}
//------------------------------------------------------------------------------
double betaQ(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt((2*l+1)/(2*l-3))*sqrt(((l+m-1)*(l-m-1))/((l-m)*(l+m)));
}
//------------------------------------------------------------------------------
double gammaQ(int lo){
   double l=lo;
   return sqrt((2*l+1)/(2*l));
}
//------------------------------------------------------------------------------
double deltaQ(int lo){
   double l=lo;
   return sqrt(2*l+1);
}
//------------------------------------------------------------------------------
// POSITIONER IN THE (l,m) LIST
//------------------------------------------------------------------------------
int jlm(int l, int m){
   if(abs(m)>l){
      return 0;
   }else{
      return l*(l+1)+m;
   }
}
//------------------------------------------------------------------------------
// POSITIONER OF ith,jth, kth element in the list  l from 0 to ni*nj*nk
//------------------------------------------------------------------------------
/*In general, the offset formula for an array with dimensions
* [d1 d2 d3 ... dn] using any subscripts (s1 s2 s3 ... sn) is
* (sn-1)(dn-1)(dn-2)...(d1)+(sn-1-1)(dn-2)...(d1)+...+(s2-1)(d1)+s1
* Because of this scheme, you can index an array using any number of subscripts.
*
* i<-iz+nz*(iy-1)+nz*ny*(ix-1) : comecando de 1
* i<-iz+nz*iy+nz*ny*ix         : comecando de 0
*/
int lijk(int ix, int ni, int iy,int ny, int iz, int nz){
   int l=iz+nz*iy+nz*ny*ix;
   return(l);

}
//------------------------------------------------------------------------------
// CONSTANTS FOR \vec{X}_{lm}
//------------------------------------------------------------------------------
double cp(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt(l*(l+1)-m*(m+1));
}
//------------------------------------------------------------------------------
double cm(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt(l*(l+1)-m*(m-1));
}
//------------------------------------------------------------------------------
// POSITION DEPENDENT CALCULATIONS - MULTIPLE LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void def_vsw(double *k,double *x, double *y, double *z,
                    int *lmax, 
                    double complex *GTE, double complex *GTM,
                    double complex *Em, double complex *Ez, double complex *Ep,
                    double complex *Hm, double complex *Hz, double complex *Hp
                    ){
//------------------------------------------------------------------------------
   int l;
//   int LMAX=*lmax*(*lmax+2);
   int LMAXE=(*lmax+1)*(*lmax+3);
   double cph; 
   double sph;
   double rho=sqrt(*x*(*x)+*y*(*y));
   double r=sqrt(rho*rho+*z*(*z));
   double sth=rho/r;
   double cth=*z/r;
//------------------------------------------------------------------------------
   if((*x==0)&&(*y==0)){
      cph=1;
      sph=0;
   }else{
      cph=*x/rho;
      sph=*y/rho;
   }
//------------------------------------------------------------------------------
   // Spherical Bessel Funtions
   double JLM[*lmax+2];
   gsl_sf_bessel_jl_steed_array(*lmax+1,*k*r,JLM);
//------------------------------------------------------------------------------
   // Qlm - First 4 terms
   double Qlm[LMAXE];
   Qlm[jlm(0, 0)]=1/sqrt(4*M_PI);
   Qlm[jlm(1, 1)]=-gammaQ(1)*sth*Qlm[jlm(0,0)]; // Q11
   Qlm[jlm(1, 0)]=sqrt(3.0)*cth*Qlm[jlm(0,0)];  // Q10
   Qlm[jlm(1,-1)]=-Qlm[jlm(1,1)];               // Q11*(-1)
//------------------------------------------------------------------------------
   // Complex Exponencial for m=-1,0,1
   double complex Eim[2*(*lmax)+3];
   Eim[*lmax-1]=(cph-I*sph);
   Eim[*lmax  ]=1+I*0;
   Eim[*lmax+1]=(cph+I*sph);
//------------------------------------------------------------------------------
   // Ylm - First 4 terms
   double complex Ylm[LMAXE];
   Ylm[jlm(0, 0)]=Qlm[jlm(0, 0)];
   Ylm[jlm(1,-1)]=Qlm[jlm(1,-1)]*Eim[*lmax-1];
   Ylm[jlm(1, 0)]=Qlm[jlm(1, 0)];
   Ylm[jlm(1, 1)]=Qlm[jlm(1, 1)]*Eim[*lmax+1];
//------------------------------------------------------------------------------
   // VECTOR SPHERICAL HARMONICS
//------------------------------------------------------------------------------
   // r
   double complex rm=sth*(cph-I*sph)/sqrt(2);
   double complex rz=cth;
   double complex rp=sth*(cph+I*sph)/sqrt(2);
   // X
   double complex XM;
   double complex XZ;
   double complex XP;
   // Y
   double complex YM;
   double complex YZ;
   double complex YP;
   // V
   double complex VM;
   double complex VZ;
   double complex VP;
//------------------------------------------------------------------------------
   // HANSEN MULTIPOLES
//------------------------------------------------------------------------------
   // M
   double complex MM;
   double complex MZ;
   double complex MP;
   // N
   double complex NM;
   double complex NZ;
   double complex NP;
//------------------------------------------------------------------------------
   // OTHERS
//------------------------------------------------------------------------------
   double kl;
//------------------------------------------------------------------------------
   // MAIN LOOP
   for(l=1;l<=(*lmax);l++){
//------------------------------------------------------------------------------
      //Qlm extremos positivos um passo a frente
      Qlm[jlm(l+1, l+1)]=-gammaQ(l+1)*sth*Qlm[jlm(l,l)];
      Qlm[jlm(l+1, l  )]= deltaQ(l+1)*cth*Qlm[jlm(l,l)];
      //Qlm extremos negativos um passo a frente
      Qlm[jlm(l+1,-l-1)]=pow(-1,l+1)*Qlm[jlm(l+1, l+1)];
      Qlm[jlm(l+1,-l  )]=pow(-1,l  )*Qlm[jlm(l+1, l  )];
      // Exponenciais um passo a frente
      Eim[*lmax+l+1]=Eim[*lmax+l]*(cph+I*sph);
      Eim[*lmax-l-1]=Eim[*lmax-l]*(cph-I*sph);
      // Harmonicos esfericos extremos um passo a frente
      Ylm[jlm(l+1, l+1)]=Qlm[jlm(l+1, l+1)]*Eim[*lmax+l+1];
      Ylm[jlm(l+1, l  )]=Qlm[jlm(l+1, l  )]*Eim[*lmax+l  ];
      Ylm[jlm(l+1,-l-1)]=Qlm[jlm(l+1,-l-1)]*Eim[*lmax-l-1];
      Ylm[jlm(l+1,-l  )]=Qlm[jlm(l+1,-l  )]*Eim[*lmax-l  ];
      // others
      kl=1/(sqrt(l*(l+1)));
//------------------------------------------------------------------------------
      for(int m=l; m>=(-l); m--){
         // CALCULATIONS OF SSH
         if(m>=0){
            Qlm[jlm(l+1, m)]=alfaQ(l+1,m)*cth*Qlm[jlm(l,m)]
               -betaQ(l+1,m)*Qlm[jlm(l-1,m)];
            Qlm[jlm(l+1,-m)]=pow(-1,m)*Qlm[jlm(l+1, m)];
            Ylm[jlm(l+1, m)]=Qlm[jlm(l+1, m)]*Eim[*lmax+m];
            Ylm[jlm(l+1,-m)]=Qlm[jlm(l+1,-m)]*Eim[*lmax-m];
         }
         // CALCULATIONS OF VSH
         // X
         XM=kl*cm(l,m)*Ylm[jlm(l,m-1)]/sqrt(2);
         XZ=kl*m*Ylm[jlm(l,m  )];
         XP=kl*cp(l,m)*Ylm[jlm(l,m+1)]/sqrt(2);
         // Y
         YM=rm*Ylm[jlm(l,m)];
         YZ=rz*Ylm[jlm(l,m)];
         YP=rp*Ylm[jlm(l,m)];
         // V
         VM=rm*XZ-rz*XM;
         VZ=rp*XM-rm*XP;
         VP=rz*XP-rp*XZ;
         // CALCULATION OF HANSEM MULTIPOLES
         // M
         MM=JLM[l]*XM;
         MZ=JLM[l]*XZ;
         MP=JLM[l]*XP;
         // N
         NM=((1.*l+1.)*JLM[l-1]-l*JLM[l+1])*VM/(2.*l+1.)
            +sqrt(l*(l+1.))*(JLM[l-1]+JLM[l+1])*YM/(2.*l+1.);
         NZ=((1.*l+1.)*JLM[l-1]-l*JLM[l+1])*VZ/(2.*l+1.)
            +sqrt(l*(l+1.))*(JLM[l-1]+JLM[l+1])*YZ/(2.*l+1.);
         NP=((1.*l+1.)*JLM[l-1]-l*JLM[l+1])*VP/(2.*l+1.)
            +sqrt(l*(l+1.))*(JLM[l-1]+JLM[l+1])*YP/(2.*l+1.);
         // CALCULATION OF THE ELECTROMAGNETIC FIELDS
         *Em=*Em+MM*GTE[jlm(l,m)]-NM*GTM[jlm(l,m)]; 
         *Ez=*Ez+MZ*GTE[jlm(l,m)]-NZ*GTM[jlm(l,m)];
         *Ep=*Ep+MP*GTE[jlm(l,m)]-NP*GTM[jlm(l,m)];
         *Hm=*Hm+MM*GTM[jlm(l,m)]+NM*GTE[jlm(l,m)]; 
         *Hz=*Hz+MZ*GTM[jlm(l,m)]+NZ*GTE[jlm(l,m)];
         *Hp=*Hp+MP*GTM[jlm(l,m)]+NP*GTE[jlm(l,m)];
      }
   }
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void pwe_vsw(double *k,double *x, double *y, double *z,
                    int *lmax, int *nx, int *ny, int *nz, 
                    double complex *GTE, double complex *GTM,
                    double *rx, double *ry, double *rz,
                    double complex *Em, double complex *Ez, double complex *Ep,
                    double complex *Hm, double complex *Hz, double complex *Hp
                    ){ 
//------------------------------------------------------------------------------
   int i=0;
   #pragma omp parallel for
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           /*In general, the offset formula for an array with dimensions 
            * [d1 d2 d3 ... dn] using any subscripts (s1 s2 s3 ... sn) is
            * (sn-1)(dn-1)(dn-2)...(d1)+(sn-1-1)(dn-2)...(d1)+...+(s2-1)(d1)+s1
            * Because of this scheme, you can index an array using any number of subscripts. 
            * 
            * i<-iz+nz*(iy-1)+nz*ny*(ix-1) : comecando de 1
            * i<-iz+nz*iy+nz*ny*ix         : comecando de 0
            */
           i=iz+(*nz)*iy+(*nz)*(*ny)*ix; 
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           def_vsw(k,
                 &x[ix],&y[iy],&z[iz],
                 lmax,GTE,GTM,
                 &Em[i],&Ez[i],&Ep[i],
                 &Hm[i],&Hz[i],&Hp[i]);
           //i++; 
         } /* for iz */    
      } /* for iy */    
   } /* for ix */    
} /*void pwe_vsw */
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
