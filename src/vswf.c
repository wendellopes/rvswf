/*------------------------------------------------------------------------------
 * VSWF: Vector Spherical Wave Functions dynamic libraries written in C.
 * author: Wendel Lopes Moreira <wendellopes@gmail.com>
 * date: 2014-06-21
 * version 1.2
 * depends: Gnu Scientific Library <http://www.gnu.org/software/gsl/>
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
// #include <float.h> Was not working, changed 1e-37 by 1e-37
#include <gsl/gsl_sf_bessel.h>
#ifdef _OPENMP
   #include <omp.h>
#endif
/*----------------------------------------------------------------------------*
 *                             BESSEL                                         *
 *----------------------------------------------------------------------------*/
//------------------------------------------------------------------------------
//# Auxliary 
//------------------------------------------------------------------------------
double lcfe_afs(/* FUNCTION */
      int n,
      double x
   ){
   return(n/x);
}
//------------------------------------------------------------------------------
double complex lcfc_afs(/* FUNCTION */
      int n,
      double complex x
   ){
   return(n/x);
}
//------------------------------------------------------------------------------
//# Logarithmic Derivative of Cylindrical Bessel [OK] == REAL
//------------------------------------------------------------------------------
void lcfe_cbl(/* FUNCTION */
      int *n,
      double *x,
      int *NMAX,
      double *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfe_afs(*n,*x);
   if(fabs(*fn)<eo){*fn=eo;}
   double Pn=*fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;
   while(fabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      int u=2*(*n+j);
      bn=lcfe_afs(u,*x);
      Pn=bn+an/Pn;
      if(fabs(Pn)<eo){Pn=eo;}// # migth be zero
      Qn=bn+an*Qn;
      if(fabs(Qn)<eo){Qn=eo;}// # migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# Logarithmic Derivative of Cylindrical Bessel [OK] == COMPLEX
//------------------------------------------------------------------------------
void lcfc_cbl(/* FUNCTION */
      int *n,
      double complex *x,
      int *NMAX,
      double complex *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfc_afs(*n,*x);
   if(cabs(*fn)<eo){*fn=eo;}
   double complex Pn=*fn;
   double complex Qn=0.0;
   // Loop Parameters
   int j=0;
   double complex Dn=10.0;
   double complex an;
   double complex bn;
   while(cabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1+0*I;
      int u=2*(*n+j);
      bn=lcfc_afs(u,*x);
      Pn=bn+an/Pn;
      if(cabs(Pn)<eo){Pn=eo;}// # migth be zero
      Qn=bn+an*Qn;
      if(cabs(Qn)<eo){Qn=eo;}// # migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# J_{n}/J_{n+1} [OK] DIRECT == REAL
//------------------------------------------------------------------------------
void lcfe_cbd(/* FUNCTION */
      int *n,
      double *x,
      int *NMAX,
      double *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfe_afs(2*(*n+1),*x);
   if(fabs(*fn)<eo){*fn=eo;} // migth be zero
   double Pn=*fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;
   while(fabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfe_afs(2*(*n+j+1),*x);
      Pn=bn+an/Pn;
      if(fabs(Pn)<eo){Pn=eo;} // migth be zero
      Qn=bn+an*Qn;
      if(fabs(Qn)<eo){Qn=eo;} // migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# J_{n}/J_{n+1} [OK] DIRECT == COMPLEX
//------------------------------------------------------------------------------
void lcfc_cbd(/* FUNCTION */
      int *n,
      double complex *x,
      int *NMAX,
      double complex *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfc_afs(2*(*n+1),*x);
   if(cabs(*fn)<eo){*fn=eo*(1+I);} // migth be zero
   double complex Pn=*fn;
   double complex Qn=0.0;
   // Loop Parameters
   int j=0;
   double complex Dn=10.0;
   double complex an;
   double complex bn;
   while(cabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfc_afs(2*(*n+j+1),*x);
      Pn=bn+an/Pn;
      if(cabs(Pn)<eo){Pn=eo;} // migth be zero
      Qn=bn+an*Qn;
      if(cabs(Qn)<eo){Qn=eo;} // migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# J_{n+1}/J_{n} [OK] BARNETT == REAL
//------------------------------------------------------------------------------
void lcfe_cbi(/* FUNCTION */
      int *n,
      double *x,
      int *NMAX,
      double *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=0.0;
   if(fabs(*fn)<eo){*fn=eo;} // migth be zero
   double Pn=*fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;
   while(fabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      if(j==0){
         an=1;
      }else{
         an=-1;
      }
      j=j+1;
      bn=lcfe_afs(2*(*n+j),*x);
      Pn=bn+an/Pn;
      if(fabs(Pn)<eo){Pn=eo;} //# migth be zero
      Qn=bn+an*Qn;
      if(fabs(Qn)<eo){Qn=eo;} //# migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# J_{n+1}/J_{n} [OK] BARNETT == COMPLEX
//------------------------------------------------------------------------------
void lcfc_cbi(/* FUNCTION */
      int *n,
      double complex *x,
      int *NMAX,
      double complex *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   if(cabs(*fn)<eo){*fn=eo;} // migth be zero
   double complex Pn=*fn;
   double complex Qn=0.0;
   // Loop Parameters
   int j=0;
   double complex Dn=10.0;
   double complex an;
   double complex bn;
   while(cabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      if(j==0){
         an=1;
      }else{
         an=-1;
      }
      j=j+1;
      bn=lcfc_afs(2*(*n+j),*x);
      Pn=bn+an/Pn;
      if(cabs(Pn)<eo){Pn=eo;} //# migth be zero
      Qn=bn+an*Qn;
      if(cabs(Qn)<eo){Qn=eo;} //# migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# Logarithmic Derivative of Riccati-Bessel [OK] == REAL
//------------------------------------------------------------------------------
void lcfe_rbl(/* FUNCTION */
      int *n,
      double *x,
      int *NMAX,
      double *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;;
   double ACC=1e-50;
   *fn=lcfe_afs(*n+1,*x);
   if(fabs(*fn)<eo){*fn=eo;} // migth be zero;
   double Pn=*fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;
   while(fabs(Dn-1.)>ACC){;
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfe_afs(2*(*n+j)+1,*x);
      Pn=bn+an/Pn;
      if(fabs(Pn)<eo){Pn=eo;} // migth be zero;
      Qn=bn+an*Qn;
      if(fabs(Qn)<eo){Qn=eo;} // migth be zero;
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# Logarithmic Derivative of Riccati-Bessel [OK] == COMPLEX
//------------------------------------------------------------------------------
void lcfc_rbl(/* FUNCTION */
      int *n,
      double complex *x,
      int *NMAX,
      double complex *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;;
   double ACC=1e-50;
   *fn=lcfc_afs(*n+1,*x);
   if(cabs(*fn)<eo){*fn=eo;} // migth be zero;
   double complex Pn=*fn;
   double complex Qn=0.0;
   // Loop Parameters
   int j=0;
   double complex Dn=10.0;
   double complex an;
   double complex bn;
   while(cabs(Dn-1.)>ACC){;
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfc_afs(2*(*n+j)+1,*x);
      Pn=bn+an/Pn;
      if(cabs(Pn)<eo){Pn=eo;} // migth be zero;
      Qn=bn+an*Qn;
      if(cabs(Qn)<eo){Qn=eo;} // migth be zero;
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# Logarithmic Derivative of Spherical Bessel [OK] == REAL
//------------------------------------------------------------------------------
void lcfe_sbl(/* FUNCTION */
      int *n,
      double *x,
      int *NMAX,
      double *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfe_afs(*n,*x);
   if(fabs(*fn)<eo){*fn=eo;} // migth be zero;
   double Pn=*fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;
   while(fabs(Dn-1.)>ACC){;
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfe_afs(2*(*n+j)+1,*x);;
      Pn=bn+an/Pn;
      if(fabs(Pn)<eo){Pn=eo;} // migth be zero;
      Qn=bn+an*Qn;
      if(fabs(Qn)<eo){Qn=eo;} // migth be zero;
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# Logarithmic Derivative of Spherical Bessel [OK] == COMPLEX
//------------------------------------------------------------------------------
void lcfc_sbl(/* FUNCTION */
      int *n,
      double complex *x,
      int *NMAX,
      double complex *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfc_afs(*n,*x);
   if(cabs(*fn)<eo){*fn=eo;} // migth be zero;
   double complex Pn=*fn;
   double complex Qn=0.0;
   // Loop Parameters
   int j=0;
   double complex Dn=10.0;
   double complex an;
   double complex bn;
   while(cabs(Dn-1.)>ACC){;
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfc_afs(2*(*n+j)+1,*x);;
      Pn=bn+an/Pn;
      if(cabs(Pn)<eo){Pn=eo;} // migth be zero;
      Qn=bn+an*Qn;
      if(cabs(Qn)<eo){Qn=eo;} // migth be zero;
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# j_{n}/j_{n+1} [OK] DIRECT == REAL
//------------------------------------------------------------------------------
void lcfe_sbd(/* FUNCTION */
      int *n,
      double *x,
      int *NMAX,
      double *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfe_afs(2*(*n+1)+1,*x);
   if(fabs(*fn)<eo){*fn=eo;} // migth be zero
   double Pn=*fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;   
   while(fabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfe_afs(2*(*n+j+1)+1,*x);
      Pn=bn+an/Pn;
      if(fabs(Pn)<eo){Pn=eo;} // migth be zero
      Qn=bn+an*Qn;
      if(fabs(Qn)<eo){Qn=eo;} // migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# j_{n}/j_{n+1} [OK] DIRECT == COMPLEX
//------------------------------------------------------------------------------
void lcfc_sbd(/* FUNCTION */
      int *n,
      double complex *x,
      int *NMAX,
      double complex *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=lcfc_afs(2*(*n+1)+1,*x);
   if(cabs(*fn)<eo){*fn=eo;} // migth be zero
   double complex Pn=*fn;
   double complex Qn=0.0;
   // Loop Parameters
   int j=0;
   double complex Dn=10.0;
   double complex an;
   double complex bn;   
   while(cabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      j=j+1;
      an=-1;
      bn=lcfc_afs(2*(*n+j+1)+1,*x);
      Pn=bn+an/Pn;
      if(cabs(Pn)<eo){Pn=eo;} // migth be zero
      Qn=bn+an*Qn;
      if(cabs(Qn)<eo){Qn=eo;} // migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# j_{n+1}/j_{n} [OK] BARNETT == REAL
//------------------------------------------------------------------------------
void lcfe_sbi(/* FUNCTION */
      int *n,
      double *x,
      int *NMAX,
      double *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=0.0;
   if(fabs(*fn)<eo){*fn=eo;} // migth be zero
   double Pn=*fn;
   double Qn=0.0;
   // Loop Parameters
   int j=0;
   double Dn=10.0;
   double an;
   double bn;
   while(fabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      if(j==0){
         an=1;
      }else{
         an=-1;
      }
      j=j+1;
      bn=lcfe_afs(2*(*n+j)+1,*x);
      Pn=bn+an/Pn;
      if(fabs(Pn)<eo){Pn=eo;} //# migth be zero
      Qn=bn+an*Qn;
      if(fabs(Qn)<eo){Qn=eo;} //# migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
//# j_{n+1}/j_{n} [OK] BARNETT == COMPLEX
//------------------------------------------------------------------------------
void lcfc_sbi(/* FUNCTION */
      int *n,
      double complex *x,
      int *NMAX,
      double complex *fn
   ){
   //-----------------------------------
   const double eo = 1e-37;
   double ACC=1e-50;
   *fn=0.0;
   if(cabs(*fn)<eo){*fn=eo;} // migth be zero
   double complex Pn=*fn;
   double complex Qn=0.0;
   // Loop Parameters
   int j=0;
   double complex Dn=10.0;
   double complex an;
   double complex bn;
   while(cabs(Dn-1.)>ACC){
      if(j>*NMAX){
         break;
      }
      if(j==0){
         an=1;
      }else{
         an=-1;
      }
      j=j+1;
      bn=lcfc_afs(2*(*n+j)+1,*x);
      Pn=bn+an/Pn;
      if(cabs(Pn)<eo){Pn=eo;} //# migth be zero
      Qn=bn+an*Qn;
      if(cabs(Qn)<eo){Qn=eo;} //# migth be zero
      Qn=1/Qn;
      Dn=Pn*Qn;
      *fn=*fn*Dn;
   }
   *NMAX=j;
}
//------------------------------------------------------------------------------
// Auxiliary func for calculation of Bessel functions by expansion
double KJ(/* FUNCTION */
      double x,
      int j){
   return(x/(1.*j));
}
//------------------------------------------------------------------------------
// Cylindrical Bessel func ratio and logarithmic derivative 
void lcfa_cyl(/* FUNCTION */
      int *nmax,
      double *x,
      double *gn,
      double *Dn,
      int *NMAX
   ){
//--------------------------------------
   lcfe_cbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfe_cbd(nmax,x,NMAX,&gn[*nmax]);
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=lcfe_afs(n  ,*x)+Dn[n];
      Dn[n-1]=lcfe_afs(n-1,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Cylindrical Bessel func ratio and logarithmic derivative -- COMPLEX
void lcfc_cyl(/* FUNCTION */
      int *nmax,
      double complex *x,
      double complex *gn,
      double complex *Dn,
      int *NMAX
   ){
//--------------------------------------
   lcfc_cbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfc_cbd(nmax,x,NMAX,&gn[*nmax]);
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=lcfc_afs(n  ,*x)+Dn[n];
      Dn[n-1]=lcfc_afs(n-1,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Spherical Bessel func ratio and logarithmic derivative 
void lcfa_sph(/* FUNCTION */
      int *nmax,
      double *x,
      double *gn,
      double *Dn,
      int *NMAX
   ){
//--------------------------------------
   lcfe_sbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfe_sbd(nmax,x,NMAX,&gn[*nmax]);
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=lcfe_afs(n+1,*x)+Dn[n];
      Dn[n-1]=lcfe_afs(n-1,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Spherical Bessel func ratio and logarithmic derivative -- COMPLEX
void lcfc_sph(/* FUNCTION */
      int *nmax,
      double complex *x,
      double complex *gn,
      double complex *Dn,
      int *NMAX
   ){
//--------------------------------------
   lcfc_sbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfc_sbd(nmax,x,NMAX,&gn[*nmax]);
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=lcfc_afs(n+1,*x)+Dn[n];
      Dn[n-1]=lcfc_afs(n-1,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Riccati-Bessel func ratio and logarithmic derivative 
void lcfa_ric(/* FUNCTION */
      int *nmax,
      double *x,
      double *gn,
      double *Dn,
      int *NMAX
   ){
//--------------------------------------
   lcfe_rbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfe_sbd(nmax,x,NMAX,&gn[*nmax]);
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=lcfe_afs(n,*x)+Dn[n];
      Dn[n-1]=lcfe_afs(n,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Riccati-Bessel func ratio and logarithmic derivative -- COMPLEX
void lcfc_ric(/* FUNCTION */
      int *nmax,
      double complex *x,
      double complex *gn,
      double complex *Dn,
      int *NMAX
   ){
//--------------------------------------
   lcfc_rbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfc_sbd(nmax,x,NMAX,&gn[*nmax]);
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=lcfc_afs(n,*x)+Dn[n];
      Dn[n-1]=lcfc_afs(n,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Spherical Bessel func j_0(x)
void bess_szr(/* FUNCTION */
      double *x, 
      double *j0
   ){
//--------------------------------------
   if(*x>1.){
      *j0=sin(*x)/(*x);
   }else{
      int j;
      *j0=1.0;
      for(j=10;j>0;j--){
         *j0=1.-*j0*KJ(*x,2*j)*KJ(*x,2*j+1);
      }
   }
}
//------------------------------------------------------------------------------
// Spherical Bessel func j_1(x)
void bess_sun(/* FUNCTION */
      double *x, 
      double *j1
   ){
//--------------------------------------
   if(*x>1.){
      *j1=(sin(*x)/(*x)-cos(*x))/(*x);
   }else{
      int j;
      *j1=0.0;
      for(j=10;j>0;j--){
         *j1=(2.*j)/(2*j+1.)-*j1*KJ(*x,2*j+1)*KJ(*x,2*j+2);
      }
      *j1=(*x)*(*j1)/2.0;
   }
}
//------------------------------------------------------------------------------
//Starting values for downward recurrence, Cylindrical Bessel func
//For calculation of J_n(x), one can calculate by downward recurrence from
//   J_{n-dig}(x), where dig is the number of precision in J_n(x).
void bess_csv(/* FUNCTION */
      int *nmax,
      int *dig,
      double *x,
      double *JN,
      double *DN
   ){
//--------------------------------------
   double Jn[*dig], Dn[*dig];
   Jn[*dig]=0.0;  // n-th Bessel func
   Dn[*dig]=1.0;  // Derivative of the n-th Bessel func
   int n,j;
   for(n=*dig;n>0;n--){
      // Downward Recursion
      Jn[n-1]=lcfe_afs(*nmax+n  ,*x)*Jn[n  ]+Dn[n];
      Dn[n-1]=lcfe_afs(*nmax+n-1,*x)*Jn[n-1]-Jn[n];
      //Renormalization may be required
      if(fabs(Jn[n-1])>1e+2){
         for(j=*dig;j>=n-1;j--){
            Dn[j]=Dn[j]/Jn[n-1];
            Jn[j]=Jn[j]/Jn[n-1];
         }
      }
   }
   *JN=Jn[0];
   *DN=Dn[0];
}
//------------------------------------------------------------------------------
// Array of cylindrical Bessel func
void bess_cyl(/* FUNCTION */
      int *nmax, 
      double *x, 
      double *Jn, 
      double *Dn, 
      int *NMAX
   ){
//--------------------------------------
   Jn[*nmax]=0.0;  // n-th Bessel func
   Dn[*nmax]=1.0;  // Derivative of the n-th Bessel func
   int dig=15;     // Number of digits of precision
   int n,j;
   bess_csv(nmax,&dig,x,&Jn[*nmax],&Dn[*nmax]);
   //lcfe_cbl(nmax,x,NMAX,&Dn[*nmax]);
   for(n=*nmax;n>0;n--){
      // Downward Recursion
      Jn[n-1]=lcfe_afs(n  ,*x)*Jn[n  ]+Dn[n];
      Dn[n-1]=lcfe_afs(n-1,*x)*Jn[n-1]-Jn[n];
      //Renormalization may be required
      if(fabs(Jn[n-1])>1e2){
         for(j=*nmax;j>=n-1;j--){
            Dn[j]=Dn[j]/Jn[n-1];
            Jn[j]=Jn[j]/Jn[n-1];
         }
      }
   }
   // Normalization factor for the func
   double pf=1.0;
   // Using j0(x) and j1(x) from math.h,
   // otherwise it should be entered from R
   // or calculated internally.
   double J0=j0(*x);
   double J1=j1(*x);
   if(fabs(J0)>1e-10){
      //pf=*J0/Jn[0];
      pf=j0(*x)/Jn[0];
   }else{
      //pf=*J1/Jn[1];
      pf=j1(*x)/Jn[1];
   }
   // Normalization factor for the derivative
   double pd=1.0;
   if(fabs(J1)>1e-10){
      pd=-(J1)/Dn[0];
   }else{
      pd=0.5*(Jn[0]-Jn[2])/Dn[1];
   }
   // Normalizations
   for(n=0;n<=*nmax;n++){
      Jn[n]=Jn[n]*pf;
      Dn[n]=Dn[n]*pd;
   }
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Starting values for downward recurrence, Spherical Bessel func
void bess_ssv(/* FUNCTION */
      int *nmax,
      int *dig,
      double *x,
      double *JN,
      double *CN
   ){
//--------------------------------------
   double jn[*dig], cn[*dig];
   jn[*dig]=1.0;  // n-th Bessel func
   cn[*dig]=0.0;  // Derivative of the n-th Bessel func
   int n,j;
   for(n=*dig;n>0;n--){
      // Downward Recursion
      jn[n-1]=lcfe_afs(*nmax+n+1,*x)*jn[n  ]+cn[n];
      cn[n-1]=lcfe_afs(*nmax+n-1,*x)*jn[n-1]-jn[n];
      //Renormalization may be required
      if(fabs(jn[n-1])>1e+2){
         for(j=n-1;j<=*dig;j--){
            cn[j]=cn[j]/jn[n-1];
            jn[j]=jn[j]/jn[n-1];
         }
      }
   }
   *JN=jn[0];
   *CN=cn[0];
}
//------------------------------------------------------------------------------
// Array of Spherical Bessel func
// NEED COMPLEX VERSION?
void bess_sph(/* FUNCTION */
      int *nmax, 
      double *x, 
      double *jn, 
      double *dn,
      int *NMAX
   ){
//--------------------------------------
   double j0=0.0;  // Normalization values
   double j1=0.0;  // Normalization values
   bess_szr(x,&j0);
   bess_sun(x,&j1);
//--------------------------------------
   int NMAY=*NMAX;
   lcfe_sbl(nmax,x,NMAX,&dn[*nmax]);
   while(*NMAX==NMAY){
      *NMAX=2*NMAY;
      NMAY=*NMAX;
      lcfe_sbl(nmax,x,NMAX,&dn[*nmax]);
   }
//--------------------------------------
   int j,n;
   jn[*nmax] = 1.0; //rn[*nmax];
   //dn[*nmax] = 0.0; //rn[*nmax];
   for(n=*nmax;n>0;n--){
      jn[n-1]=lcfe_afs(n+1,*x)*jn[n  ]+dn[n  ];
      dn[n-1]=lcfe_afs(n-1,*x)*jn[n-1]-jn[n  ];
      // Renormalization
      if(fabs(jn[n-1])>1e2){
         for(j=*nmax;j<=n-1;j--){
            dn[j]=dn[j]/jn[n-1];
            jn[j]=jn[j]/jn[n-1];
         }
      }
   }
//--------------------------------------
   // Normalization factor for the func
   double pj=1.0;
   if(fabs(j0)>1e-10){
      pj=j0/jn[0];
   }else{
      pj=j1/jn[1];
   }
//--------------------------------------
   // Normalization factor for the derivative
   double pd=1.0;
   if(fabs(j1)>1e-10){
      pd=-(j1)/dn[0];
   }else{
      pd=(1./3.)*(jn[0]-2*jn[2])/dn[1];
   }
//--------------------------------------
   // Normalizations
   for(n=0;n<=*nmax;n++){
      jn[n]=jn[n]*pj;
      dn[n]=dn[n]*pd;
   }
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Starting values for downward recurrence, Spherical Bessel func
void bess_rsv(/* FUNCTION */
      int *nmax,
      int *dig,
      double *x,
      double *RN,
      double *CN
   ){
//--------------------------------------
   double Rn[*dig], Cn[*dig];
   Rn[*dig]=0.0;  // n-th Bessel func
   Cn[*dig]=1.0;  // Derivative of the n-th Bessel func
   int n,j;
   for(n=*dig;n>0;n--){
      // Downward Recursion
      Rn[n-1]=lcfe_afs(n+*nmax,*x)*Rn[n  ]+Cn[n];
      Cn[n-1]=lcfe_afs(n+*nmax,*x)*Rn[n-1]-Rn[n];
      //Renormalization may be required
      if(fabs(Rn[n-1])>1e2){
         for(j=*dig;j>=n-1;j--){
            Cn[j]=Cn[j]/Rn[n-1];
            Rn[j]=Rn[j]/Rn[n-1];
         }
      }
   }
   *RN=Rn[0];
   *CN=Cn[0];
}
//------------------------------------------------------------------------------
// Array of Riccati-Bessel func
void bess_ric(/* FUNCTION */
      int *nmax, 
      double *x, 
      double *Rn, 
      double *Cn
   ){
//--------------------------------------
   double R0=0.0;  // Initialization values
   double R1=0.0;  // Initialization values
   double j1=0.0;  // Initialization values
   bess_szr(x,&R0);
   bess_sun(x,&R1);
   bess_sun(x,&j1);
   R0=*x*R0;
   R1=*x*R1;
   int dig=15;
   bess_rsv(nmax,&dig,x,&Rn[*nmax],&Cn[*nmax]);
   int n,j;
   for(n=*nmax;n>0;n--){
      // Downward Recursion
      Rn[n-1]=lcfe_afs(n,*x)*Rn[n  ]+Cn[n];
      Cn[n-1]=lcfe_afs(n,*x)*Rn[n-1]-Rn[n];
      //Renormalization may be required
      if(fabs(Rn[n-1])>1e2){
         for(j=*nmax;j>=n-1;j--){
            Cn[j]=Cn[j]/Rn[n-1];
            Rn[j]=Rn[j]/Rn[n-1];
         }
      }
   }
   // Normalization factor for the func
   double pf=1.0;
   if(fabs(R0)>1e-10){
      pf=(R0)/Rn[0];
   }else{
      pf=(R1)/Rn[1];
   }
   // Normalization factor for the derivative
   double pd=1.0;
   if(fabs(R1)>1e-10){
      pd=cos(*x)/Cn[0];
   }else{
      pd=(-j1+R0)/Cn[1];
   }
   // Normalizations
   for(n=0;n<=*nmax;n++){
      Rn[n]=Rn[n]*pf;
      Cn[n]=Cn[n]*pd;
   }
}
/*------------------------------------------------------------------------------
 *                  VECTOR SPHERICAL WAVE FUNC                                *
------------------------------------------------------------------------------*/
//------------------------------------------------------------------------------
// Psi_m(\vec{k},\vec{r}): Basic func for cylindrical simetries
//------------------------------------------------------------------------------
double complex vswf_psi(/* FUNCTION */
      int *m,int *s,
      double *gamma,double *kz,
      double *x,double *y,double *z
   ){
   //-----------------------------------
   double rho=sqrt(*x*(*x)+*y*(*y));
   double cph,sph;
   cph=*x/rho;
   sph=*y/rho;
   double complex eiph=cpow(cph+I*sph,*m*(*s));
   double complex eikz=cexp(I*(*kz)*(*z)); 
   // gsl
   double jn=gsl_sf_bessel_Jn(*m,rho*(*gamma));
   // math.h
   //double jn=jn(*m,rho*(*gamma));
   double complex u=jn*eiph*eikz;
   return(u);
}
//------------------------------------------------------------------------------
// RECTANGULAR WAVE GUIDE - TM MODE
//------------------------------------------------------------------------------
void drwg_tmm(/* FUNCTION */
      double *kx, double *ky, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){
   //-----------------------------------
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
   *Em=(ETMx-I*ETMy)/sqrt(2.0);
   *Ez=sin(*kx*(*x))*sin(*ky*(*y))*cexp(I*(*kz)*(*z));
   *Ep=(ETMx+I*ETMy)/sqrt(2.0);
   //
   double complex HTMx=-I*(k*(*ky)/gamma2)*skx*cky*cexp(I*(*kz)*(*z));
   double complex HTMy= I*(k*(*kx)/gamma2)*ckx*sky*cexp(I*(*kz)*(*z));
   //
   *Hm=(HTMx-I*HTMy)/sqrt(2.0);
   *Hz=0.0+I*0.0;
   *Hp=(HTMx+I*HTMy)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// RECTANGULAR WAVE GUIDE - TE MODE
//------------------------------------------------------------------------------
void drwg_tem(/* FUNCTION */
      double *kx, double *ky, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){
   //-----------------------------------
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
void vwfd_rwg(/* FUNCTION */
      int *TE,
      int *nx, int *ny, int *nz,
      double *kx, double *ky, double *kz,
      double *x,  double *y,  double *z,
      double *rx, double *ry, double *rz,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){ 
   //-----------------------------------
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           if(*TE==1){
              drwg_tem(kx,ky,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Hm[i],&Hz[i],&Hp[i],
                    &Em[i],&Ez[i],&Ep[i]);
           }
           if(*TE==0){
              drwg_tmm(kx,ky,kz,
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
void dcwg_all(/* FUNCTION */
      int *MD,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){
   //-----------------------------------
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double ctg=*kz/(*gamma);
   double csc=k/(*gamma);
   int msm=*M-*S;
   int msp=*M+*S;
   //
   *Em= vswf_psi(&msm,S,gamma,kz,x,y,z)*( I*(*S)*ctg)/sqrt(2.0);
   *Ez= vswf_psi(M,   S,gamma,kz,x,y,z);
   *Ep= vswf_psi(&msp,S,gamma,kz,x,y,z)*(-I*(*S)*ctg)/sqrt(2.0);
   // If TM, H = -H, E = E
   // If TE, E =  H, H = E
   *Hm=-(*MD)*vswf_psi(&msm,S,gamma,kz,x,y,z)*(*S*csc)/sqrt(2.0);
   *Hz=0;
   *Hp=-(*MD)*vswf_psi(&msp,S,gamma,kz,x,y,z)*(*S*csc)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void vwfd_cwg(/* FUNCTION */
      int *TE, int *m, int *s,
      int *nx, int *ny, int *nz,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double *rx, double *ry, double *rz,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){ 
   //-----------------------------------
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
              dcwg_all(&MD,m,s,
                    gamma,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Em[i],&Ez[i],&Ep[i],
                    &Hm[i],&Hz[i],&Hp[i]);
           }
           if(*TE==0){
              // MODO TM
              MD=-1;
              dcwg_all(&MD,m,s,
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
void fdef_bbz(/* FUNCTION */
      int *MD,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){
   //-----------------------------------
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double cth=*kz/k;
   double sth=(*gamma)/k;
   int msm=*M-*S;
   int msp=*M+*S;
   // TM MODE
   *Em=(*MD)*vswf_psi(&msm,S,gamma,kz,x,y,z)*( I*(*S)*cth*sth)/sqrt(2.0);
   *Ez=(*MD)*vswf_psi(M,   S,gamma,kz,x,y,z)*sth*sth;
   *Ep=(*MD)*vswf_psi(&msp,S,gamma,kz,x,y,z)*(-I*(*S)*cth*sth)/sqrt(2.0);
   //
   *Hm=vswf_psi(&msm,S,gamma,kz,x,y,z)*(*S*sth)/sqrt(2.0);
   *Hz=0;
   *Hp=vswf_psi(&msp,S,gamma,kz,x,y,z)*(*S*sth)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void vwfd_bbz(/* FUNCTION */
      int *TM, int *m, int *s,
      int *nx, int *ny, int *nz,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double *rx, double *ry, double *rz,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){ 
   //-----------------------------------
   int MD;
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           // MESMA ORDEM DO PARSER
           if(*TM==0){
              MD=1;
              fdef_bbz(&MD,m,s,
                    gamma,kz,
                    &x[ix],&y[iy],&z[iz],
                    &Hm[i],&Hz[i],&Hp[i],
                    &Em[i],&Ez[i],&Ep[i]);
           }
           if(*TM==1){
              MD=-1;
              fdef_bbz(&MD,m,s,
                    gamma,kz,&x[ix],
                    &y[iy],&z[iz],
                    &Em[i],&Ez[i],&Ep[i],
                    &Hm[i],&Hz[i],&Hp[i]);
           }
           i++; 
         }    
      }    
   }    
}
//------------------------------------------------------------------------------
// BESSEL BEAMS P ORIENTED
//------------------------------------------------------------------------------
void fdef_bbp(/* FUNCTION */
      int *P,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){
   //-----------------------------------
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double cth=*kz/k;
   double sth=(*gamma)/k;

   int pp1=1+*P;
   int pm1=1-*P;
   int mm=*M;
   int mspm1=*M+*S*(*P-1);
   int mspm0=*M+*S*(*P);
   int mspp1=*M+*S*(*P+1);

   *Em=pp1*vswf_psi(&mm,S,gamma,kz,x,y,z)/2
      -*P*(sth*sth)*vswf_psi(&mspm1,S,gamma,kz,x,y,z)/2;
   *Ez= -I*(*P)*(*S)*cth*sth*vswf_psi(&mspm0,S,gamma,kz,x,y,z)/sqrt(2.0);
   *Ep=pm1*vswf_psi(&mm,S,gamma,kz,x,y,z)/2
      +*P*(sth*sth)*vswf_psi(&mspp1,S,gamma,kz,x,y,z)/2;
   *Hm=-I*(*P)*cth*vswf_psi(&mm,S,gamma,kz,x,y,z)*pp1/2.0;
   *Hz=-(*S)*sth*vswf_psi(&mspm0,S,gamma,kz,x,y,z)/sqrt(2.0);
   *Hp=-I*(*P)*cth*vswf_psi(&mm,S,gamma,kz,x,y,z)*pm1/2.0;
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void vwfd_bbp(/* FUNCTION */
      int *p, int *m, int *s, 
      int *nx, int *ny, int *nz,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double *rx, double *ry, double *rz,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
   ){ 
   //-----------------------------------
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           fdef_bbp(p,m,s,
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
// VSWF - CALCULATIONS
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// CONSTANTS FOR Qlm (NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS)
//------------------------------------------------------------------------------
double alfaQ(/* FUNCTION */
      int lo, 
      int mo
   ){
   //-----------------------------------
   double l=lo;
   double m=mo;
   return sqrt(((2*l-1)*(2*l+1))/((l-m)*(l+m)));
}
//------------------------------------------------------------------------------
double betaQ(/* FUNCTION */
      int lo, 
      int mo
   ){
   //-----------------------------------
   double l=lo;
   double m=mo;
   return sqrt((2*l+1)/(2*l-3))*sqrt(((l+m-1)*(l-m-1))/((l-m)*(l+m)));
}
//------------------------------------------------------------------------------
double gammaQ(/* FUNCTION */
      int lo
   ){
   //-----------------------------------
   double l=lo;
   return sqrt((2*l+1)/(2*l));
}
//------------------------------------------------------------------------------
double deltaQ(/* FUNCTION */
      int lo
   ){
   //-----------------------------------
   double l=lo;
   return sqrt(2*l+1);
}
//------------------------------------------------------------------------------
// POSITIONER IN THE (l,m) LIST
//------------------------------------------------------------------------------
int jlm(/* FUNCTION */
      int l, 
      int m
   ){
   //-----------------------------------
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
* i=iz+nz*(iy-1)+nz*ny*(ix-1) : comecando de 1
* i=iz+nz*iy+nz*ny*ix         : comecando de 0
*/
int lijk(/* FUNCTION */
      int ix, 
      int ni, 
      int iy,
      int ny, 
      int iz, 
      int nz
   ){
   //-----------------------------------
   return(iz+nz*iy+nz*ny*ix);
}
//------------------------------------------------------------------------------
// CONSTANTS FOR \vec{X}_{lm}
//------------------------------------------------------------------------------
double clmp(/* FUNCTION */
      int l, 
      int m,
      int p
   ){
   return sqrt(1.*l*(l+1.)-1.*m*(m+p));
}
//------------------------------------------------------------------------------
double Klmqp(/* FUNCTION */
      int l,
      int m,
      int q,
      int p
   ){
   return ((1.*l+p*m)*(1.*l+p*q))/((2.*l-1.)*(2.*l+1.));
}
//------------------------------------------------------------------------------
// NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
//------------------------------------------------------------------------------
void vswf_qlm(/* FUNCTION */
      double *x,
      int *lmax, 
      double *Qlm, double *dQlm,
      int *lo, int *mo
   ){
   //-----------------------------------
   int l,m;
   double kl;
   double cth = *x;
   double sth = sqrt(1-cth*cth);
   double cs2 = pow(sth,-2);
   //-----------------------------------
   // Qlm - First 4 terms
   Qlm[jlm(0, 0)]=1/sqrt(4*M_PI);
   Qlm[jlm(1, 1)]=-gammaQ(1)*sth*Qlm[jlm(0,0)]; // Q11
   Qlm[jlm(1, 0)]=sqrt(3.0)*cth*Qlm[jlm(0,0)];  // Q10
   Qlm[jlm(1,-1)]=-Qlm[jlm(1,1)];               // Q11*(-1)
   //-----------------------------------
   // dQlm - First 4 terms
   dQlm[jlm(0, 0)]=0.0;
   dQlm[jlm(1, 1)]= -cth*Qlm[jlm(1, 1)]*cs2;
   dQlm[jlm(1, 0)]=-(cth*Qlm[jlm(1, 0)]-sqrt(3.)*Qlm[jlm(0,0)])*cs2;
   dQlm[jlm(1,-1)]= -cth*Qlm[jlm(1,-1)]*cs2;
   //-----------------------------------
   // lo - First 4 terms
   lo[jlm(0, 0)] = 0;
   lo[jlm(1,-1)] = 1;
   lo[jlm(1, 0)] = 1;
   lo[jlm(1, 1)] = 1;
   // mo - First 4 terms
   mo[jlm(0, 0)] =  0;
   mo[jlm(1,-1)] = -1;
   mo[jlm(1, 0)] =  0;
   mo[jlm(1, 1)] =  1;
   //-----------------------------------
   for(l=2;l<=(*lmax);l++){
      //Qlm positives extremes 
      Qlm[jlm(l, l  )]=-gammaQ(l)*sth*Qlm[jlm(l-1,l-1)];
      Qlm[jlm(l, l-1)]= deltaQ(l)*cth*Qlm[jlm(l-1,l-1)];
      //dQlm positives extremes 
      Qlm[jlm(l,-l  )]=pow(-1,l  )*Qlm[jlm(l,l  )];
      Qlm[jlm(l,-l+1)]=pow(-1,l-1)*Qlm[jlm(l,l-1)];
      
      kl=1/(sqrt(l*(l+1)));
      //--------------------------------
      for(m=l; m>=0; m--){
         lo[jlm(l, m)] =  l;
         lo[jlm(l,-m)] =  l;
         mo[jlm(l, m)] =  m;
         mo[jlm(l,-m)] = -m;

         if(m<l-1){
            Qlm[jlm(l, m)]=alfaQ(l,m)*cth*Qlm[jlm(l-1,m)]
               -betaQ(l,m)*Qlm[jlm(l-2,m)];
            Qlm[jlm(l,-m)]=pow(-1,m)*Qlm[jlm(l, m)];
         }

         dQlm[jlm(l, m)] = -l*cth*cs2*Qlm[jlm(l, m)]+
             sqrt((2*l+1.)/(2*l-1.))*sqrt((l+m)*(l-m))*cs2*Qlm[jlm(l-1, m)];

         dQlm[jlm(l,-m)] = -l*cth*cs2*Qlm[jlm(l,-m)]+
             sqrt((2*l+1.)/(2*l-1.))*sqrt((l+m)*(l-m))*cs2*Qlm[jlm(l-1,-m)];
      }
   }
}
//------------------------------------------------------------------------------
// POSITION DEPENDENT CALCULATIONS - MULTIPLE LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void vswf_def(/* FUNCTION */
      double *k,double *x, double *y, double *z,
      int *lmax, 
      double complex *GTE, double complex *GTM,
      double complex *Em, double complex *Ez, double complex *Ep,
      double complex *Hm, double complex *Hz, double complex *Hp
   ){
   //-----------------------------------
   int l;
//   int LMAX=*lmax*(*lmax+2); # MUST CALCULATE ONE MORE POINT #
   int LMAXE=(*lmax+1)*(*lmax+3);
   double cph; 
   double sph;
   double rho=sqrt(*x*(*x)+*y*(*y));
   double r=sqrt(rho*rho+*z*(*z));
   double sth=rho/r;
   double cth=*z/r;
   //-----------------------------------
   if((*x==0)&&(*y==0)){
      cph=1;
      sph=0;
   }else{
      cph=*x/rho;
      sph=*y/rho;
   }
   /*---------------------------------*/
   /* SPHERICAL BESSEL                */
   /*---------------------------------*/
   // GSL
   //gsl_sf_bessel_jl_steed_array(*lmax+1,*k*r,JLM);
   // INTERNAL
   double  JLM[*lmax+2]; //Spherical Bessel func.
   double dJLM[*lmax+2]; //Derivative of Spherical Bessel func
   int lp1=*lmax+1;
   double kr =*k*r;
   int maxn=2000;
   bess_sph(&lp1,&kr,(void *)&JLM,(void *)&dJLM,&maxn);
   //-----------------------------------
   // Qlm - First 4 terms
   double Qlm[LMAXE];
   Qlm[jlm(0, 0)]=1/sqrt(4*M_PI);
   Qlm[jlm(1, 1)]=-gammaQ(1)*sth*Qlm[jlm(0,0)]; // Q11
   Qlm[jlm(1, 0)]=sqrt(3.0)*cth*Qlm[jlm(0,0)];  // Q10
   Qlm[jlm(1,-1)]=-Qlm[jlm(1,1)];               // Q11*(-1)
   //-----------------------------------
   // Complex Exponencial for m=-1,0,1
   double complex Eim[2*(*lmax)+3];
   Eim[*lmax-1]=(cph-I*sph);
   Eim[*lmax  ]=1+I*0;
   Eim[*lmax+1]=(cph+I*sph);
   //-----------------------------------
   // Ylm - First 4 terms
   double complex Ylm[LMAXE];
   Ylm[jlm(0, 0)]=Qlm[jlm(0, 0)];
   Ylm[jlm(1,-1)]=Qlm[jlm(1,-1)]*Eim[*lmax-1];
   Ylm[jlm(1, 0)]=Qlm[jlm(1, 0)];
   Ylm[jlm(1, 1)]=Qlm[jlm(1, 1)]*Eim[*lmax+1];
   /*---------------------------------*/
   /* VECTOR SPHERICAL HARMONICS      */
   /*---------------------------------*/
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
   /*---------------------------------*/
   /* HANSEN MULTIPOLES               */
   /*---------------------------------*/
   // M
   double complex MM;
   double complex MZ;
   double complex MP;
   // N
   double complex NM;
   double complex NZ;
   double complex NP;
   /*---------------------------------*/
   /* OTHERS                          */
   /*---------------------------------*/
   double kl;
   //-----------------------------------
   // MAIN LOOP
   //-----------------------------------
   for(l=1;l<=(*lmax);l++){
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
      //--------------------------------
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
         XM=kl*clmp(l,m, 1)*Ylm[jlm(l,m-1)]/sqrt(2);
         XZ=kl*m*Ylm[jlm(l,m  )];
         XP=kl*clmp(l,m,-1)*Ylm[jlm(l,m+1)]/sqrt(2);
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
void vswf_pwe(/* FUNCTION */
      double *k,double *x, double *y, double *z,
      int *lmax, int *nx, int *ny, int *nz, 
      double complex *GTE, double complex *GTM,
      double *rx, double *ry, double *rz,
      double complex *Em, double complex *Ez, double complex *Ep,
      double complex *Hm, double complex *Hz, double complex *Hp
   ){ 
   //-----------------------------------
   int i=0;
   #pragma omp parallel for
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           /*In general, the offset formula for an array with dimensions 
            * [d1 d2 d3 ... dn] using any subscripts (s1 s2 s3 ... sn) is
            * (sn-1)(dn-1)(dn-2)...(d1)+(sn-1-1)(dn-2)...(d1)+...+(s2-1)(d1)+s1
            * Because of this scheme, you can index an array 
            * using any number of subscripts. 
            * 
            * i=iz+nz*(iy-1)+nz*ny*(ix-1) : comecando de 1
            * i=iz+nz*iy+nz*ny*ix         : comecando de 0
            */
           i=iz+(*nz)*iy+(*nz)*(*ny)*ix; 
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           vswf_def(k,
                 &x[ix],&y[iy],&z[iz],
                 lmax,GTE,GTM,
                 &Em[i],&Ez[i],&Ep[i],
                 &Hm[i],&Hz[i],&Hp[i]);
         } /* end for iz */    
      } /* end for iy */    
   } /* end for ix */    
}
//------------------------------------------------------------------------------
void vswf_vsh(/* FUNCTION */
      double *x, double *y, double *z,
      int *lmax, 
      double complex *Y,
      double complex *VM, double complex *VZ, double complex *VP,
      double complex *XM, double complex *XZ, double complex *XP,
      double complex *YM, double complex *YZ, double complex *YP
   ){
   //-----------------------------------
//   int LMAX=*lmax*(*lmax+2); # MUST CALCULATE ONE MORE POINT #
   int LMAXE=(*lmax+1)*(*lmax+3);
   int l,m;
   double cph; 
   double sph;
   double rho=sqrt(*x*(*x)+*y*(*y));
   double r=sqrt(rho*rho+*z*(*z));
   double sth=rho/r;
   double cth=*z/r;
   //-----------------------------------
   if((*x==0)&&(*y==0)){
      cph=1;
      sph=0;
   }else{
      cph=*x/rho;
      sph=*y/rho;
   }
   //-----------------------------------
   // Qlm - First 4 terms
   double Qlm[LMAXE];
   Qlm[jlm(0, 0)]=1/sqrt(4*M_PI);
   Qlm[jlm(1, 1)]=-gammaQ(1)*sth*Qlm[jlm(0,0)]; // Q11
   Qlm[jlm(1, 0)]=sqrt(3.0)*cth*Qlm[jlm(0,0)];  // Q10
   Qlm[jlm(1,-1)]=-Qlm[jlm(1,1)];               // Q11*(-1)
   //-----------------------------------
   // Complex Exponencial for m=-1,0,1
   double complex Eim[2*(*lmax)+3];
   Eim[*lmax-1]=(cph-I*sph);
   Eim[*lmax  ]=1+I*0;
   Eim[*lmax+1]=(cph+I*sph);
   //-----------------------------------
   // Ylm - First 4 terms
   double complex Ylm[LMAXE];
   Ylm[jlm(0, 0)]=Qlm[jlm(0, 0)];
   Ylm[jlm(1,-1)]=Qlm[jlm(1,-1)]*Eim[*lmax-1];
   Ylm[jlm(1, 0)]=Qlm[jlm(1, 0)];
   Ylm[jlm(1, 1)]=Qlm[jlm(1, 1)]*Eim[*lmax+1];
   /*---------------------------------*/
   /* VECTOR SPHERICAL HARMONICS      */
   /*---------------------------------*/
   // r
   double complex rm=sth*(cph-I*sph)/sqrt(2);
   double complex rz=cth;
   double complex rp=sth*(cph+I*sph)/sqrt(2);
   double kl;
   //-----------------------------------
   // MAIN LOOP
   //-----------------------------------
   Y[jlm(0,0)]=Ylm[jlm(0, 0)];

   YM[jlm(0,0)]=rm*Ylm[jlm(0,0)];
   YZ[jlm(0,0)]=rz*Ylm[jlm(0,0)];
   YP[jlm(0,0)]=rp*Ylm[jlm(0,0)];
   for(l=1;l<=(*lmax);l++){
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
      kl=1./(sqrt(1.*l*(l+1.)));
      //--------------------------------
      for(m=l; m>=(-l); m--){
         // CALCULATIONS OF SSH
         if(m>=0){
            Qlm[jlm(l+1, m)]=alfaQ(l+1,m)*cth*Qlm[jlm(l,m)]
               -betaQ(l+1,m)*Qlm[jlm(l-1,m)];
            Qlm[jlm(l+1,-m)]=pow(-1,m)*Qlm[jlm(l+1, m)];
            Ylm[jlm(l+1, m)]=Qlm[jlm(l+1, m)]*Eim[*lmax+m];
            Ylm[jlm(l+1,-m)]=Qlm[jlm(l+1,-m)]*Eim[*lmax-m];
         }
         // CALCULATIONS OF VSH
         Y[jlm(l,m)]=Ylm[jlm(l,m)];
         // X
         XM[jlm(l,m)]=kl*clmp(l,m,-1)*Ylm[jlm(l,m-1)]/sqrt(2);
         XZ[jlm(l,m)]=kl*m*Ylm[jlm(l,m  )];
         XP[jlm(l,m)]=kl*clmp(l,m, 1)*Ylm[jlm(l,m+1)]/sqrt(2);
         // Y
         YM[jlm(l,m)]=rm*Ylm[jlm(l,m)];
         YZ[jlm(l,m)]=rz*Ylm[jlm(l,m)];
         YP[jlm(l,m)]=rp*Ylm[jlm(l,m)];
         // V
         VM[jlm(l,m)]=rm*XZ[jlm(l,m)]-rz*XM[jlm(l,m)];
         VZ[jlm(l,m)]=rp*XM[jlm(l,m)]-rm*XP[jlm(l,m)];
         VP[jlm(l,m)]=rz*XP[jlm(l,m)]-rp*XZ[jlm(l,m)];
      }
   }
}
//------------------------------------------------------------------------------
// MIE COEFFICIENTS CALCULATION
//------------------------------------------------------------------------------
void lmie_rho(/* FUNCTION */
      double complex  *m, double complex *x, 
      int *lmax, int *NMAX,
      double complex *an, double complex *bn
   ){
   // initial values
   double complex p0 = csin(*x);
   double complex q0 = ccos(*x);
   double complex p1 = p0/(*x)-q0;
   double complex q1 = q0/(*x)+p0;
   double complex z0 = p0+I*q0;
   double complex z1 = p1+I*q1;
   // starting series
   double complex g,Cn;
   g  = z0/z1;
   Cn = p0/z0;
   //------------------------------------
   // CALCULATION OF THE RATIOS
   double complex mx = *m * *x;
   double complex rho1[*lmax+1],rhom[*lmax+1],Dn1[*lmax+1],Dnm[*lmax+1];
   lcfc_ric(lmax,  x,rho1,Dn1,NMAX);
   lcfc_ric(lmax,&mx,rhom,Dnm,NMAX);
   //------------------------------------
   // UPWARD RECURRENCE
   int l;
   double complex k, Ta, Tb;

   for(l=0; l<*lmax; l++){
      Cn=Cn*g/rho1[l];

      k = (1-1/((*m)*(*m)))*(l+1)/(*x);
      Ta = (rhom[l]/(*m)-rho1[l]+k)/(rhom[l]/(*m)-g+k);
      Tb = (rhom[l]*(*m)-rho1[l]  )/(rhom[l]*(*m)-g  );
      an[l] = Cn*Ta; 
      bn[l] = Cn*Tb;

      g=1/(lcfc_afs(2*(l+1)+1,*x)-g);
   }
}
//------------------------------------------------------------------------------
void lmie_log(/* FUNCTION */
      double complex  *m, double complex *x, 
      int *lmax, int *NMAX,
      double complex *an, double complex *bn
   ){
   //------------------------------------
   // CALCULATION OF THE RATIOS
   double complex mx = *m * *x;
   double complex rho1[*lmax+1],rhom[*lmax+1],An1[*lmax+1],Anm[*lmax+1];
   lcfc_ric(lmax,  x,rho1,An1,NMAX);
   lcfc_ric(lmax,&mx,rhom,Anm,NMAX);
   double complex Cn,Bn,Ta,Tb;
   // STARTING VALUES
   Cn = 1/(1+I*(ccos(*x)+*x*csin(*x))/(csin(*x)-*x*ccos(*x)));
   Bn = -lcfc_afs(1,*x)+1/(lcfc_afs(1,*x)+I);

   Ta = (Anm[1]/(*m)-An1[1])/(Anm[1]/(*m)-Bn);
   Tb = (Anm[1]*(*m)-An1[1])/(Anm[1]*(*m)-Bn);

   an[0] = Cn*Ta;
   bn[0] = Cn*Tb;
   
   // UPWARD RECURRENCE
   int n;
   for(n=1;n<*lmax;n++){
      Bn = -lcfc_afs(n+1,*x)+1/(lcfc_afs(n+1,*x)-Bn);
      Cn = Cn*(Bn+(n+1)/(*x))/(An1[n+1]+(n+1)/(*x));

      Ta = (Anm[n+1]/(*m)-An1[n+1])/(Anm[n+1]/(*m)-Bn);
      Tb = (Anm[n+1]*(*m)-An1[n+1])/(Anm[n+1]*(*m)-Bn);
      
      an[n] = Cn*Ta;
      bn[n] = Cn*Tb;
   }
}
/*----------------------------------------------------------------------------*/
/*                    MIE OPTICAL FORCE CALCULATIONS                          */
/*----------------------------------------------------------------------------*/
void lmie_ofc(/* FUNCTION */
      double complex *M,double complex *x, 
      int *lmax,int *NMAX,
      double complex *GTE, double complex *GTM,
      double *Fx, double *Fy, double *Fz
   ){
   double complex an[*lmax+1],bn[*lmax+1],An,Bn,Cn;
   int lmaxe = *lmax+1;
   lmie_log(M,x,&lmaxe,NMAX,an,bn);

   double K1,K2,K3,K4;
   double complex U,V,W,Y,Fxy,Fzz;

   Fxy=0.0+0.0*I;
   Fzz=0.0+0.0*I;

   int l,m;
   int ll;
   for(l=0;l<*lmax;l++){

      An = an[l]+conj(an[l+1])-2*an[l]*conj(an[l+1]);
      Bn = bn[l]+conj(bn[l+1])-2*bn[l]*conj(bn[l+1]);
      Cn = an[l]+conj(bn[l  ])-2*an[l]*conj(bn[l  ]);

      ll=l+1;

      for(m=-ll;m<=ll;m++){

         K1 = sqrt((ll*(ll+2.))/((2*ll+1.)*(2*ll+3.)));
         K2 = sqrt((ll+m+2.)*(ll+m+1.));
         K3 = sqrt((ll-m)*(ll+m+1.));
         K4 = sqrt((ll+m+1.)*(ll-m+1.));

         U = An*GTM[jlm(ll,m)]*conj(GTM[jlm(ll+1,m+1)])
           + Bn*GTE[jlm(ll,m)]*conj(GTE[jlm(ll+1,m+1)])
           + conj(An*GTM[jlm(ll,-m)])*GTM[jlm(ll+1,-m-1)]
           + conj(Bn*GTE[jlm(ll,-m)])*GTE[jlm(ll+1,-m-1)];

         V = Cn*GTM[jlm(ll,m)]*conj(GTE[jlm(ll,m+1)])
           -  conj(Cn)*GTE[jlm(ll,m)]*conj(GTM[jlm(ll,m+1)]);

         W = An*GTM[jlm(ll,m)]*conj(GTM[jlm(ll+1,m)])
           - Bn*GTE[jlm(ll,m)]*conj(GTE[jlm(ll+1,m)]);

         Y = Cn*GTM[jlm(ll,m)]*conj(GTE[jlm(ll,m)]);

         Fxy = (K1*K2*U-K3*V/(1.*ll))/(ll+1.);
         Fzz = (K1*K4*W-m*Y/(1.*ll))/(ll+1.);

         *Fx = *Fx + creal(Fxy);
         *Fy = *Fy + cimag(Fxy);
         *Fz = *Fz + creal(I*Fzz);

        // printf("Fx = %f, Fy = %f, Fz = %f\n",*Fx,*Fy,*Fz);
      }
   }
}
/*----------------------------------------------------------------------------*/
void vswf_gpw(/* FUNCTION */ 
      double *kx, double *ky, double *kz,
      double complex *ux, double complex *uy, double complex *uz,
      int *lmax,
      double complex *GTE, double complex *GTM
   ){
   int LMAX=(*lmax)*(*lmax+2)+1;
   double k = sqrt(pow(*kx,2)+pow(*ky,2)+pow(*kz,2));
   double u = sqrt(conj(*ux)**ux+conj(*uy)**uy+conj(*uz)**uz);
   

   double complex Y[LMAX];
   double complex VM[LMAX],VZ[LMAX],VP[LMAX];
   double complex XM[LMAX],XZ[LMAX],XP[LMAX];
   double complex YM[LMAX],YZ[LMAX],YP[LMAX]; 

   vswf_vsh(kx,ky,kz,lmax,Y,VM,VZ,VP,XM,XZ,XP,YM,YZ,YP);
   // k normalization
   double hkx = *kx/k;
   double hky = *ky/k;
   double hkz = *kz/k;
   double complex hkp = (hkx+I*hky)/sqrt(2.);
   double complex hkm = (hkx-I*hky)/sqrt(2.);
   // u normalization
   double complex hux = *ux/u;
   double complex huy = *uy/u;
   double complex huz = *uz/u;
   double complex hup = (hux+I*huy)/sqrt(2.);
   double complex hum = (hux-I*huy)/sqrt(2.);
   // ALFA a = k x u
   double complex ham = I*(hkm*huz-hkz*hum);
   double complex haz = I*(hkp*hum-hkm*hup);
   double complex hap = I*(hkz*hup-hkp*huz);

   int l,m;
   for(l=1;l<=*lmax;l++){
      for(m=-l;m<=l;m++){
         GTE[jlm(l,m)] = 4*M_PI*cpow(I,l)*(conj(XM[jlm(l,m)])*hum+conj(XZ[jlm(l,m)])*huz+conj(XP[jlm(l,m)])*hup);
         GTM[jlm(l,m)] = 4*M_PI*cpow(I,l)*(conj(XM[jlm(l,m)])*ham+conj(XZ[jlm(l,m)])*haz+conj(XP[jlm(l,m)])*hap);
      }
   }
}
//------------------------------------------------------------------------------
void vswf_gwg(/* FUNCTION */ 
      double *gamma, double *kz,
      int *lmax,
      double complex *A, double complex *B
   ){
   int LMAX=(*lmax)*(*lmax+2)+1;

   double k = sqrt(pow(*gamma,2)+pow(*kz,2));
   double cth = (*kz)/k;

   double Qlm[LMAX],dQlm[LMAX];
   int ll[LMAX],mm[LMAX];

   vswf_qlm(&cth,lmax,Qlm,dQlm,ll,mm);

   int l,m;
   double llp1;
   for(l=1;l<=*lmax;l++){
      for(m=-l;m<=l;m++){
         llp1 = 1/sqrt((1.*l)*(1.*l+1.));
         A[jlm(l,m)] = 2*cpow(I,l)*pow(k/(*gamma),2)*Qlm[jlm(l,m)]*m*llp1;
         B[jlm(l,m)] = 2*cpow(I,l-1)*dQlm[jlm(l,m)]*llp1;
      }
   }
}
//------------------------------------------------------------------------------
void vswf_cwg(/* FUNCTION */ 
      int *M,int *s,int *lmax,int *NMAX,int *TM,
      double *gamma, double *kz,
      double *x, double *y, double *z,
      double complex *GTE, double complex *GTM
   ){
   int LMAX=(*lmax)*(*lmax+2)+1;

   double rho = sqrt(pow(*x,2)+pow(*y,2));
   double k = sqrt(pow(*gamma,2)+pow(*kz,2));
   double gr = *gamma*rho;

   double cth = (*kz)/k;
   double cph = *x/rho;
   double sph = *y/rho;
   double zz = (*kz)*(*z);

   double Qlm[LMAX],dQlm[LMAX],Jn[LMAX],Dn[LMAX];
   int ll[LMAX],mm[LMAX];

   int ls,le,mo;
   ls = abs(*M-*lmax);
   le = abs(*M+*lmax);
   if(ls>le){
      le=ls;
      ls=0;
   }
   vswf_qlm(&cth,lmax,Qlm,dQlm,ll,mm);
   bess_cyl(&le,&gr,Jn,Dn,NMAX); 

   double complex A,B,A0,psi;
   double llp1;
   int l,m;
   for(l=1;l<=*lmax;l++){
      for(m=-l;m<=l;m++){
         mo = *M-*s*m;
         A0 = 2*M_PI*cpow(-I*(*s),m);
         if(mo<0){
            psi = Jn[abs(mo)]*cpow((cph+*s*I*sph),mo)*cexp(I*zz)*pow(-1,abs(mo));
         }else{
            psi = Jn[abs(mo)]*cpow((cph+*s*I*sph),mo)*cexp(I*zz);
         }

         llp1 = 1/sqrt((1.*l)*(1.*l+1.));
         A = 2*cpow(I,l)*pow(k/(*gamma),2)*Qlm[jlm(l,m)]*m*llp1;
         B = 2*cpow(I,l-1)*dQlm[jlm(l,m)]*llp1;

         if(*TM==1){
            GTE[jlm(l,m)] =  A*A0*psi;
            GTM[jlm(l,m)] = -B*A0*psi;
         }else{
            GTE[jlm(l,m)] =  B*A0*psi;
            GTM[jlm(l,m)] =  A*A0*psi;
         }
      }
   }
}
//------------------------------------------------------------------------------
void vswf_bbp(/* FUNCTION */ 
      int *M,int *S,int *P,int *lmax,int *NMAX,
      double *gamma, double *kz,
      double *x, double *y, double *z,
      double complex *GTE, double complex *GTM
   ){
   int LMAX=(*lmax)*(*lmax+2)+1;

   double rho = sqrt(pow(*x,2)+pow(*y,2));
   double k = sqrt(pow(*gamma,2)+pow(*kz,2));
   double gr = *gamma*rho;

   double cth = (*kz)/k;
   double sth = (*gamma)/k;
   double cph = *x/rho;
   double sph = *y/rho;
   double zz = (*kz)*(*z);

   double Qlm[LMAX],dQlm[LMAX],Jn[LMAX],Dn[LMAX];
   int ll[LMAX],mm[LMAX];

   int ls,le,mo;
   ls = abs(*M-*lmax-1);
   le = abs(*M+*lmax+1);
   if(ls>le){
      le=ls;
      ls=0;
   }
   vswf_qlm(&cth,lmax,Qlm,dQlm,ll,mm);
   bess_cyl(&le,&gr,Jn,Dn,NMAX); 

   double complex A0,psi;
   int l,m;
   for(l=1;l<=*lmax;l++){
      for(m=-l;m<=l;m++){
         mo = *M-*S*(m-*P);
         A0 = 4*M_PI*cpow(I,l-*S*(m-*P))/sqrt(2.*l*(l+1.));
         if(mo<0){
            psi = Jn[abs(mo)]*cpow((cph+*S*I*sph),mo)*cexp(I*zz)*pow(-1,abs(mo));
         }else{
            psi = Jn[abs(mo)]*cpow((cph+*S*I*sph),mo)*cexp(I*zz);
         }
         GTE[jlm(l,m)] =   A0*psi*clmp(l,m,-*P)*Qlm[jlm(l,m-*P)];
         GTM[jlm(l,m)] =-I*A0*psi*(cth*sth*dQlm[jlm(l,m)]-*P*m*Qlm[jlm(l,m)]/sth);
      }
   }
}
//------------------------------------------------------------------------------
void vswf_bbz(/* FUNCTION */ 
      int *M,int *S,int *lmax,int *NMAX,int *TM,
      double *gamma, double *kz,
      double *x, double *y, double *z,
      double complex *GTE, double complex *GTM
   ){
   int LMAX=(*lmax)*(*lmax+2)+1;

   double rho = sqrt(pow(*x,2)+pow(*y,2));
   double k = sqrt(pow(*gamma,2)+pow(*kz,2));
   double gr = *gamma*rho;

   double cth = (*kz)/k;
   double sth = (*gamma)/k;
   double cph = *x/rho;
   double sph = *y/rho;
   double zz = (*kz)*(*z);

   double Qlm[LMAX],dQlm[LMAX],Jn[LMAX],Dn[LMAX];
   int ll[LMAX],mm[LMAX];

   int ls,le,mo;
   ls = abs(*M-*lmax-1);
   le = abs(*M+*lmax+1);
   if(ls>le){
      le=ls;
      ls=0;
   }
   vswf_qlm(&cth,lmax,Qlm,dQlm,ll,mm);
   bess_cyl(&le,&gr,Jn,Dn,NMAX); 

   double complex A0,psi;
   int l,m;
   for(l=1;l<=*lmax;l++){
      for(m=-l;m<=l;m++){
         mo = *M-*S*m;
         A0 = 4*M_PI*cpow(I,l-*S*m)/sqrt(l*(l+1.));
         if(mo<0){
            psi = Jn[abs(mo)]*cpow((cph+*S*I*sph),mo)*cexp(I*zz)*pow(-1,abs(mo));
         }else{
            psi = Jn[abs(mo)]*cpow((cph+*S*I*sph),mo)*cexp(I*zz);
         }
         if(*TM==1){
            GTE[jlm(l,m)] =   A0*psi*m*Qlm[jlm(l,m)];
            GTM[jlm(l,m)] = I*A0*psi*pow(sth,2)*dQlm[jlm(l,m)];
         }else{
            GTE[jlm(l,m)] = I*A0*psi*pow(sth,2)*dQlm[jlm(l,m)];
            GTM[jlm(l,m)] =  -A0*psi*m*Qlm[jlm(l,m)];

         }
      }
   }
}
//------------------------------------------------------------------------------
void vswf_rwg(/* FUNCTION */ 
      int *TM,int *lmax,
      double *kx, double *ky, double *kz,
      double  *x, double  *y, double  *z,
      double complex *GTE, double complex *GTM
   ){
   int LMAX=(*lmax)*(*lmax+2)+1;

   double k = sqrt(pow(*kx,2)+pow(*ky,2)+pow(*kz,2));
   double gamma = sqrt(pow(*kx,2)+pow(*ky,2));

   double cth = (*kz)/k;
   double kzz = (*kz)*(*z);

   double Qlm[LMAX],dQlm[LMAX];
   int ll[LMAX],mm[LMAX];

   vswf_qlm(&cth,lmax,Qlm,dQlm,ll,mm);

   double complex A,B,eimz,f,g;
   double llp1;
   int l,m;
   int S = 1-2*(*TM);

   double complex EXmY = cexp(I*(*kx*(*x)-*ky*(*y)));
   double complex EXpY = cexp(I*(*kx*(*x)+*ky*(*y)));
   double czt  = *kx/gamma;
   double szt  = *ky/gamma;
   

//   #pragma omp parallel for private (m)
   for(l=1;l<=*lmax;l++){
      for(m=-l;m<=l;m++){

         llp1 = 1/sqrt((1.*l)*(1.*l+1.));
         A = 2*cpow(I,l)*pow(k/(gamma),2)*Qlm[jlm(l,m)]*m*llp1;
         B = 2*cpow(I,l-1)*dQlm[jlm(l,m)]*llp1;

         eimz = cpow(czt+I*szt,m);
         f = pow(-1,m);

         g = .5*M_PI*cexp(I*kzz)*((EXmY+f*conj(EXmY))*eimz+S*((EXpY+f*conj(EXpY))*conj(eimz)));
         
         if(*TM==1){
            GTE[jlm(l,m)] =  A*g;
            GTM[jlm(l,m)] = -B*g;
         }else{
            GTE[jlm(l,m)] =  B*g;
            GTM[jlm(l,m)] =  A*g;
         }
      }
   }
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void bscf_rwg(/* FUNCTION */
      int *TM, int *lmax, int *NMAX,
      double *kx, double *ky, double *kz,
      double complex *M, double complex *X,
      double *x, double *y, double *z,
      int *nx, int *ny, int *nz,
      double *rx, double *ry, double *rz,
      double *Fx, double *Fy, double *Fz
   ){
   //-----------------------------------
   int lmaxe=*lmax+1;
   int LMAXE=(lmaxe)*(lmaxe+2)+1;
   #pragma omp parallel for
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
            double complex GTE[LMAXE],GTM[LMAXE];
            int i=iz+(*nz)*iy+(*nz)*(*ny)*ix;
            rx[i]=x[ix];
            ry[i]=y[iy];
            rz[i]=z[iz];
            vswf_rwg(TM,&lmaxe,kx,ky,kz,&x[ix],&y[iy],&z[iz],GTE,GTM);
            lmie_ofc(M,X,lmax,NMAX,GTE,GTM,&Fx[i],&Fy[i],&Fz[i]);
         }
      }
   }
}
//------------------------------------------------------------------------------
void bscf_cwg(/* FUNCTION */
      int *TM, int *m, int *s,int *lmax, int *NMAX,
      double *gamma, double *kz,
      double complex *M, double complex *X,
      double *x, double *y, double *z,
      int *nx, int *ny, int *nz,
      double *rx, double *ry, double *rz,
      double *Fx, double *Fy, double *Fz
   ){
   //-----------------------------------
   int lmaxe=*lmax+1;
   int LMAXE=(lmaxe)*(lmaxe+2)+1;
   #pragma omp parallel for
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
            double complex GTE[LMAXE],GTM[LMAXE];
            int i=iz+(*nz)*iy+(*nz)*(*ny)*ix;
            rx[i]=x[ix];
            ry[i]=y[iy];
            rz[i]=z[iz];
            vswf_cwg(m,s,&lmaxe,NMAX,TM,gamma,kz,&x[ix],&y[iy],&z[iz],GTE,GTM);
            lmie_ofc(M,X,lmax,NMAX,GTE,GTM,&Fx[i],&Fy[i],&Fz[i]);
         } /* end for iz */    
      } /* end for iy */    
   } /* end for ix */    
}
//------------------------------------------------------------------------------
void bscf_bbz(/* FUNCTION */
      int *TM, int *m, int *s, int *lmax, int *NMAX,
      double *gamma, double *kz,
      double complex *M, double complex *X,
      double *x, double *y, double *z,
      int *nx, int *ny, int *nz,
      double *rx, double *ry, double *rz,
      double *Fx, double *Fy, double *Fz
   ){
   //-----------------------------------
   int lmaxe=*lmax+1;
   int LMAXE=(lmaxe)*(lmaxe+2)+1;
   #pragma omp parallel for
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
            double complex GTE[LMAXE],GTM[LMAXE];
            int i=iz+(*nz)*iy+(*nz)*(*ny)*ix;
            rx[i]=x[ix];
            ry[i]=y[iy];
            rz[i]=z[iz];
            vswf_bbz(m,s,&lmaxe,NMAX,TM,gamma,kz,&x[ix],&y[iy],&z[iz],GTE,GTM);
            lmie_ofc(M,X,lmax,NMAX,GTE,GTM,&Fx[i],&Fy[i],&Fz[i]);
         } /* end for iz */    
      } /* end for iy */    
   } /* end for ix */    
}
//------------------------------------------------------------------------------
void bscf_bbp(/* FUNCTION */
      int *TM, int *m, int *s, int *p, int *lmax, int *NMAX,
      double *gamma, double *kz,
      double complex *M, double complex *X,
      double *x, double *y, double *z,
      int *nx, int *ny, int *nz,
      double *rx, double *ry, double *rz,
      double *Fx, double *Fy, double *Fz
   ){
   //-----------------------------------
   int lmaxe=*lmax+1;
   int LMAXE=(lmaxe)*(lmaxe+2)+1;
   #pragma omp parallel for
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
            double complex GTE[LMAXE],GTM[LMAXE];
            int i=iz+(*nz)*iy+(*nz)*(*ny)*ix;
            rx[i]=x[ix];
            ry[i]=y[iy];
            rz[i]=z[iz];
            vswf_bbp(m,s,p,&lmaxe,NMAX,gamma,kz,&x[ix],&y[iy],&z[iz],GTE,GTM);
            lmie_ofc(M,X,lmax,NMAX,GTE,GTM,&Fx[i],&Fy[i],&Fz[i]);
         } /* end for iz */    
      } /* end for iy */    
   } /* end for ix */    
}
//------------------------------------------------------------------------------
void bscf_gpw(/* FUNCTION */
      int *lmax, int *NMAX,
      double *kx, double *ky, double *kz,
      double complex *ux, double complex *uy, double complex *uz,
      double complex *M, double complex *X,
      double *x, double *y, double *z,
      int *nx, int *ny, int *nz,
      double *rx, double *ry, double *rz,
      double *Fx, double *Fy, double *Fz
   ){
   //-----------------------------------
   int lmaxe=*lmax+1;
   int LMAXE=(lmaxe)*(lmaxe+2)+1;
   #pragma omp parallel for
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
            double complex GTE[LMAXE],GTM[LMAXE];
            int i=iz+(*nz)*iy+(*nz)*(*ny)*ix;
            rx[i]=x[ix];
            ry[i]=y[iy];
            rz[i]=z[iz];
            vswf_gpw(kx,ky,kz,ux,uy,uz,&lmaxe,GTE,GTM);
            lmie_ofc(M,X,lmax,NMAX,GTE,GTM,&Fx[i],&Fy[i],&Fz[i]);
         } /* end for iz */    
      } /* end for iy */    
   } /* end for ix */    
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
