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
/*------------------------------------------------------------------------------
 *                             BESSEL FUNCTIONS                                *
------------------------------------------------------------------------------*/
//------------------------------------------------------------------------------
//# Auxliary Functions
//------------------------------------------------------------------------------
double lcfe_afs(/* FUNCTION */
      int n,
      double x
   ){
   //-----------------------------------
   return(n/x);
}
//------------------------------------------------------------------------------
double complex lcfc_afs(/* FUNCTION */
      int n,
      double complex x
   ){
   //-----------------------------------
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
// Auxiliary function for calculation of Bessel functions by expansion
double KJ( /* FUNCTION */
      int j, 
      double x){
//--------------------------------------
   return(x/j);
}
//------------------------------------------------------------------------------
// Auxiliary function for calculation of Bessel functions by downward recurrency
double SN( /* FUNCTION */
      int j, 
      double x
   ){
//--------------------------------------
   return(j/x);
}
//------------------------------------------------------------------------------
// Cylindrical Bessel function ratio and logarithmic derivative 
// \gamma_n = J_{n}/J_{n+1} == gn;
// D_{n}    = J_{n}'/J_{n}  == Dn;
void lcfa_cyl( /* FUNCTION */
      int *nmax,
      double *x,
      double *gn,
      double *Dn,
      int *NMAX
   ){
//--------------------------------------
   int NMAY=*NMAX;
   // starting values calculated by Lentz
   lcfe_cbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfe_cbd(nmax,x,NMAX,&gn[*nmax]);
   // testing if Lentz reached max number of iterations
   // if yes, so we must increase the number of iteratations
   while(*NMAX==NMAY){
      *NMAX=2*NMAY;
      NMAY=*NMAX;
      lcfe_cbl(nmax,x,NMAX,&Dn[*nmax]);
      lcfe_cbd(nmax,x,NMAX,&gn[*nmax]);
      printf("Precision increased to %d\n",*NMAX);
   }
   // Downward Recursion
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=SN(n  ,*x)+Dn[n];
      Dn[n-1]=SN(n-1,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Spherical Bessel function ratio and logarithmic derivative 
// \rho_{n} = j_{n}/j_{n+1} == gn;
// c_{n}    = j_{n}'/j_{n}  == Dn;
void lcfa_sph( /* FUNCTION */
      int *nmax,
      double *x,
      double *gn,
      double *Dn,
      int *NMAX
   ){
//--------------------------------------
   int NMAY=*NMAX;
   // starting values calculated by Lentz
   lcfe_sbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfe_sbd(nmax,x,NMAX,&gn[*nmax]);
   // testing if Lentz reached max number of iterations
   // if yes, so we must increase the number of iteratations
   while(*NMAX==NMAY){
      *NMAX=2*NMAY;
      NMAY=*NMAX;
      lcfe_sbl(nmax,x,NMAX,&Dn[*nmax]);
      lcfe_sbd(nmax,x,NMAX,&gn[*nmax]);
      printf("Precision increased to %d\n",*NMAX);
   }
   // Downward Recursion
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=SN(n+1,*x)+Dn[n];
      Dn[n-1]=SN(n-1,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Riccati-Bessel function ratio and logarithmic derivative 
// \rho_n = R_{n}/R_{n+1} =j_{n}/j_{n+1} == gn;
// C_{n}    = R_{n}'/R_{n}  == Dn;
void lcfa_ric( /* FUNCTION */
      int *nmax,
      double *x,
      double *gn,
      double *Dn,
      int *NMAX
   ){
//--------------------------------------
   int NMAY=*NMAX;
   // starting values calculated by Lentz
   lcfe_rbl(nmax,x,NMAX,&Dn[*nmax]);
   lcfe_sbd(nmax,x,NMAX,&gn[*nmax]);
   // testing if Lentz reached max number of iterations
   // if yes, so we must increase the number of iteratations
   while(*NMAX==NMAY){
      *NMAX=2*NMAY;
      NMAY=*NMAX;
      lcfe_rbl(nmax,x,NMAX,&Dn[*nmax]);
      lcfe_sbd(nmax,x,NMAX,&gn[*nmax]);
      printf("Precision increased to %d\n",*NMAX);
   }
   // Downward Recursion
   int n;
   for(n=*nmax;n>0;n--){
      gn[n-1]=SN(n,*x)+Dn[n];
      Dn[n-1]=SN(n,*x)-1/gn[n-1];
   }
}
//------------------------------------------------------------------------------
// Spherical Bessel function j_0(x)
void bess_zro( /* FUNCTION */
      double *x, 
      double *j0
   ){
//--------------------------------------
   if(*x>1.){
      *j0=sin(*x)/(*x);
   }else{
      int j;
      for(j=10;j>0;j--){
         *j0=1.-*j0*KJ(*x,2*j)*KJ(*x,2*j+1);
      }
   }
}
//------------------------------------------------------------------------------
// Spherical Bessel function j_1(x)
void bess_uno( /* FUNCTION */
      double *x, 
      double *j1
   ){
//--------------------------------------
   if(*x>1.){
      *j1=(sin(*x)/(*x)-cos(*x))/(*x);
   }else{
      int j;
      for(j=10;j>0;j--){
         *j1=(2*j*1.)/(2*j+1.)-*j1*KJ(*x,2*j+1)*KJ(*x,2*j+2);
      }
      *j1=(*x)*(*j1)/2.0;
   }
}
//------------------------------------------------------------------------------
//Starting values for downward recurrence, Cylindrical Bessel functions
void bess_csv( /* FUNCTION */
      int *nmax,
      int *dig,
      double *x,
      double *JN,
      double *DN
   ){
//--------------------------------------
   double Jn[*dig], Dn[*dig];
   Jn[*dig]=0.0;  // n-th Bessel function
   Dn[*dig]=1.0;  // Derivative of the n-th Bessel function
   int n,j;
   for(n=*dig;n>0;n--){
      // Downward Recursion
      Jn[n-1]=SN(*nmax+n  ,*x)*Jn[n  ]+Dn[n];
      Dn[n-1]=SN(*nmax+n-1,*x)*Jn[n-1]-Jn[n];
      //Renormalization may be required
      if(fabs(Jn[n-1])>1e+2){
         for(j=n-1;j<=*dig;j++){
            Dn[j]=Dn[j]/Jn[n-1];
            Jn[j]=Jn[j]/Jn[n-1];
         }
      }
   }
   *JN=Jn[0];
   *DN=Dn[0];
}
//------------------------------------------------------------------------------
// Array of cylindrical Bessel functions
void bess_cyl( /* FUNCTION */
      int *nmax, 
      double *x, 
      double *Jn, 
      double *Dn, 
      double *J0, 
      double *J1
   ){
//--------------------------------------
   Jn[*nmax]=0.0;  // n-th Bessel function
   Dn[*nmax]=1.0;  // Derivative of the n-th Bessel function
   int dig=15;     // Number of digits of precision
   int n,j;
   bess_csv(nmax,&dig,x,&Jn[*nmax],&Dn[*nmax]);
   for(n=*nmax;n>0;n--){
      // Downward Recursion
      Jn[n-1]=SN(n  ,*x)*Jn[n  ]+Dn[n];
      Dn[n-1]=SN(n-1,*x)*Jn[n-1]-Jn[n];
      //Renormalization may be required
      if(fabs(Jn[n-1])>1e2){
         for(j=n-1;j<=*nmax;j++){
            Dn[j]=Dn[j]/Jn[n-1];
            Jn[j]=Jn[j]/Jn[n-1];
         }
      }
   }
   // Normalization factor for the function
   double pf=1.0;
   if(fabs(*J0)>1e-10){
      pf=*J0/Jn[0];
   }else{
      pf=*J1/Jn[1];
   }
   // Normalization factor for the derivative
   double pd=1.0;
   if(fabs(*J1)>1e-10){
      pd=-(*J1)/Dn[0];
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
//Starting values for downward recurrence, Spherical Bessel functions
void bess_ssv( /* FUNCTION */
      int *nmax,
      int *dig,
      double *x,
      double *JN,
      double *CN
   ){
//--------------------------------------
   //printf("passei\n");
   double jn[*dig], cn[*dig];
   jn[*dig]=1.0;  // n-th Bessel function
   cn[*dig]=0.0;  // Derivative of the n-th Bessel function
   int n,j;
   for(n=*dig;n>0;n--){
      // Downward Recursion
      jn[n-1]=SN(*nmax+n+1,*x)*jn[n  ]+cn[n];
      cn[n-1]=SN(*nmax+n-1,*x)*jn[n-1]-jn[n];
      //Renormalization may be required
      if(fabs(jn[n-1])>1e+2){
         for(j=n-1;j<=*dig;j++){
            cn[j]=cn[j]/jn[n-1];
            jn[j]=jn[j]/jn[n-1];
         }
      }
   }
   *JN=jn[0];
   *CN=cn[0];
}
//------------------------------------------------------------------------------
// Array of Spherical Bessel functions
void bess_sph( /* FUNCTION */
      int *nmax, 
      double *x, 
      double *jn, 
      double *dn
   ){
//--------------------------------------
   double jo=0.0;  // Initialization values
   double ju=0.0;  // Initialization values
   double *j0=&jo;
   double *j1=&ju;
   bess_zro(x,j0);
   bess_uno(x,j1);
   //beginning of the array calculation
   int dig=15;     // Number of digits of precision
   bess_ssv(nmax,&dig,x,&jn[*nmax],&dn[*nmax]);
   int n,j;
   for(n=*nmax;n>0;n--){
      // Downward Recursion
      jn[n-1]=SN(n+1,*x)*jn[n  ]+dn[n];
      dn[n-1]=SN(n-1,*x)*jn[n-1]-jn[n];
      //Renormalization may be required
      if(fabs(jn[n-1])>1e10){
         for(j=n-1;j<=*nmax;j++){
            dn[j]=dn[j]/jn[n-1];
            jn[j]=jn[j]/jn[n-1];
         }
      }
   }
   // Normalization factor for the function
   double pf=1.0;
   if(fabs(*j0)>1e-10){
      pf=*j0/jn[0];
   }else{
      pf=*j1/jn[1];
   }
   // Normalization factor for the derivative
   double pd=1.0;
   if(fabs(*j1)>1e-10){
      pd=-(*j1)/dn[0];
   }else{
      pd=(1./3.)*(jn[0]-2*jn[2])/dn[1];
   }
   // Normalizations
   for(n=0;n<=*nmax;n++){
      jn[n]=jn[n]*pf;
      dn[n]=dn[n]*pd;
   }
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Starting values for downward recurrence, Spherical Bessel functions
void bess_rsv( /* FUNCTION */
      int *nmax,
      int *dig,
      double *x,
      double *RN,
      double *CN
   ){
//--------------------------------------
   double Rn[*dig], Cn[*dig];
   Rn[*dig]=0.0;  // n-th Bessel function
   Cn[*dig]=1.0;  // Derivative of the n-th Bessel function
   int n,j;
   for(n=*dig;n>0;n--){
      // Downward Recursion
      Rn[n-1]=SN(n+*nmax,*x)*Rn[n  ]+Cn[n];
      Cn[n-1]=SN(n+*nmax,*x)*Rn[n-1]-Rn[n];
      //Renormalization may be required
      if(fabs(Rn[n-1])>1e2){
         for(j=n-1;j<=*dig;j++){
            Cn[j]=Cn[j]/Rn[n-1];
            Rn[j]=Rn[j]/Rn[n-1];
         }
      }
   }
   *RN=Rn[0];
   *CN=Cn[0];
}
//------------------------------------------------------------------------------
// Array of Riccati-Bessel functions
void bess_ric( /* FUNCTION */
      int *nmax, 
      double *x, 
      double *Rn, 
      double *Cn
   ){
//--------------------------------------
   double j0=0.0;  // Initialization values
   double j1=0.0;  // Initialization values
   double ro=0.0;
   double ru=0.0;
   double *R0=&ro;
   double *R1=&ru;
   bess_zro(x,&j0);
   bess_uno(x,&j1);
   *R0=*x*j0;
   *R1=*x*j1;
   int dig=15;
   bess_rsv(nmax,&dig,x,&Rn[*nmax],&Cn[*nmax]);
   //printf("Rn = %f, Cn = %f\n",Rn[*nmax],Cn[*nmax]);
   int n,j;
   for(n=*nmax;n>0;n--){
      // Downward Recursion
      Rn[n-1]=SN(n,*x)*Rn[n  ]+Cn[n];
      Cn[n-1]=SN(n,*x)*Rn[n-1]-Rn[n];
      //Renormalization may be required
      if(fabs(Rn[n-1])>1e2){
         for(j=n-1;j<=*nmax;j++){
            Cn[j]=Cn[j]/Rn[n-1];
            Rn[j]=Rn[j]/Rn[n-1];
         }
      }
   }
   // Normalization factor for the function
   double pf=1.0;
   if(fabs(*R0)>1e-10){
      pf=(*R0)/Rn[0];
   }else{
      pf=(*R1)/Rn[1];
   }
   // Normalization factor for the derivative
   double pd=1.0;
   if(fabs(*R1)>1e-10){
      pd=cos(*x)/Cn[0];
   }else{
      pd=(-j1+*R0)/Cn[1];
   }
   // Normalizations
   for(n=0;n<=*nmax;n++){
      Rn[n]=Rn[n]*pf;
      Cn[n]=Cn[n]*pd;
   }
}
/*------------------------------------------------------------------------------
 *                  VECTOR SPHERICAL WAVE FUNCTIONS                            *
------------------------------------------------------------------------------*/
//------------------------------------------------------------------------------
// Psi_m(\vec{k},\vec{r}): Basic function for cylindrical simetries
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
   double jn=gsl_sf_bessel_Jn(*m,rho*(*gamma));
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
// VECTOR SPHERICAL WAVE FUNCTIONS CALCULATIONS
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
double cp(/* FUNCTION */
      int lo, 
      int mo
   ){
   //-----------------------------------
   double l=lo;
   double m=mo;
   return sqrt(l*(l+1)-m*(m+1));
}
//------------------------------------------------------------------------------
double cm(/* FUNCTION */
      int lo, 
      int mo
   ){
   //-----------------------------------
   double l=lo;
   double m=mo;
   return sqrt(l*(l+1)-m*(m-1));
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
   /* SPHERICAL BESSEL FUNCTIONS      */
   /*---------------------------------*/
   double JLM[*lmax+2];
   gsl_sf_bessel_jl_steed_array(*lmax+1,*k*r,JLM);
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
/*----------------------------------------------------------------------------*/
