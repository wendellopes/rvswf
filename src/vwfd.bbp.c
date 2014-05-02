/*------------------------------------------------------------------------------
 * VSWF: Vector Spherical Wave Functions dynamic libraries written in C.
 * author: Wendel Lopes Moreira <wendellopes@gmail.com>
 * date: 2013-02-18
 * version 1.1
 * depends: Gnu Scientific Library <http://www.gnu.org/software/gsl/>
------------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
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
   int ix, iy, iz;
   for(ix=0; ix<*nx; ix++){
      for(iy=0; iy<*ny; iy++){
         for(iz=0; iz<*nz; iz++){
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
