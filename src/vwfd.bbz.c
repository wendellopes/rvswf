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
   int ix, iy, iz;
   for(ix=0; ix<*nx; ix++){
      for(iy=0; iy<*ny; iy++){
         for(iz=0; iz<*nz; iz++){
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
