/*------------------------------------------------------------------------------
#include <gsl/gsl_sf_bessel.h>
#include <omp.h>
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
