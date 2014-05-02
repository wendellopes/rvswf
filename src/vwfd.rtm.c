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
#include <gsl/gsl_sf_bessel.h>
#include <omp.h>
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
