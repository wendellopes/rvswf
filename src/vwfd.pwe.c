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
   int ix, iy, iz;
   #pragma omp parallel for
   for(ix=0; ix<*nx; ix++){
      for(iy=0; iy<*ny; iy++){
         for(iz=0; iz<*nz; iz++){
           //i=*ny*ix+*nz*iy+iz;
           i=iz+*nz*(iy+*ny*ix);
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
