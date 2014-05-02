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
