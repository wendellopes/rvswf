\name{vswf.gpw}
\alias{vswf.gpw}
\alias{bsc}
\alias{vswf}
\title{Generic Plane Wave}
\description{Beam Shape Coefficients for the Generic Plane Wave}
\usage{
vswf.gpw(kx,ky,kz,ux,uy,uz,lmax,xo=0,yo=0,zo=0)
}
\arguments{
   \item{kx,ky,kz}{Wavenumbers components in the x, y and z directions 
   (direction of propagation)}
   \item{ux,uy,uz}{Components of the polarization of the electric field}
   \item{lmax}{Maximum number of partial waves}
   \item{xo,yo,zo}{Origin at which one wants to calculate}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\note{The result is a data.frame with columns GTE/GTM calculated by
   implicit method (general expression using Vector Spherical Harmonics) or
   explicit method (general expression using Scalar Spherical Harmonics).
}
\examples{
# Plane Wave propagating along x direction, with polarization 
# of electric field along x direction: 
# k=kx, E=Ex
c1<-vswf.gpw(1,0,0,1,0,0,4)
# k=ky, E=Ez
c2<-vswf.gpw(0,1,0,0,0,1,4)
# k=kz, E=Ex
c3<-vswf.gpw(0,1,0,1,0,0,4)
\dontrun{
ylim<-range(abs(c3))
plot(abs(c3$GTE.in),type='l',ylim=ylim)
points(abs(c3$GTM.in),type='l',col='red')
}
}
