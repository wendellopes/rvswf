\name{vwfd.bbp}
\alias{vwfd.bbp}
\alias{bbp}
\alias{Bessel}
\alias{Beam}
\title{Vector Wave Field Definition for the Bessel Beam with only vector p 
   component for the vector potential, where p=x+iy.}
\description{Given a position $x,y,z$, it returns the value
of the field at that position for a Bessel Beam P.}
\usage{
vwfd.bbp(P,M,S,g,kz,x,y,z)
}
\arguments{
   \item{P}{Plus or minus one,the direction of the vector potential.}
   \item{M,S}{M indicates the order of the Bessel function and
   S the sign of the solution.}
   \item{g,kz}{wave number in radial (g^2=kx^2+ky^2) and z direction}
   \item{x,y,z}{positions at which one wants to perform calculations.
   x, y and z can be vectors of diferent sizes.}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vswf.pwe}}
}
