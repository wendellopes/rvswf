\name{vwfd.bbz}
\alias{vwfd.bbz}
\title{Wave Field Definition for the Bessel Beam with vector only z component
   for the vector potential.}
\description{Given a position $x,y,z$, it returns the value
of the field at that position for a Bessel Beam Z.}
\usage{
vwfd.bbz(TE=TRUE,M,S,g,kz,x,y,z)
}
\arguments{
   \item{TE}{boolean with 'TRUE' meaning the Transversal
   Eletric Mode and 'FALSE' the Transversal Magnetica Mode.}
   \item{M,S}{M indicates the order of the Bessel function and
   S the sign of the solution.}
   \item{g,kz}{wave number in radial (g^2=kx^2+ky^2) and z direction}
   \item{x,y,z}{positions at which one wants to perform calculations.
   x, y and z can be vectors of diferent sizes.}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbp},\link{vswf.pwe}}
}

