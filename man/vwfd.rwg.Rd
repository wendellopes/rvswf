\name{vwfd.rwg}
\alias{vwfd.rwg}
\title{Wave Field Definition for the Rectangular Wave Guide}
\description{Given a position x,y,z, it returns the value
of the field at that position for a Rectangular Wave Guide.}
\usage{
vwfd.rwg(TE=TRUE,kx,ky,kz,x,y,z)
}
\arguments{
  \item{TE}{boolean with 'TRUE' meaning the Transversal
   Eletric Mode and 'FALSE' the Transversal Magnetica Mode.}
  \item{kx,ky,kz}{wavenumber in the x,y,z directions.}
  \item{x,y,z}{positions at which one wants to perform calculations.
   x, y and z can be vectors of diferent sizes.}
}
\seealso{
   \code{\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp},\link{vswf.pwe}}
}

