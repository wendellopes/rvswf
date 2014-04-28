\name{vswf.pwe}
\alias{vswf.pwe}
\title{Partial Wave Expansion}
\description{Partial Wave Expansion in Vector Spherical Wave Functions by
   means of the Beam Shape Coefficients.}
\usage{
vswf.pwe(k,x,y,z,lmax,gte,gtm)
}
\arguments{
  \item{k}{the modulus of the wavenumber}
  \item{x,y,z}{the positions. x, y and z can be
   vectors of independent sizes}
  \item{lmax}{maximum l of the summation}
  \item{gte,gtm}{Beam Shape Coefficients}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{

#Plane Wave

#Rectangular Wave Guide

#Cylindrical Wave Guide

#Bessel Beam Z

#Bessel Beam P
}
