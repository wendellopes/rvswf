\name{vswf.vsh}
\alias{vswf.vsh}
\title{Scalar and Vector Spherical Harmonics}
\description{Calculate Spherical Harmonics}
\usage{
vswf.vsh(x,y,z,lmax)
}
\arguments{
  \item{x,y,z}{Position in cartesian coordinates}
  \item{lmax}{Maximum number of l terms}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
th<-pi/3
ph<-pi/6
lmax<-4
u<-vswf.vsh(sin(th)*cos(ph),sin(th)*sin(ph),cos(th),lmax)
#Comparisons of results
}

