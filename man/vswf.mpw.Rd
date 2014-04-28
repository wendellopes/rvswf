\name{vswf.mpw}
\alias{vswf.mpw}
\title{Mie Plane Wave}
\description{Beam Shape Coefficients for Plane Wave propagating along z axis,
   using the results of Jackon's Classical Electrodynamics book.}
\usage{
vswf.mpw(lmax,norm=FALSE,s=1)
}
\arguments{
   \item{lmax}{Maximum number of l terms}
   \item{norm}{If 'FALSE' we get the same result of Jackson's book; if
      'TRUE' we obtain the result divided by sqrt(2).}
   \item{s}{Chirality of the plane wave}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
lmax<-4
c3<-vswf.mpw(lmax)
\dontrun{
ylim<-range(abs(c3))
plot(abs(c3$GTE.in),type='l',ylim=ylim)
points(abs(c3$GTM.in),type='l',col='red')
}
}
