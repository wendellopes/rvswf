\name{vswf.gwg}
\alias{vswf.gwg}
\title{Beam Shape Coefficients for Generic Wave Guides}
\description{This generic function calculate the Beam Shape Coefficients for
   the Wave Guides. It depends of the specific functions for each Wave Guide 
   type, which calculate the Fourier Series coefficients of the tangential
   component of the fields.}
\usage{
vswf.gwg(gama, kz, lmax)
}

\arguments{
   \item{gama,kz}{Components of the wave vector}
   \item{lmax}{Maximum order of l terms to sum}
}
\seealso{
   \code{
      \link{vswf.rwg},
      \link{vswf.cwg},
      \link{vwfd.rwg},
      \link{vwfd.cwg},
      \link{vwfd.bbz},
      \link{vwfd.bbp}
      }
}
\examples{
lmax=5
u<-vswf.gwg(sin(3),cos(3),lmax)
print(u)
}
