\name{vswf.hmp}
\alias{vswf.hmp}
\title{Hansen Multipoles}
\description{Calculates Hansen Multipoles (Partial Waves M_{lm} and N_{lm})}
\usage{
vswf.hmp(k,x,y,z,lmax)
}
\arguments{
  \item{k}{Wavenumber}
  \item{x,y,z}{Position}
  \item{lmax}{Maximum number of l terms}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
lmax=5
u<-vswf.hmp(1,1,1,1,lmax)
print(u)
}
