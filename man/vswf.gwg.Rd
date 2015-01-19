\name{vswf.gwg}
\alias{vswf.gwg}
\title{General expression for wave guides.}
\usage{
vswf.gwg(gama, kz, lmax)
}
\arguments{
  \item{gamma}{The transversal component of the wave vector
  \eqn{\gamma}.}

  \item{kz}{The longitudinal component of the wave vector
  \eqn{k_z}.}

  \item{lmax}{The maximum value of \eqn{l} of the
  expansion.}
}
\value{
The values \eqn{A} and \eqn{B}.
}
\description{
General expression for wave guides.
}
\details{
This expression does not take in account the value of the
Fourier transform of azimuthal component. It does not
depend on the position of the expansion, only on the wave
vector.
}
\examples{
th<-pi/3
u<-vswf.gwg(sin(th),cos(th),10)
plot(Re(u$A),type='b')
points(Im(u$A),type='b',pch=4,col='red')
plot(Re(u$B),type='b')
points(Im(u$B),type='b',pch=4,col='red')
}
\seealso{
\code{\link{vswf.cwg}}, \code{\link{vswf.rwg}}.
}

