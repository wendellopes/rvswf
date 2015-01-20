\name{vswf.hmp}
\alias{vswf.hmp}
\title{Hansen Multipoles.}
\usage{
vswf.hmp(k, x, y, z, lmax)
}
\arguments{
  \item{k}{The module of the wave vector.}

  \item{x}{The component \eqn{x} of the position vector.}

  \item{y}{The component \eqn{y} of the position vector.}

  \item{z}{The component \eqn{z} of the position vector.}

  \item{lmax}{The maximum value of \eqn{l}.}
}
\value{
A list with the complex components of the Hansen Multipoles
\eqn{M_{lm}} and \eqn{N_{lm}}.
}
\description{
Hansen Multipoles.
}
\details{
The Hansen Multipoles combine the Vector Spherical
Harmonics with the spherical Bessel functions.
}
\examples{
th<-pi/3
ph<-pi/4
VSH<-vswf.hmp(1,sin(th)*cos(ph),sin(th)*sin(ph),cos(th),5)
print(as.data.frame(VSH))
}
\seealso{
\code{\link{vswf.sbf}}, \code{\link{vswf.vsh}}.
}

