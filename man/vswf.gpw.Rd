\name{vswf.gpw}
\alias{vswf.gpw}
\title{Beam Shape Coefficients for a Generic Plane Wave}
\usage{
vswf.gpw(kx, ky, kz, ux, uy, uz, lmax, xo = 0, yo = 0, zo = 0,
  implicit = TRUE)
}
\arguments{
  \item{kx}{The component \eqn{x} of the wave vector.}

  \item{ky}{The component \eqn{y} of the wave vector.}

  \item{kz}{The component \eqn{z} of the wave vector.}

  \item{x}{The component \eqn{x} of the position vector.}

  \item{y}{The component \eqn{y} of the position vector.}

  \item{z}{The component \eqn{z} of the position vector.}

  \item{xo}{The component \eqn{x} of the origin vector.}

  \item{yo}{The component \eqn{y} of the origin vector.}

  \item{zo}{The component \eqn{z} of the origin vector.}

  \item{implicit}{The way it calculates the Beam Shape
  Coefficients. If TRUE, so we use the function
  \code{vswf.vsh} to calculate. If FALSE, we use the
  explicit expression for the components.}
}
\description{
Beam Shape Coefficients for a Generic Plane Wave
}
\details{
Calculate the Beam Shape Coefficients for a plane wave with
arbitrary polarization and propagating in arbitrary
direction. It can be compared with the specific case
calculated by Mie or that done in the book of Jackson,
called here by Mie Plane Wave.
}
\examples{
lm<-5
a<-vswf.mpw(lm,norm=TRUE)
b<-vswf.gpw(0,0,1,1/sqrt(2),1i/sqrt(2),0,lm)
#
plot(Re(a$GTE),type='b')
points(Re(b$GTE),pch=4,col='red',type='b')
#
plot(Re(a$GTM),type='b')
points(Re(b$GTM),pch=4,col='red',type='b')
#
plot(Im(a$GTE),type='b')
points(Im(b$GTE),pch=4,col='red',type='b')
#
plot(Im(a$GTM),type='b')
points(Im(b$GTM),pch=4,col='red',type='b')
}
\seealso{
\code{\link{vswf.vsh}}, \code{\link{vswf.mpw}},
\code{\link{vswf.qlm}}, \code{\link{vswf.jlm}}.
}

