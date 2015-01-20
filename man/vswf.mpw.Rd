\name{vswf.mpw}
\alias{vswf.mpw}
\title{Beam Shape Coefficients for a Plane Wave, as done by Mie}
\usage{
vswf.mpw(lmax, norm = FALSE, s = 1)
}
\arguments{
  \item{lmax}{The maximum value of \eqn{l}.}

  \item{norm}{If TRUE, the Beam Shape Coefficient will be
  divided by \eqn{\sqrt{2}}.}

  \item{s}{The polarity of the plane wave.}
}
\description{
Beam Shape Coefficients for a Plane Wave, as done by Mie
}
\details{
We use the expression showed on Jackson's Classical
Electrodynamics book.
}
\examples{
lm<-5
a<-vswf.mpw(lm)
plot(Re(a$GTE),type='b')
points(Im(a$GTE),pch=4,col='red',type='b')
plot(Re(a$GTM),type='b')
points(Im(a$GTM),pch=4,col='red',type='b')
}
\seealso{
\code{\link{vswf.gpw}}, \code{\link{vswf.jlm}}.
}

