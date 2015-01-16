\name{bess.sph}
\alias{bess.sph}
\title{Calculates Spherical Bessel functions from 0 to nmax.}
\usage{
bess.sph(nmax, x, code = "C")
}
\arguments{
  \item{nmax}{The maximum order of \eqn{j_n(x)}}

  \item{x}{The argument of \eqn{j_n(x)}}

  \item{code}{If you prefer to use native R or C language.
  The algorithm is the same.}
}
\value{
An array of Spherical Bessel functions and its derivatives
from 0 to \code{nmax} at point \code{x}
}
\description{
Calculates Spherical Bessel functions from 0 to nmax.
}
\details{
\code{bess.sph} calculates the Spherical Bessel functions
using downward recurrence, from \eqn{j_nmax(x)} to
\eqn{j_0(x)}. The system of equations is given by
\eqn{S_n(x)=n/x},
\eqn{\rho_n=j_n(x)/j_{n+1}(x)}{r[n]=j_n/j_{n+1}} and
\eqn{c_n=j_n'(x)/j_n(x)}. The system can be solved by means
of the recurrence relations of the Spherical Bessel
functions \deqn{ \rho_{n-1}+\frac{1 }{\rho_n}=S_{2n+1} }{
r[n-1]+ 1/r[n]=S[2n+1]}
\deqn{n\rho_{n-1}-\frac{n+1}{\rho_n}=(2n+1)c_{n}}{nr[n-1]-(n+1)/r[n]=(2n+1)c[n]}
that can be rewriten \deqn{\rho_n=S_{n+2}+c_{n+1} }{
r[n]=S[n+2]+c[n+1]} \deqn{\frac{1}{\rho_n}=S_n-c_n.
}{1/r[n]=S[n ]-c[n ].} The logarithmic derivatives obeys
the relation, \deqn{(S_{n+2}-c_{n+1})(S_n+c_n)=1.
}{(S[n+2]-c[n])(S[n]+C[n])=1.} The values can be calculated
upward or downward.
}
\examples{
x<-5
nmax<-50
a<-bess.sph(nmax,x,code="C")
b<-bess.sph(nmax,x,code="R")
d<-sqrt(pi/(2*x))*besselJ(x=x,nu=.5+(0:nmax))
plot(a$jn)
points(b$jn,col='red',pch=4)
points(d,col='blue',pch=3)
}

