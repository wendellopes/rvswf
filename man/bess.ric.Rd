\name{bess.ric}
\alias{bess.ric}
\title{Calculates Ricatti-Bessel functions from 0 to nmax.}
\usage{
bess.ric(nmax, x, code = "C")
}
\arguments{
  \item{nmax}{The maximum order of \eqn{R_n(x)}}

  \item{x}{The argument of \eqn{R_n(x)}}

  \item{code}{If you prefer to use native R or C language.
  The algorithm is the same.}
}
\value{
An array of Ricatti-Bessel functions and its derivatives.
from 0 to \code{nmax} at point \code{x}
}
\description{
Calculates Ricatti-Bessel functions from 0 to nmax.
}
\details{
\code{bess.ric} calculates the Ricatti-Bessel functions
using downward recurrence, from \eqn{R_nmax(x)} to
\eqn{R_0(x)}. The system of equations is given by
\eqn{S_n(x)=n/x},
\eqn{\rho_n=R_n(x)/R_{n+1}(x)}{r[n]=R_n/R_{n+1}} and
\eqn{C_n=R_n'(x)/R_n(x)}. The system can be solved by means
of the recurrence relations of the Ricatti-Bessel functions
\deqn{ \rho_{n-1}+\frac{1}{\rho_n}=S_{2n+1}
}{r[n-1]+1/r[n]=S[2n+1]}
\deqn{(n+1\rho_{n-1}-\frac{n}{\rho_n}=(2n+1)C_{n}}{r[n-1]-n/r[n]=(2n+1)C[n]}
that can be rewriten \deqn{\rho_{n}=S_{n+1}+C_{n+1} }{
r[n]=S[n+1]+C[n+1]} \deqn{\frac{1}{\rho_n}=S_{n+1}-C_n.
}{1/r[n]=S[n+1]-C[n ].} The logarithmic derivatives obeys
the relation, \deqn{(S_{n+1}-C_n)(S_{n+1}+C_{n+1})=1.
}{(S[n+1]-C[n])(S[n+1]+C[n+1])=1.} The values can be
calculated upward or downward.
}
\examples{
x<-30
nmax<-50
a<-bess.ric(nmax,x,code="C")
b<-bess.ric(nmax,x,code="R")
d<-sqrt(x*pi/2)*besselJ(x=x,nu=.5+(0:nmax))
plot(a$Rn,type='b')
points(b$Rn,col='red',pch=4)
points(d,col='blue',pch=3)
}

