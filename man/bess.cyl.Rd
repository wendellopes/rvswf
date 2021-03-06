% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/bess.cyl.r
\name{bess.cyl}
\alias{bess.cyl}
\title{Calculates Cylindrical Bessel functions from 0 to nmax.}
\usage{
bess.cyl(nmax, x, code = "C")
}
\arguments{
\item{nmax}{The maximum order of \eqn{J_n(x)}}

\item{x}{The argument of \eqn{J_n(x)}}

\item{code}{If you prefer to use native R or C language.
The algorithm is the same.}
}
\value{
An array of Cylindrical Bessel functions and its derivatives.
from 0 to \code{nmax} at point \code{x}
}
\description{
Calculates Cylindrical Bessel functions from 0 to nmax.
}
\details{
\code{bess.cyl} calculates the Cylindrical Bessel
function using downward recurrence, from \eqn{J_nmax(x)} to \eqn{J_0(x)}.
The system of equations is given by \eqn{S_n(x)=n/x},
\eqn{\gamma_n=J_n(x)/J_{n+1}(x)}{g[n]=J_n/J_{n+1}} and
\eqn{D_n=J_n'(x)/J_n(x)}. The system can be solved by means of
the recurrence relations of the Cylindrical Bessel functions
\deqn{\gamma_{n-1}+\frac{1}{\gamma_n}=2S_{n}}{g[n-1]+1/g[n]=2 S[n]}
\deqn{\gamma_{n-1}-\frac{1}{\gamma_n}=2D_{n}}{g[n-1]-1/g[n]=2 D[n]}
that can be rewriten
\deqn{\gamma_{n}=S_{n+1}+D_{n+1}            }{  g[n]=S[n+1]+D[n+1]}
\deqn{\frac{1}{\gamma_n}=S_n-D_n.           }{1/g[n]=S[n  ]-D[n  ].}
The logarithmic derivatives obeys the relation,
\deqn{(S_n-D_n)(S_{n+1}+D_{n+1})=1.         }{(S[n]-D[n])(S[n+1]+D[n+1])=1.}
The values can be calculated upward or downward.
}
\examples{
x<-30
nmax<-50
a<-bess.cyl(nmax,x,code="C")
b<-bess.cyl(nmax,x,code="R")
d<-besselJ(x=x,nu=0:nmax)
plot(a$Jn,type='b')
points(b$Jn,col='red',pch=4)
points(d,col='blue',pch=3)
}

