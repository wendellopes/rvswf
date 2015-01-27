\name{lcfa.sph}
\alias{lcfa.sph}
\title{Array of logarithmic derivative and ratio of Spherical Bessel functions.}
\usage{
lcfa.sph(nmax, x, code = "C", NMAX = 200)
}
\arguments{
  \item{nmax}{Maximun value of \code{n}.}

  \item{x}{The argument of the functions. Can be complex.}

  \item{code}{If the \code{C} code or \code{R} built in.}

  \item{NMAX}{Maximum number of iterations on the Lentz
  algorithm.}
}
\description{
Array of logarithmic derivative and ratio of Spherical
Bessel functions.
}
\details{
Calculate by downward recurrence the derivative of
Spherical Bessel functions \eqn{c_n=j_n'(x)/j_n(x)} and
simultaneously the ratio \eqn{\rho_n=j_n(x)/j_{n+1}(x)},
from 0 to \code{nmax}. The starting values are calculated
by Lentz continued fraction method.
}
\examples{
nmax<-10
x<-5
u.c<-lcfa.sph(nmax,x,code="C")
u.r<-lcfa.sph(nmax,x,code="R")
u<-data.frame(
   # Logarithmic Derivatives
   C.LogDev=u.c$cn,R.LogDev=u.r$cn,
   # Ratio between Spherical Bessel functions
   C.SphRat=u.c$rn,R.SphRat=u.r$rn)
print(u)
}
\seealso{
\code{\link{lcfe.sbl}}, \code{\link{lcfe.sbd}},
\code{\link{lcfe.sbi}}.
}

