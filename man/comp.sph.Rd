\name{comp.sph}
\alias{comp.sph}
\title{Compare results for Spherical Bessel Functions.}
\usage{
comp.sph(nmax, x)
}
\arguments{
  \item{x}{The argument of \eqn{j_n(x)}.}

  \item{n}{The order of the Spherical Bessel function.}
}
\description{
Compare results for Spherical Bessel Functions.
}
\details{
Compare results using vswf, built in \code{R} and
\code{gsl} algorithms.
}
\examples{
x<-5
nmax<-10
print(comp.sph(5,3))
}

