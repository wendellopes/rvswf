\name{comp.cyl}
\alias{comp.cyl}
\title{Compare results for Cylindrical Bessel Functions.}
\usage{
comp.cyl(nmax, x)
}
\arguments{
  \item{x}{The argument of \eqn{J_n(x)}.}

  \item{n}{The order of the Cylindrical Bessel function.}
}
\description{
Compare results for Cylindrical Bessel Functions.
}
\details{
Compare results using vswf, built in \code{R} and
\code{gsl} algorithms.
}
\examples{
x<-5
nmax<-10
print(comp.cyl(5,3))
}

