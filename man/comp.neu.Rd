\name{comp.neu}
\alias{comp.neu}
\title{Compare results for Cylindrical Neumann Functions.}
\usage{
comp.neu(n, x)
}
\arguments{
  \item{x}{The argument of \eqn{\psi_n(x)=xj_n(x)}.}

  \item{n}{The order of the Ricatti-Bessel function.}
}
\description{
Compare results for Cylindrical Neumann Functions.
}
\details{
Compare results using vswf, built in \code{R} and
\code{gsl} algorithms.
}
\examples{
x<-5
nmax<-10
print(comp.neu(5,3))
}

