\name{comp.ric}
\alias{comp.ric}
\title{Compare results for Ricatti-Bessel Functions.}
\usage{
comp.ric(nmax, x)
}
\arguments{
  \item{x}{The argument of \eqn{\psi_n(x)=xj_n(x)}.}

  \item{n}{The order of the Ricatti-Bessel function.}
}
\value{
Table comparing built-in \code{R} functions, \code{gsl} and
\code{rvswf}.
}
\description{
Compare results for Ricatti-Bessel Functions.
}
\details{
Compare results using vswf, built in \code{R} and
\code{gsl} algorithms.
}
\examples{
x<-5
nmax<-10
print(comp.ric(5,3))
}

