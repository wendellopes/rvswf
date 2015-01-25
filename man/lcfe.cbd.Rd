\name{lcfe.cbd}
\alias{lcfe.cbd}
\title{Calculates the ratio of Cylindrical Bessel Functions.}
\usage{
lcfe.cbd(n, x, NMAX = 2000, code = "C")
}
\arguments{
  \item{n}{The order of \eqn{\gamma_n(x)}.}

  \item{x}{The argument of \eqn{\gamma_n(x)} of type
  complex.}

  \item{NMAX}{The maximum number of iterations.}

  \item{code}{Choice between \code{C} or native {R} code.}
}
\value{
The value of \eqn{\gamma_n} for complex arguments.
}
\description{
Calculates the ratio of Cylindrical Bessel Functions.
}
\details{
Using Lentz method is possible to calculate the ratio
\eqn{\gamma_n(x)=J_{n}(x)/J_{n+1}(x)}. By downward
recurrence one can calculate \eqn{\gamma_n}.
}
\seealso{
\code{\link{lcfa.cyl}}, \code{\link{cfe.cbd}}.
}

