\name{lcfe.sbd}
\alias{lcfe.sbd}
\title{Calculates the inverse of the ratio of Spherical Bessel Functions.}
\usage{
lcfe.sbd(n, x, NMAX = 2000, code = "C")
}
\arguments{
  \item{n}{The order of \eqn{\rho_n(x)}.}

  \item{x}{The argument of \eqn{\rho_n(x)} of type
  complex.}

  \item{NMAX}{The maximum number of iterations.}

  \item{code}{Choice between \code{C} or native {R} code.}
}
\value{
The value of \eqn{\rho_n} for complex arguments.
}
\description{
Calculates the inverse of the ratio of Spherical Bessel
Functions.
}
\details{
Using Lentz method is possible to calculate the ratio
\eqn{\rho_n(x)=j_{n}(x)/j_{n+1}(x)}. By downward recurrence
one can calculate \eqn{\rho_n}, as shown in Barnett's
paper.
}
\seealso{
\code{\link{lcfa.cyl}}, \code{\link{lcfe.cbd}}.
}

