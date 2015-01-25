\name{lcfe.sbl}
\alias{lcfe.sbl}
\title{Calculates the logarithmic derivative of Spherical Bessel Functions.}
\usage{
lcfe.sbl(n, x, NMAX = 2000, code = "C")
}
\arguments{
  \item{n}{The order of \eqn{c_n(x)}.}

  \item{x}{The argument of \eqn{c_n(x)} of type complex.}

  \item{NMAX}{The maximum number of iterations.}

  \item{code}{Choice between \code{C} or native {R} code.}
}
\value{
The value of \eqn{c_n} for complex arguments.
}
\description{
Calculates the logarithmic derivative of Spherical Bessel
Functions.
}
\details{
Using Lentz method is possible to calculate the ratio
\eqn{c_n(x)=j_{n}'(x)/j_{n}(x)}.
}
\seealso{
\code{\link{lcfa.cyl}}, \code{\link{cfe.cbd}}.
}

