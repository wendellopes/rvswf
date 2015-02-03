\name{lcfe.cbl}
\alias{lcfe.cbl}
\title{Calculates the logarithmic derivative of Cylindrical Bessel Functions.}
\usage{
lcfe.cbl(n, x, NMAX = 2000, code = "C")
}
\arguments{
  \item{n}{The order of \eqn{D_n(x)}.}

  \item{x}{The argument of \eqn{D_n(x)} of type complex.}

  \item{NMAX}{The maximum number of iterations.}

  \item{code}{Choice between \code{C} or native {R} code.}
}
\value{
The value of \eqn{D_n} for complex arguments.
}
\description{
Calculates the logarithmic derivative of Cylindrical Bessel
Functions.
}
\details{
Using Lentz method is possible to calculate the ratio
\eqn{D_n(x)=J_{n}'(x)/J_{n}(x)}. By downward recurrence one
can calculate \eqn{D_n}.
}
\seealso{
\code{\link{lcfa.cyl}}, \code{\link{lcfe.cbd}}.
}

