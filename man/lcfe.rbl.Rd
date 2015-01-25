\name{lcfe.rbl}
\alias{lcfe.rbl}
\title{Calculates the logarithmic derivative of Ricatti-Bessel Functions.}
\usage{
lcfe.rbl(n, x, NMAX = 2000, code = "C")
}
\arguments{
  \item{n}{The order of \eqn{C_n(x)}.}

  \item{x}{The argument of \eqn{C_n(x)} of type complex.}

  \item{NMAX}{The maximum number of iterations.}

  \item{code}{Choice between \code{C} or native {R} code.}
}
\value{
The value of \eqn{C_n} for complex arguments.
}
\description{
Calculates the logarithmic derivative of Ricatti-Bessel
Functions.
}
\details{
Using Lentz method is possible to calculate the ratio
\eqn{C_n(x)=(xj_{n}(x))'/xj_{n}(x)}.
}
\seealso{
\code{\link{lcfa.cyl}}, \code{\link{cfe.cbd}}.
}

