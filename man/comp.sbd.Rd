\name{comp.sbd}
\alias{comp.sbd}
\title{Checks the results of Lentz method and internal R Bessel function}
\usage{
comp.sbd(n, x, code = "C")
}
\arguments{
  \item{n}{Order of the ratio given by
  \eqn{\gamma_n=j_n/j_{n+1}}{r[n]=j_n/j_{n+1}}.}

  \item{x}{Argument of Bessel functions.}

  \item{code}{If C or native R function.}
}
\value{
Data frame with the values calculated by the algorithms.
}
\description{
Checks the results of Lentz method and internal R Bessel
function
}
\details{
In order to check the results of Lentz method for
calculation of Cylindrical Bessel ratio of consecutive
function \eqn{\rho_n}{r[n]}.
}
\examples{
comp.sbd(5,4,code="C")
comp.sbd(5,4,code="R")
}
\seealso{
\code{\link{lcfe.sbi}}, \code{\link{lcfe.sbl}},
\code{\link{lcfe.sbd}}, \code{\link{lcfe.afs}},
\code{\link{besselJ}}.
}

