\name{comp.sbl}
\alias{comp.sbl}
\title{Checks the results of Lentz method and internal R Bessel function}
\usage{
comp.sbl(n, x, code = "C")
}
\arguments{
  \item{n}{Order of the logarithmic derivative given by
  \eqn{c_n=j_n'/j_n}.}

  \item{x}{Argument of Bessel functions.}

  \item{code}{If C or native R function.}
}
\value{
Data frame with the values calculated by the algorithm.
}
\description{
Checks the results of Lentz method and internal R Bessel
function
}
\details{
In order to check the results of Lentz method for
calculation of Cylindrical Bessel logarithmic derivative
\eqn{c_n}.
}
\examples{
comp.sbl(5,4,code="C")
comp.sbl(5,4,code="R")
}
\seealso{
\code{\link{lcfe.sbi}}, \code{\link{lcfe.sbl}},
\code{\link{lcfe.sbd}}, \code{\link{lcfe.afs}},
\code{\link{besselJ}}.
}

