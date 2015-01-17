\name{comp.cbl}
\alias{comp.cbl}
\title{Checks the results of Lentz method and internal R Bessel function}
\usage{
comp.cbl(n, x, code = "C")
}
\arguments{
  \item{n}{Order of the logarithmic derivative given by
  \eqn{D_n=J_n'/J_n}.}

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
\eqn{D_n}.
}
\examples{
comp.cbl(5,4,code="C")
comp.cbl(5,4,code="R")
}
\seealso{
\code{\link{lcfe.cbi}}, \code{\link{lcfe.cbl}},
\code{\link{lcfe.cbd}}, \code{\link{lcfe.afs}},
\code{\link{besselJ}}.
}

