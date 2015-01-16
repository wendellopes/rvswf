\name{comp.cbd}
\alias{comp.cbd}
\title{Checks the results of Lentz method and internal R Bessel function}
\usage{
comp.cbd(n, x)
}
\arguments{
  \item{n}{Order of the ratio given by
  \eqn{\gamma_n=J_n/J_{n+1}}{g[n]=J_n/J_{n+1}}.}

  \item{x}{Argument of Bessel functions.}
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
function \eqn{\gamma_n}{g[n]}.
}
\examples{
comp.cbd(5,4)
}
\seealso{
\code{\link{lcfe.cbi}}, \code{\link{lcfe.cbl}},
\code{\link{lcfe.cbd}}, \code{\link{lcfe.afs}},
\code{\link{besselJ}}.
}

