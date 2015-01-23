\name{vswf.sbf}
\alias{vswf.sbf}
\title{Calculation of Spherical Bessel Functions.}
\usage{
vswf.sbf(lmax, xo, compare = FALSE, verbose = FALSE)
}
\arguments{
  \item{lmax}{The maximum value of \eqn{l}.}

  \item{xo}{The argument of \eqn{j_l}. Can be real or
  complex.}

  \item{compare}{Make the comparison with the GSL (Gnu
  Scientific Library).}

  \item{verbose}{Shows the normalizations on downward
  recurrence.}
}
\value{
An array with the Spherical Bessel functions and its
derivatives.
}
\description{
Calculation of Spherical Bessel Functions.
}
\details{
We use Lentz method for calculating the starting value
j_lmax(x)=bo+(a1/b1+)(a2/b2+)(a3/b3+)...(an/bn), and then
we use downward recurrence to calculate the other values,
from lmax-1 to 0.
}
\examples{
vswf.sbf(5,3,compare=TRUE)
}
\seealso{
\code{\link{bess.sph}}, \code{\link{bess.cyl}},
\code{\link{bess.ric}}, \code{\link{bess.czr}},
\code{\link{bess.cun}}, \code{\link{bess.szr}},
\code{\link{bess.sun}}.
}

