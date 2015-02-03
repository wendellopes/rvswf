\name{lmie.log}
\alias{lmie.log}
\title{Calculates the Lorentz-Mie coefficients \eqn{a_n} and \eqn{b_n}.}
\usage{
lmie.log(m, x, NMAX = floor(abs(x + 7.5 * x^(1/3)) + 2))
}
\arguments{
  \item{m}{The ratio between the refractive indices.}

  \item{x}{The form factor value.}

  \item{NMAX}{The maximum value of \code{n}. Default value
  from paper.}
}
\value{
The coefficients \eqn{a_n=C_nT_a} and \eqn{b_n=C_nT_b} and
also the coefficients \eqn{C_n}, \eqn{T_a} and \eqn{T_b}.
The maximum value of \code{n}.
}
\description{
Calculates the Lorentz-Mie coefficients \eqn{a_n} and
\eqn{b_n}.
}
\details{
Calculation done using the logarithmic derivative of
Ricatti-Bessel functions.
}
\examples{
m<-1.2+.1i
x<-3
print(lmie.rho(1.2,5))
}

