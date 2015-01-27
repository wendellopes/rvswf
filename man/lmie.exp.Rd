\name{lmie.exp}
\alias{lmie.exp}
\title{Calculates the Lorentz-Mie coefficients \eqn{a_n} and \eqn{b_n}.}
\usage{
lmie.exp(m, x, by = "LD", ...)
}
\arguments{
  \item{m}{The ratio between the refractive indices.}

  \item{x}{The form factor value.}

  \item{by}{\code{LD} using logarithmic derivative;
  \code{RB} using the ratio between Ricatti-Bessel
  functions.}
}
\value{
The coefficients \eqn{a_n=C_nT_a} and \eqn{b_n=C_nT_b} and
also the coefficients \eqn{C_n}, \eqn{T_a} and \eqn{T_b}.
}
\description{
Calculates the Lorentz-Mie coefficients \eqn{a_n} and
\eqn{b_n}.
}
\details{
Parser to the two ways to calculate the Lorentz-Mie
coefficents.
}
\examples{
m<-5
x<-3
a<-lmie.exp(m,x,by="LD")
b<-lmie.exp(m,x,by="RB")
print(a-b) # Must be close to zero
}

