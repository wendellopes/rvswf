\name{vwfd.bbz}
\alias{vwfd.bbz}
\title{Vector Wave Field Definition: Transversal Bessel Beam}
\usage{
vwfd.bbz(TM = TRUE, M, S, g, kz, x, y, z)
}
\arguments{
  \item{TE}{Type of the wave field.}

  \item{M}{Order of the Bessel Beam.}

  \item{S}{Chirality of the Bessel Beam (\eqn{S=\pm 1})}

  \item{g}{Transversal component of the Bessel Beam
  (\eqn{\gamma}).}

  \item{kz}{Component \eqn{z} of the wave vector (single
  value).}

  \item{x}{Component \eqn{x} of the position vector
  (vector).}

  \item{y}{Component \eqn{y} of the position vector
  (vector).}

  \item{z}{Component \eqn{z} of the position vector
  (vector).}
}
\value{
A list with the input values, the number of points and the
six complex values of components of the electromagnetic
fields.
}
\description{
Vector Wave Field Definition: Transversal Bessel Beam
}
\details{
Calculates the vector components of Bessel Beam of type
Transversal Electric (TE), that means that \eqn{E_z=0}, or
Transversal Magnetic (TM), that means \eqn{H_z=0}, in terms
of the components of the wave vector and the position
vector.
}

