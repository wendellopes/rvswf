\name{vwfd.bbp}
\alias{vwfd.bbp}
\title{Vector Wave Field Definition: Circularly polarized Bessel Beam}
\usage{
vwfd.bbp(P, M, S, g, kz, x, y, z)
}
\arguments{
  \item{P}{Polarity of the Bessel beam (\eqn{S=\pm 1}).}

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
Vector Wave Field Definition: Circularly polarized Bessel
Beam
}
\details{
Calculates the vector components of Bessel Beam circularly
polarized.
}

