\name{vwfd.rwg}
\alias{vwfd.rwg}
\title{Vector Wave Field Definition: Rectangular Wave Guide}
\usage{
vwfd.rwg(TE = TRUE, kx, ky, kz, x, y, z)
}
\arguments{
  \item{TE}{Type of the wave field.}

  \item{kx}{Component \eqn{x} of the wave vector(vector).}

  \item{ky}{Component \eqn{y} of the wave vector(vector).}

  \item{kz}{Component \eqn{z} of the wave vector(vector).}

  \item{kx}{Component \eqn{x} of the position vector(single
  value).}

  \item{ky}{Component \eqn{y} of the position vector(single
  value).}

  \item{kz}{Component \eqn{z} of the position vector(single
  value).}
}
\value{
A list with the input values, the number of points and the
complex six values of components of the electromagnetic
fields.
}
\description{
Vector Wave Field Definition: Rectangular Wave Guide
}
\details{
Calculates the vector components of a rectangular wave
guide of type Transversal Electric (TE), that means that
\eqn{E_z=0}, or Transversal Magnetic (TM), that means
\eqn{H_z=0}, in terms of the components of the wave vector
and the position vector.
}

