% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/reff.sjn.r
\name{reff.sjn}
\alias{reff.sjn}
\title{Reference function for Spherical Bessel Functions.}
\usage{
reff.sjn(n, x)
}
\arguments{
\item{n}{The order of the function}

\item{x}{The argument of the function}
}
\value{
The value of the spherical Bessel function \eqn{j_n(x)}
}
\description{
Reference function for Spherical Bessel Functions.
}
\details{
It uses de built in Cylindrical Bessel function to calculate the
Spherical Bessel function, taking care about the zero.
}

