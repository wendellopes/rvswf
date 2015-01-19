\name{vswf.jlm}
\alias{vswf.jlm}
\title{Convert the value of (l,m) to a single value}
\usage{
vswf.jlm(l, m)
}
\arguments{
  \item{l}{The orbital angular momentum eigenvalue}

  \item{m}{The eigenvalue of the projection of angular
  momentum on \eqn{z} axis}
}
\value{
A unique value used to transform a double summation on a
single one.
}
\description{
Convert the value of (l,m) to a single value
}
\examples{
for(l in 0:5){
   for(m in -l:l){
      print(c(l,m,vswf.jlm(l,m)))
   }
}
}

