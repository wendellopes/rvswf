\name{lcfa.cyl}
\alias{lcfa.cyl}
\title{Calculates the logarithmic derivative of the Cylindrical Bessel functions.}
\usage{
lcfa.cyl(nmax, x, code = "C", NMAX = 200)
}
\arguments{
  \item{nmax}{The maximum order of \eqn{J_n(x)}}

  \item{x}{The argument of the functions. Can be complex.}

  \item{code}{If you prefer to use native R or C language.}

  \item{NMAX}{Maximum number of iterations}
}
\value{
An array of the logarithmic derivative of Cylindrical
Bessel functions and the ratio between two consecutive
Cylindrical Bessel functions. from 0 to \code{nmax} at
point \code{x}
}
\description{
Calculates the logarithmic derivative of the Cylindrical
Bessel functions.
}
\details{
\code{lcfa.cyl} calculates the logarithmic derivative of
Cylindrical Bessel functions
\eqn{D_n(x)=J_n'(x)/J_n(x)}{D[n]=J_n'/J_n} using downward
recurrence. The system of equations is given by
\eqn{S_n(x)=n/x},
\eqn{\gamma_n=J_n(x)/J_{n+1}(x)}{g[n]=J_n/J_{n+1}} and
\eqn{D_n=J_n'(x)/J_n(x)}. The system can be solved by means
of the recurrence relations of the Cylindrical Bessel
functions
\deqn{\gamma_{n-1}+\frac{1}{\gamma_n}=2S_{n}}{g[n-1]+1/g[n]=2
S[n]}
\deqn{\gamma_{n-1}-\frac{1}{\gamma_n}=2D_{n}}{g[n-1]-1/g[n]=2
D[n]} that can be rewriten \deqn{\gamma_{n}=S_{n+1}+D_{n+1}
}{ g[n]=S[n+1]+D[n+1]} \deqn{\frac{1}{\gamma_n}=S_n-D_n.
}{1/g[n]=S[n ]-D[n ].} The logarithmic derivatives obeys
the relation, \deqn{(S_n-D_n)(S_{n+1}+D_{n+1})=1.
}{(S[n]-D[n])(S[n+1]+D[n+1])=1.} The values are calculated
by downward recurrence, and the inicial values calculated
by Lentz method.
}
\examples{
nmax<-10
x<-5
u.c<-lcfa.cyl(nmax,x,code="C")
u.r<-lcfa.cyl(nmax,x,code="R")
u<-data.frame(
   # Logarithmic Derivatives
   C.LogDev=u.c$Dn,R.LogDev=u.r$Dn,
   # Ratio between Cylindrical Bessel functions
   C.CylRat=u.c$gn,R.CylRat=u.r$gn)
print(u)
}
\seealso{
\code{\link{lcfe.cbi}}, \code{\link{lcfe.cbl}},
\code{\link{lcfe.cbd}}.
}

