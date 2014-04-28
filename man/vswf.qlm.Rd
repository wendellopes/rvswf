\name{vswf.qlm}
\alias{vswf.qlm}
\title{Normalized Associated Legendre Polynomials Qlm & Qlm'}
\description{Gives the Normalized Legendre Polinomials Qlm
   which allows one to calculate the Scalar Spherical
   Harmonics. The result is a data.frame containing
   Qlm(cth) and its derivative Qlm'(cth) from 0 to lmax.}
\usage{
vswf.qlm(cth,lmax)
}
\arguments{
  \item{cth}{Argument to the Associated Legendre Polynomials -1<=cth<=1}
  \item{lmax}{Maximum value of l. The result goes from 0 to lmax}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
lmax=4                    # Maximum number of l
th<-seq(-pi,pi,by=pi/300) # Generic list of values theta
tho<-sample(th,1)         # Aleatory value of theta
u<-vswf.qlm(cos(tho),lmax)
print(u)
# To obtain an element of order l,m, one must do
(l<-sample(0:lmax,1))
(m<-sample(-l:l,1))
u[vswf.jlm(l,m),]
}
