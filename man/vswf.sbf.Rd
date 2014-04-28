\name{vswf.sbf}
\alias{vswf.sbf}
\title{Spherical Bessel Functions}
\description{Spherical Bessel Functions are calculated vi downward recurrence
   and Lentz method and its derivative}
\usage{
vswf.sbf(NN,xo,compare=FALSE,verbose=FALSE)
}
\arguments{
  \item{NN}{Maximum order to be calculated}
  \item{xo}{Point at which it's calculated}
  \item{compare}{Compare with bessel_jl_steed_array(lmax=NN,x=xo)}
  \item{verbose}{Shows renormalizations.}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
lmax=5
xo<-3
u<-vswf.sbf(lmax,xo,compare=TRUE,verbose=TRUE)
print(u)
}
