\name{ridx.h2o}
\alias{ridx.h2o}
\title{Parameters: Refraction index of water}
\usage{
ridx.h2o(l)
}
\arguments{
  \item{l}{The wavelength in micrometers}
}
\value{
The index of refraction of the water as function of the
wavelength
}
\description{
Parameters: Refraction index of water
}
\examples{
l<-seq(.38,.75,.01)
u<-ridx.h2o(l)
plot(l,u,type='l')
}
\seealso{
\code{\link{ridx.h2o}}
}

