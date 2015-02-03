\name{lmie.sct}
\alias{lmie.sct}
\title{Scattering coeficients in the Lorentz-Mie}
\usage{
lmie.sct(m, x, dist = FALSE)
}
\arguments{
  \item{m}{Ratio of refraction indices.}

  \item{x}{Form factor}

  \item{dist}{If is to distinguish contributions of
  \eqn{a_n} and \eqn{b_n}.}
}
\value{
Extinction and Scattering coefficients \eqn{Qsca} and
\eqn{Qext} if \code{dist=FALSE}. If \code{dist=TRUE}, so
the contribution of \eqn{a_n} and \eqn{b_n} are
distinguished in \eqn{Q_{sca}\to Q_{sca}+Q_{scb}} and
\eqn{Q_{ext}\to Q_{exa}+Q_{exb}}.
}
\description{
Scattering coeficients in the Lorentz-Mie
}
\details{
Scattering and Extinction coefficients.
}
\examples{
dx<-.1   # FAST VALUE
#dx<-.001 # HIGH RESOLUTION
m<-1.33
x<-seq(dx,50,dx)
u<-lmie.sct(m,x)
v<-lmie.sct(m,x,dist=TRUE)
plot(x,u$Qsca,type='l')
points(x,v$Qsca,type='l',col='red')
points(x,v$Qscb,type='l',col='blue')
}

