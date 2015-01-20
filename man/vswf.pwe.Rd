\name{vswf.pwe}
\alias{vswf.pwe}
\title{Partial Wave Expansion in Vector Spherical Harmonics.}
\usage{
vswf.pwe(k, x, y, z, lmax, gte, gtm)
}
\arguments{
  \item{k}{The module of the wave vector.}

  \item{x}{The component \eqn{x} of the position vector.}

  \item{y}{The component \eqn{y} of the position vector.}

  \item{z}{The component \eqn{z} of the position vector.}

  \item{lmax}{The maximum value of \eqn{l}.}

  \item{gte}{One of the components of the expansion.}

  \item{gtm}{The oder component of the expansion.}
}
\value{
A list with the input values, the position and the values
of the electromagnetic fields.
}
\description{
Partial Wave Expansion in Vector Spherical Harmonics.
}
\details{
Calculate the electromagnetic fields \eqn{\bm{E}} and
\eqn{Z\bm{H}} as functions of \eqn{G^{TE}_{lm}} and
\eqn{G^{TM}_{lm}} for a given position and from 0 to lmax.
}
\examples{
lambda=.5e-6
k=2*pi/lambda
lmax<-20
a=lambda
NP=200
z<-seq(-a,a,2*a/(NP-1))
xo<-sample(z,1)
yo<-sample(z,1)
zo<-sample(z,1)
Em<-exp(1i*k*zo)
Ez<-0
Ep<-0
Hm<--1i*Em
Hz<--1i*Ez
Hp<--1i*Ep
v<-vswf.hmp(k,xo,yo,zo,lmax)
u<-vswf.mpw(lmax,norm=TRUE)

Em.pweZ<-sum(u$GTE*v$M.m-u$GTM*v$N.m)
Ez.pweZ<-sum(u$GTE*v$M.z-u$GTM*v$N.z)
Ep.pweZ<-sum(u$GTE*v$M.p-u$GTM*v$N.p)

Hm.pweZ<-sum(u$GTM*v$M.m+u$GTE*v$N.m)
Hz.pweZ<-sum(u$GTM*v$M.z+u$GTE*v$N.z)
Hp.pweZ<-sum(u$GTM*v$M.p+u$GTE*v$N.p)
#-------------------------------------------------------------------------------
# PLANE WAVE
#-------------------------------------------------------------------------------
w<-vswf.gpw(0,0,k,1/sqrt(2),1i/sqrt(2),0,lmax)

Em.pweK<-sum(w$GTE*v$M.m-w$GTM*v$N.m)
Ez.pweK<-sum(w$GTE*v$M.z-w$GTM*v$N.z)
Ep.pweK<-sum(w$GTE*v$M.p-w$GTM*v$N.p)

Hm.pweK<-sum(w$GTM*v$M.m+w$GTE*v$N.m)
Hz.pweK<-sum(w$GTM*v$M.z+w$GTE*v$N.z)
Hp.pweK<-sum(w$GTM*v$M.p+w$GTE*v$N.p)
#-------------------------------------------------------------------------------
# CONFERINDO
#-------------------------------------------------------------------------------
b<-data.frame(DEF=c(Em,Ez,Ep,Hm,Hz,Hp),
              row.names=c("Em","Ez","Ep","Hm","Hz","Hp"))
b$MPW<-c(Em.pweZ,Ez.pweZ,Ep.pweZ,Hm.pweZ,Hz.pweZ,Hp.pweZ)
b$GPW<-c(Em.pweK,Ez.pweK,Ep.pweK,Hm.pweK,Hz.pweK,Hp.pweK)
print(b)
}

