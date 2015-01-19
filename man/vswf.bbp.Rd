\name{vswf.bbp}
\alias{vswf.bbp}
\title{Beam Shape Coefficient for Circularly polarized Bessel Beam}
\usage{
vswf.bbp(gama, kz, xo, yo, zo, lmax, M, s, p)
}
\arguments{
  \item{gama}{The transversal component of the wave vector
  \eqn{\gamma}.}

  \item{kz}{The longitudinal component of the wave vector
  \eqn{k_z}.}

  \item{xo}{The position \eqn{x} at which the expansion is
  performed.}

  \item{yo}{The position \eqn{y} at which the expansion is
  performed.}

  \item{zo}{The position \eqn{z} at which the expansion is
  performed.}

  \item{lmax}{The maximum value of \eqn{l} in the
  expansion.}

  \item{M}{Order of the Bessel Beam.}

  \item{s}{Chirality of the Bessel Beam (\eqn{S=\pm 1}).}

  \item{p}{Polarity of the Bessel beam (\eqn{S=\pm 1}).}
}
\value{
The Beam Shape Coefficients for a Circularly polarized
Bessel Beam.
}
\description{
Beam Shape Coefficient for Circularly polarized Bessel Beam
}
\details{
Calculate the Beam Shape Coefficients for the partial wave
expansion in Vector Spherical Wave Functions.
}
\examples{
lambda=.5e-6            # Propagating wavelength
a<-8*lambda             # Size x of the waveguide (Rectangular Wave Guide)
b<-8*lambda             # Size y of the waveguide (Rectangular Wave Guide)
l<-6*lambda             # Size z of the waveguide (Rectangular Wave Guide)
M<-5                    # x wavefield mode
N<-4                    # y wavefield mode
S<--1                   # Chirality
P<- 1                    # Mode (+/- 1)
# Wave Field Parameters
k=2*pi/lambda           # Propagating wavenumber
kx<-M*pi/a              # x component of the wavevector
ky<-N*pi/b              # y component of the wavevector
gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
kz<-sqrt(k^2-gama^2)    # z component of the wavevector
TM<-ifelse(P==-1,TRUE,FALSE)
#-------------------------------------------------------------------------------
# Geometry of the calculations
#-------------------------------------------------------------------------------
NPX<-4*50+1
NPY<-4*60+1
NPZ<-4*20+1
dx<-2*a/(NPX-1)
dy<-2*b/(NPY-1)
dz<-2*l/(NPZ-1)
x<-seq(-a,a,by=dx)
y<-seq(-b,b,by=dy)
z<-0
#-------------------------------------------------------------------------------
lmax<- 80 # Number of partial waves to sum (for some l we have 2l+1 values of m)
#-------------------------------------------------------------------------------
# POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#-------------------------------------------------------------------------------
# ARBITRARY: We select two points, the is for the origin of the expansion,
#            the second is for the test, that must have be move to the new
#            system of coordinates x'=x-xo
xo<-yo<-zo<-0
x<-sample(x,1)
y<-sample(y,1)
#-------------------------------------------------------------------------------
# BSC CALCULATIONS
#-------------------------------------------------------------------------------
BBP<-vwfd.bbp(P,M,S,gama,kz,x+xo,y+yo,z+zo)  # WAVE FIELD DEFINITION
BSC<-vswf.bbp(gama,kz,xo,yo,zo,lmax,M,S,P)   # BEAM SHAPE COEFFICIENTS
PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)  # PARTIAL WAVE EXPANSION
#-------------------------------------------------------------------------------
# VALUES
#-------------------------------------------------------------------------------
cat("Distance x from origin in wavelength (from ",-a/lambda,"to ",a/lambda,"):",x/lambda,"\\n")
cat("Distance y from origin in wavelength (from ",-b/lambda,"to ",b/lambda,"):",y/lambda,"\\n")
df<-data.frame(
   PWE=c(PWE$Em,PWE$Ez,PWE$Ep,PWE$Hm,PWE$Hz,PWE$Hp),
   BBP=c(BBP$Em,BBP$Ez,BBP$Ep,BBP$Hm,BBP$Hz,BBP$Hp),
   row.names=c("Em","Ez","Ep","Hm","Hz","Hp")
   )
df$DIF<-df$PWE-df$BBP
print(df)
}
\seealso{
\code{\link{vwfd.bbp}}, \code{\link{vswf.bbp}},
\code{\link{vswf.pwd}}.
}

