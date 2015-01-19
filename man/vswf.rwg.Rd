\name{vswf.rwg}
\alias{vswf.rwg}
\title{Beam Shape Coefficients for Rectangular Wave Guides.}
\usage{
vswf.rwg(TM = TRUE, kx, ky, kz, x, y, z, lmax)
}
\arguments{
  \item{TM}{Type of the wave field.}

  \item{kx}{Component \eqn{x} of the wave vector (single
  value).}

  \item{ky}{Component \eqn{y} of the wave vector (single
  value).}

  \item{kz}{Component \eqn{z} of the wave vector (single
  value).}

  \item{x}{Component \eqn{x} of the origin of the expansion
  (vector).}

  \item{y}{Component \eqn{y} of the origin of the expansion
  (vector).}

  \item{z}{Component \eqn{z} of the origin of the expansion
  (vector).}
}
\value{
The Beam Shape Coefficients \eqn{G^{TE}_{lm}} and
\eqn{G^{TM}_{lm}}.
}
\description{
Beam Shape Coefficients for Rectangular Wave Guides.
}
\details{
Calculates the Beam Shape Coefficients used to do the
Partial Wave Expansion on Vector Spherical Wave Functions.
}
\examples{
lambda<-.5e-6           # Propagating wavelength
a<-7*lambda             # Size x of the waveguide
b<-5*lambda             # Size y of the waveguide
M<-6                    # x wavefield mode
N<-5                    # y wavefield mode
#-------------------------------------------------------------------------------
# Wave Field Parameters
#-------------------------------------------------------------------------------
k<-2*pi/lambda          # Propagating wavenumber
kx<-M*pi/a              # x component of the wavevector
ky<-N*pi/b              # y component of the wavevector
gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
kz<-sqrt(k^2-gama^2)    # z component of the wavevector
#-------------------------------------------------------------------------------
# Geometry of the calculations
#-------------------------------------------------------------------------------
NPX=200                  # Number of points in each direction (all equal)
NPY=200                  # Number of points in each direction (all equal)
#-------------------------------------------------------------------------------
# Vectors
#-------------------------------------------------------------------------------
dx<-a/(NPX-1)
dy<-b/(NPY-1)
x<-seq(0,a,dx)          # x vector of positions
y<-seq(0,b,dy)          # y vector of positions
z<-0
#-------------------------------------------------------------------------------
TM<-FALSE
lmax<- 40
#-------------------------------------------------------------------------------
# POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#-------------------------------------------------------------------------------
# ARBITRARY
xo<-a/2
yo<-b/2
zo<-0
#-------------------------------------------------------------------------------
# CHANGE THE REFERENCE SYSTEM TO THE NEW POSITIONS
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
#-------------------------------------------------------------------------------
# ARBITRARY POINT FOR CALCULATIONS
#-------------------------------------------------------------------------------
x<-sample(x,1)
y<-sample(y,1)
#-------------------------------------------------------------------------------
# BSC CALCULATIONS
#-------------------------------------------------------------------------------
RWG<-vwfd.rwg(TE=!TM,kx,ky,kz,x+xo,y+yo,z+zo)
BSC<-vswf.rwg(TM,kx,ky,kz,xo,yo,zo,lmax)
PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)
#-------------------------------------------------------------------------------
# VALUES
#-------------------------------------------------------------------------------
cat("Distance x from origin in wavelength (from ",-a/(2*lambda),"to ",a/(2*lambda),"):",x/lambda,"\\n")
cat("Distance y from origin in wavelength (from ",-b/(2*lambda),"to ",b/(2*lambda),"):",y/lambda,"\\n")
df<-data.frame(
   PWE=c(PWE$Em,PWE$Ez,PWE$Ep,PWE$Hm,PWE$Hz,PWE$Hp),
   RWG=c(RWG$Em,RWG$Ez,RWG$Ep,RWG$Hm,RWG$Hz,RWG$Hp),
   row.names=c("Em","Ez","Ep","Hm","Hz","Hp")
   )
df$DIF<-df$PWE-df$RWG
print(df)
}
\seealso{
\code{\link{vwfd.rwg}}, \code{\link{vswf.pwe}},
\code{\link{vswf.gwg}}.
}

