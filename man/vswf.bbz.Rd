\name{vswf.bbz}
\alias{vswf.bbz}
\title{Beam Shape Coefficients for the Bessel Beam Z}
\description{Beam Shape Coefficients for the Bessel Beam Z,
   whose vector potential A have only the z component.}
\usage{
vswf.bbz(gama,kz,xo,yo,zo,lmax,M,s)
}
\arguments{
   \item{gama,kz}{Wave number components}
   \item{xo,yo,zo}{Origin of the expansion system}
   \item{lmax}{Maximum number of l terms}
   \item{M}{Mode of the beam}
   \item{s}{Chirality of the beam}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
#-------------------------------------------------------------------------------
# BEAM PARAMETERS
#-------------------------------------------------------------------------------
# We take as base the Rectangular Wave Guide parameters.
# Basic Parameters
lambda=.5e-6            # Propagating wavelength
a<-5*lambda             # Size x of the waveguide (Rectangular Wave Guide)
b<-8*lambda             # Size y of the waveguide (Rectangular Wave Guide)
M<-3                    # x wavefield mode
N<-5                    # y wavefield mode
S<--1                   # Chirality
TM<-TRUE                # Mode (TM-TRUE/TE-FALSE)
# Wave Field Parameters
k=2*pi/lambda           # Propagating wavenumber
kx<-M*pi/a              # x component of the wavevector
ky<-N*pi/b              # y component of the wavevector
gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
kz<-sqrt(k^2-gama^2)    # z component of the wavevector
# Geometry of the calculations
NP<-200                 # Number of points in each direction (all equal)
dx=2*a/(NP-1)           # Sampling in x direction
dy=2*b/(NP-1)           # Sampling in y direction
x<-seq(-a,a,by=dx)      # x vector of positions
y<-seq(-b,b,by=dy)      # y vector of positions
z<-x                    # z vector of positions
lmax<-20                # Number of partial waves to sum (for l we have 2l+1 m)
#lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y)))))
#-------------------------------------------------------------------------------
# RANDOM POSITION (INDEX)
xo<-sample(x,1)         # one can choose arbitrary number, example 42
yo<-sample(y,1)
#zo<-sample(z,1)
zo<-0
#-------------------------------------------------------------------------------
# CHANGE POSITIONS TO THE REFERENCE SYSTEM
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
z<-z-zo
#-------------------------------------------------------------------------------
# TESTING POSITION
#-------------------------------------------------------------------------------
xe<-sample(x,1)
ye<-sample(y,1)
#ze<-sample(z,1)
ze<-0
#-------------------------------------------------------------------------------
# MODE: TE
# MODE: TM -> E <- H, H <- -E
#-------------------------------------------------------------------------------
Em<- 1i*S*((kz/k)*(gama/k))*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,M-S,S)/sqrt(2)
Ez<-((gama/k)^2)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,M,S)
Ep<--1i*S*((kz/k)*(gama/k))*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,M+S,S)/sqrt(2)
#
Hm<-S*(gama/k)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,M-S,S)/sqrt(2)
Hz<-0
Hp<-S*(gama/k)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,M+S,S)/sqrt(2)
#-------------------------------------------------------------------------------
if(!TM){
   E<-c(Em,Ez,Ep)
   H<-c(Hm,Hz,Hp)
}else{
   E<- c(Hm,Hz,Hp)
   H<--c(Em,Ez,Ep)
}
#-------------------------------------------------------------------------------
u<-vswf.bbz(gama,kz,xo,yo,zo,lmax,M,S)
v<-vswf.hmp(k,xe,ye,ze,lmax)

Em.PWE<-sum(u$GTE*v$M.m-u$GTM*v$N.m)
Ez.PWE<-sum(u$GTE*v$M.z-u$GTM*v$N.z)
Ep.PWE<-sum(u$GTE*v$M.p-u$GTM*v$N.p)

Hm.PWE<-sum(u$GTM*v$M.m+u$GTE*v$N.m)
Hz.PWE<-sum(u$GTM*v$M.z+u$GTE*v$N.z)
Hp.PWE<-sum(u$GTM*v$M.p+u$GTE*v$N.p)

E.PWE<-c(Em.PWE,Ez.PWE,Ep.PWE)
H.PWE<-c(Hm.PWE,Hz.PWE,Hp.PWE)
#-------------------------------------------------------------------------------
# COMPARING 
#-------------------------------------------------------------------------------
cat("WAVE FIELD DEFINITION\n")
print(cbind(E,H))
cat("PWE EXPANSION\n")
print(cbind(E.PWE,H.PWE))
\dontrun{
#-------------------------------------------------------------------------------
# LARGE SCALE CALCULATIONS (OVER ALL x,y,z POINTS)
#-------------------------------------------------------------------------------
# Wave Field Definitions
BBZ<-vwfd.bbz(TE=TRUE,M,S,gama,kz,x+xo,y+yo,ze+zo)
#
tem.BBZ<-array(BBZ$Em,c(NP,NP,1))[,,1]
#tez.BBZ<-array(BBZ$Ez,c(NP,NP,1))[,,1]
#tep.BBZ<-array(BBZ$Ep,c(NP,NP,1))[,,1]
#thm.BBZ<-array(BBZ$Hm,c(NP,NP,1))[,,1]
#thz.BBZ<-array(BBZ$Hz,c(NP,NP,1))[,,1]
#thp.BBZ<-array(BBZ$Hp,c(NP,NP,1))[,,1]
#
x11();image(x+xo,y+yo,Re(tem.BBZ),main="BBZ Em",col=cm.colors(1024))
abline(v=xo,h=yo)
abline(v=xe+xo,h=ye+yo,col='green')
#x11();image(x+xo,y+yo,Re(tez.BBZ),main="BBZ Ez",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(tep.BBZ),main="BBZ Ep",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thm.BBZ),main="BBZ Hm",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thz.BBZ),main="BBZ Hz",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thp.BBZ),main="BBZ Hp",col=cm.colors(1024))
#
#x11();image(x+xo,y+yo,Im(tem.BBZ),main="BBZ Em",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tez.BBZ),main="BBZ Ez",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tep.BBZ),main="BBZ Ep",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thm.BBZ),main="BBZ Hm",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thz.BBZ),main="BBZ Hz",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thp.BBZ),main="BBZ Hp",col=cm.colors(1024))
#-------------------------------------------------------------------------------
# Partial Wave Expansion
PWE<-vswf.pwe(k,x,y,ze,lmax,u$GTE,u$GTM)
#
tem.PWE<-array(PWE$Em,c(NP,NP,1))[,,1]
#tez.PWE<-array(PWE$Ez,c(NP,NP,1))[,,1]
#tep.PWE<-array(PWE$Ep,c(NP,NP,1))[,,1]
#thm.PWE<-array(PWE$Hm,c(NP,NP,1))[,,1]
#thz.PWE<-array(PWE$Hz,c(NP,NP,1))[,,1]
#thp.PWE<-array(PWE$Hp,c(NP,NP,1))[,,1]
# Parte Real
x11();image(x+xo,y+yo,Re(tem.PWE),main=paste("PWE Em",lmax),col=cm.colors(1024))
abline(v=xo,h=yo)
abline(v=xe+xo,h=ye+yo,col='green')
#x11();image(x+xo,y+yo,Re(tez.PWE),main=paste("PWE Ez",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(tep.PWE),main=paste("PWE Ep",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thm.PWE),main=paste("PWE Hm",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thz.PWE),main=paste("PWE Hz",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thp.PWE),main=paste("PWE Hp",lmax),col=cm.colors(1024))
# Parte Imaginaria
#x11();image(x+xo,y+yo,Im(tem.PWE),main=paste("PWE Em",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tez.PWE),main=paste("PWE Ez",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tep.PWE),main=paste("PWE Ep",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thm.PWE),main=paste("PWE Hm",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thz.PWE),main=paste("PWE Hz",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thp.PWE),main=paste("PWE Hp",lmax),col=cm.colors(1024))
#-------------------------------------------------------------------------------
}
}

