\name{vswf.bbp}
\alias{vswf.bbp}
\title{Beam Shape Coefficients for the Bessel Beam type P}
\description{Beam Shape Coefficients for Bessel beam type P,
   whose vector potential A have only the p component.}
\usage{
vswf.bbp(gama,kz,xo,yo,zo,lmax,M,s,p)
}
\arguments{
   \item{gama,kz}{Wave number components}
   \item{xo,yo,zo}{Origin of the expansion system}
   \item{lmax}{Maximum number of l terms}
   \item{M}{Mode of the beam}
   \item{s}{Chirality of the beam, s=1 or s=-1}
   \item{p}{Non zero component of the vector potential A}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbp},\link{vwfd.bbp}}
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
M<-4                    # x wavefield mode
N<-1                    # y wavefield mode
S<-1                    # Chirality
P<-1                    # Mode (+/- 1)
# Wave Field Parameters
k=2*pi/lambda           # Propagating wavenumber
kx<-M*pi/a              # x component of the wavevector
ky<-N*pi/b              # y component of the wavevector
gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
kz<-sqrt(k^2-gama^2)    # z component of the wavevector
# Geometry of the calculations
NP=200                  # Number of points in each direction (all equal)
dx<-2*a/(NP-1)          # Sampling in x direction
dy<-2*b/(NP-1)          # Sampling in y direction
x<-seq(-a,a,by=dx)      # x vector of positions
y<-seq(-b,b,by=dy)      # y vector of positions
z<-x                    # z vector of positions
lmax<-20                # Number of partial waves to sum (for l we have 2l+1 m)
#lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y)))))
#-------------------------------------------------------------------------------
# RANDOM POSITION (INDEX)
xo<-sample(x,1)+dx/2
yo<-sample(y,1)+dy/2
zo<-0
#-------------------------------------------------------------------------------
# CHANGE POSITIONS TO THE REFERENCE SYSTEM
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
z<-z-zo
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# TESTING POSITION
#-------------------------------------------------------------------------------
xe<-sample(x,1)
ye<-sample(y,1)
ze<-0
#-------------------------------------------------------------------------------
# WAVEFIELD
Em<-((1+P)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,(M-1),S)/2-
   P*((gama/k)^2)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*(P-1),S)/2)
Ez<--1i*P*S*((gama/k)*(kz/k))*
   vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*P,S)/sqrt(2)
Ep<-((1-P)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,(M-1),S)/2+
   P*((gama/k)^2)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*(P+1),S)/2)
#
Hm<--1i*P*(kz/k)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,M-1,S)*(1+P)/2
Hz<--S*(gama/k)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*P,S)/sqrt(2)
Hp<--1i*P*(kz/k)*vswf.psi(gama,kz,xe+xo,ye+yo,ze+zo,M-1,S)*(1-P)/2
#-------------------------------------------------------------------------------
E<-c(Em,Ez,Ep)
H<-c(Hm,Hz,Hp)
#-------------------------------------------------------------------------------
u<-vswf.bbp(gama,kz,xo,yo,zo,lmax,M,S,P)
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
print(cbind(E.PWE,H.PWE))
cat("PWE EXPANSION\n")
print(cbind(E,H))
\dontrun{
#-------------------------------------------------------------------------------
# LARGE SCALE CALCULATIONS (OVER ALL x,y,z POINTS)
#-------------------------------------------------------------------------------
# Wave Field Definitions
BBP<-vwfd.bbp(P,M,S,gama,kz,x+xo,y+yo,ze+zo)
#
tem.BBP<-array(BBP$Em,c(NP,NP,1))[,,1]
#tez.BBP<-array(BBP$Ez,c(NP,NP,1))[,,1]
#tep.BBP<-array(BBP$Ep,c(NP,NP,1))[,,1]
#thm.BBP<-array(BBP$Hm,c(NP,NP,1))[,,1]
#thz.BBP<-array(BBP$Hz,c(NP,NP,1))[,,1]
#thp.BBP<-array(BBP$Hp,c(NP,NP,1))[,,1]
# Parte Real
x11();image(x+xo,y+yo,Re(tem.BBP),main="BBP Em",col=cm.colors(1024))
abline(v=xo,h=yo)
abline(v=xe+xo,h=ye+yo,col='green')
#x11();image(x+xo,y+yo,Re(tez.BBP),main="BBP Ez",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(tep.BBP),main="BBP Ep",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thm.BBP),main="BBP Hm",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thz.BBP),main="BBP Hz",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thp.BBP),main="BBP Hp",col=cm.colors(1024))
# ParteTmaBnaWa
#x11();image(x+xo,y+yo,Im(tem.BBP),main="BBP Em",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tez.BBP),main="BBP Ez",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tep.BBP),main="BBP Ep",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thm.BBP),main="BBP Hm",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thz.BBP),main="BBP Hz",col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thp.BBP),main="BBP Hp",col=cm.colors(1024))
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
#x11();image(x+xo,y+yo,Re(tez.PWE),main=paste("PWE Ez",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(tep.PWE),main=paste("PWE Ep",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thm.PWE),main=paste("PWE Hm",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thz.PWE),main=paste("PWE Hz",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Re(thp.PWE),main=paste("PWE Hp",lmax),col=cm.colors(1024))
# ParteTmaBnaWa
#x11();image(x+xo,y+yo,Im(tem.PWE),main=paste("PWE Em",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tez.PWE),main=paste("PWE Ez",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(tep.PWE),main=paste("PWE Ep",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thm.PWE),main=paste("PWE Hm",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thz.PWE),main=paste("PWE Hz",lmax),col=cm.colors(1024))
#x11();image(x+xo,y+yo,Im(thp.PWE),main=paste("PWE Hp",lmax),col=cm.colors(1024))
#-------------------------------------------------------------------------------
}
}
