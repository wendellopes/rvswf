\name{vswf.rwg}
\alias{vswf.rwg}
\title{Beam Shape Coefficient for the Rectangular Wave Guide}
\description{Calculate the specific term for Rectangular Wave Guides. It's
   called from vswf.gwg with mode=1(TM) or mode=2(TE).}
\usage{
vswf.rwg(TM = TRUE, kx, ky, kz, x, y, z, lmax)
}
\arguments{
   \item{TM}{Boolean, indicates the mode. If 'TRUE', TM, else, TE.}
   \item{kx,ky,kz}{Wave number components}
   \item{x,y,z}{Origin of the expansion}
   \item{lmax}{Maximum number of l terms}
}
\seealso{
   \code{\link{vwfd.rwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
#-------------------------------------------------------------------------------
# WAVE GUIDE PARAMETERS
#-------------------------------------------------------------------------------
# Basic Parameters
lambda<-.5e-6           # Propagating wavelength
a<-8*lambda             # Size x of the waveguide
b<-5*lambda             # Size y of the waveguide
M<-7                    # x wavefield mode
N<-6                    # y wavefield mode
TM<-TRUE
# Wave Field Parameters
k<-2*pi/lambda          # Propagating wavenumber
kx<-M*pi/a              # x component of the wavevector
ky<-N*pi/b              # y component of the wavevector
gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
kz<-sqrt(k^2-gama^2)    # z component of the wavevector
# Geometry of the calculations
NP=200                  # Number of points in each direction (all equal)
dx<-a/(NP-1)            # Sampling in x direction
dy<-b/(NP-1)            # Sampling in y direction
x<-seq(0,a,dx)          # x vector of positions
y<-seq(0,b,dy)          # y vector of positions
z<-x                    # z vector of positions
lmax<-20                # Number of partial waves to sum (for l we have 2l+1 m)
#lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y))))) 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#-------------------------------------------------------------------------------
# RANDOM POSITION (INDEX)
nxo<-sample(1:NP,1)     # one can choose arbitrary number, example 42  
nyo<-sample(1:NP,1)
nzo<-sample(1:NP,1)
# RANDOM POSITION (VALUE)
xo<-x[nxo]
yo<-y[nyo]
zo<-z[nzo]
# avoid zeros in the position vectors (BUG)
xo<-xo+dx/2
yo<-yo+dy/2
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
nx<-sample(1:NP,1)
ny<-sample(1:NP,1)
nz<-sample(1:NP,1)
xe<-x[nx]
ye<-y[ny]
ze<-0
#xe<-x[73]
#ye<-y[80]
#ze<-z[173]
#-------------------------------------------------------------------------------
# CALCULATIONS OF THE WAVEFIELDS (EXACT VALUES)
#-------------------------------------------------------------------------------
# TRANSVERSAL MAGNETIC MODE
E.TM.x<-1i*((kz*kx)/gama^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TM.y<-1i*((kz*ky)/gama^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TM.z<-sin(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
#
H.TM.x<--1i*((k*ky)/gama^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TM.y<- 1i*((k*kx)/gama^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TM.z<- 0
# reference vectorial system (e-,eo,e+)
E.TM.m<-(E.TM.x-1i*E.TM.y)/sqrt(2)
E.TM.p<-(E.TM.x+1i*E.TM.y)/sqrt(2)
E.TM<-c(E.TM.m,E.TM.z,E.TM.p)
#
H.TM.m<-(H.TM.x-1i*H.TM.y)/sqrt(2)
H.TM.p<-(H.TM.x+1i*H.TM.y)/sqrt(2)
H.TM<-c(H.TM.m,H.TM.z,H.TM.p)
#-------------------------------------------------------------------------------
# TRANSVERSAL ELETRIC MODE 
H.TE.x<--1i*((kz*kx)/gama^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TE.y<--1i*((kz*ky)/gama^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TE.z<- cos(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
#
E.TE.x<--1i*((k*ky)/gama^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TE.y<- 1i*((k*kx)/gama^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TE.z<- 0
# reference vectorial system (e-,eo,e+)
E.TE.m<-(E.TE.x-1i*E.TE.y)/sqrt(2)
E.TE.p<-(E.TE.x+1i*E.TE.y)/sqrt(2)
E.TE<-c(E.TE.m,E.TE.z,E.TE.p)
#
H.TE.m<-(H.TE.x-1i*H.TE.y)/sqrt(2)
H.TE.p<-(H.TE.x+1i*H.TE.y)/sqrt(2)
H.TE<-c(H.TE.m,H.TE.z,H.TE.p)
#-------------------------------------------------------------------------------
# BSC CALCULATIONS
#-------------------------------------------------------------------------------
#MODE = 1, TM
#MODE = 2, TE
u<-vswf.rwg(TM,kx,ky,kz,xo,yo,zo,lmax)
v<-vswf.hmp(k,xe,ye,ze,lmax)
#
Em.PWE<-sum(u$GTE*v$M.m-u$GTM*v$N.m)
Ez.PWE<-sum(u$GTE*v$M.z-u$GTM*v$N.z)
Ep.PWE<-sum(u$GTE*v$M.p-u$GTM*v$N.p)
#
Hm.PWE<-sum(u$GTM*v$M.m+u$GTE*v$N.m)
Hz.PWE<-sum(u$GTM*v$M.z+u$GTE*v$N.z)
Hp.PWE<-sum(u$GTM*v$M.p+u$GTE*v$N.p)
#
E.PWE<-c(Em.PWE,Ez.PWE,Ep.PWE)
H.PWE<-c(Hm.PWE,Hz.PWE,Hp.PWE)
#-------------------------------------------------------------------------------
# COMPARING
#-------------------------------------------------------------------------------
cat("PWE EXPANSION\n")
print(cbind(E.PWE,H.PWE))
if(TM){
   cat("TM MODE\n")
   print(cbind(E.TM,H.TM))
}else{
   cat("TE MODE\n")
   print(cbind(E.TE,H.TE))
}
\dontrun{
#-------------------------------------------------------------------------------
# LARGE SCALE CALCULATIONS (OVER ALL x,y,z POINTS)
#-------------------------------------------------------------------------------
RWG<-vwfd.rwg(TE=TRUE,kx,ky,kz,x+xo,y+yo,ze+zo)
#
tem.RWG<-array(RWG$Em,c(NP,NP,1))[,,1]
#tez.RWG<-array(RWG$Ez,c(NP,NP,1))[,,1]
#tep.RWG<-array(RWG$Ep,c(NP,NP,1))[,,1]
#thm.RWG<-array(RWG$Hm,c(NP,NP,1))[,,1]
#thz.RWG<-array(RWG$Hz,c(NP,NP,1))[,,1]
#thp.RWG<-array(RWG$Hp,c(NP,NP,1))[,,1]
# Parte Real
X11();image(x+xo,y+yo,Re(tem.RWG),main="RWG Re Em",col=cm.colors(1024));grid()
abline(v=xo,h=yo)
#X11();image(x+xo,y+yo,Re(tez.RWG),main="RWG Re Ez",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(tep.RWG),main="RWG Re Ep",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(thm.RWG),main="RWG Re Hm",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(thz.RWG),main="RWG Re Hz",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(thp.RWG),main="RWG Re Hp",col=cm.colors(1024));grid()
# Parte Imaginaria
#X11();image(x+xo,y+yo,Im(tem.RWG),main="RWG Im Em",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(tez.RWG),main="RWG Im Ez",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(tep.RWG),main="RWG Im Ep",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(thm.RWG),main="RWG Im Hm",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(thz.RWG),main="RWG Im Hz",col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(thp.RWG),main="RWG Im Hp",col=cm.colors(1024));grid()
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
X11();image(x+xo,y+yo,Re(tem.PWE),main=paste("PWE Re Em",lmax),col=cm.colors(1024));grid()
abline(v=xo,h=yo)
#X11();image(x+xo,y+yo,Re(tez.PWE),main=paste("PWE Re Ez",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(tep.PWE),main=paste("PWE Re Ep",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(thm.PWE),main=paste("PWE Re Hm",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(thz.PWE),main=paste("PWE Re Hz",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Re(thp.PWE),main=paste("PWE Re Hp",lmax),col=cm.colors(1024));grid()
# Parte Imaginaria
#X11();image(x+xo,y+yo,Im(tem.PWE),main=paste("PWE Im Em",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(tez.PWE),main=paste("PWE Im Ez",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(tep.PWE),main=paste("PWE Im Ep",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(thm.PWE),main=paste("PWE Im Hm",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(thz.PWE),main=paste("PWE Im Hz",lmax),col=cm.colors(1024));grid()
#X11();image(x+xo,y+yo,Im(thp.PWE),main=paste("PWE Im Hp",lmax),col=cm.colors(1024));grid()
#-------------------------------------------------------------------------------
}
}
