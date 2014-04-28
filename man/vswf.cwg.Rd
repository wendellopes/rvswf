\name{vswf.cwg}
\alias{vswf.cwg}
\title{Beam Shape Coefficient for the Cylindrical Wave Guide}
\description{Calculate the specific term for Cylindrical Wave Guides. It's
   called from vswf.gwg with mode=1(TM) or mode=2(TE).}
\usage{
vswf.cwg(TM=TRUE,gama, kz, x, y, z, lmax, m, s = 1)
}
\arguments{
   \item{TM}{Boolean, indicates the mode. If 'TRUE', TM, else, TE.}
   \item{gama,kz}{Wave number components}
   \item{x,y,z}{Origin of the expansion}
   \item{lmax}{Maximum number of l terms}
   \item{m}{The mode of the Bessel function}
   \item{s}{Chirality of the beam}
}
\seealso{
   \code{\link{vwfd.cwg},\link{vwfd.cwg},\link{vwfd.bbz},\link{vwfd.bbp}}
}
\examples{
#-------------------------------------------------------------------------------
# WAVE GUIDE PARAMETERS
#-------------------------------------------------------------------------------
# Zeros of Bessel Functions
#-------------------------------------------------------------------------------
ZJ <-matrix(c( 2.4048, 3.8317, 5.1356, 6.3802, 7.5883, 8.7715,
               5.5201, 7.0156, 8.4172, 9.7610,11.0647,12.3386,
               8.6537,10.1735,11.6198,13.0152,14.3725,15.7002,
              11.7915,13.3237,14.7960,16.2235,17.6160,18.9801,
              14.9309,16.4706,17.9598,19.4094,20.8269,22.2178),
     ncol=6,byrow=TRUE)
#-------------------------------------------------------------------------------
# Zeros of derivatives of Bessel Functions
#-------------------------------------------------------------------------------
ZdJ<-matrix(c( 3.8317, 1.8412, 3.0542, 4.2012, 5.3175, 6.4156,
               7.0156, 5.3314, 6.7061, 8.0152, 9.2824,10.5199,
              10.1735, 8.5363, 9.9695,11.3459,12.6819,13.9872,
              13.3237,11.7060,13.1704,14.5858,15.9641,17.3128,
              16.4706,14.8636,16.3475,17.7887,19.1960,20.5755),
     ncol=6,byrow=TRUE)
#-------------------------------------------------------------------------------
cat("Nth zeros (lines) of the Mth Bessel functions (column)\n")
print(ZJ)
cat("Nth zeros (lines) of the derivative of Mth Bessel functions (column)\n")
print(ZdJ)
#-------------------------------------------------------------------------------
ux<-seq(0,12,.1)
plot(ux,besselJ(ux,0),type='l')
points(ZJ[1,1],0)
points(ZJ[2,1],0)
points(ZJ[3,1],0)
points(ux,besselJ(ux,1),type='l',col='red')
points(ZJ[1,2],0,col='red')
points(ZJ[2,2],0,col='red')
points(ZJ[3,2],0,col='red')
points(ux,besselJ(ux,2),type='l',col='blue')
points(ZJ[1,3],0,col='blue')
points(ZJ[2,3],0,col='blue')
points(ZJ[3,3],0,col='blue')
grid()
#-------------------------------------------------------------------------------
# Wave Guide Parameters
#-------------------------------------------------------------------------------
# Basic Parameters
lambda=.5e-6            # Propagating wavelength
R<-5*lambda             # Radius of the waveguide
M<-5  # 0<=M<=5         # Order of the (derivative of the) Bessel function
N<-4  # 1<=N<=5         # Order of the zero of the derivative/Bessel function
TM<-FALSE               # Mode: TM (Bessel function), TE (derivative)
#TM<-TRUE               # Mode: TM (Bessel function), TE (derivative)
# Wave field Parameters
k=2*pi/lambda           # Propagating wavenumber
S=-1                    # Chirality of the beam
#-------------------------------------------------------------------------------
# MODES: M <- M+1 because J_0(x) ocupies column 1
# TE - ZdJ (derivative)
# TM - ZJ  (Bessel function)
#-------------------------------------------------------------------------------
if(TM){
   g<-ZJ[N,M+1]/R
   MODE<--1
   MODE.PWE<-4
}else{
   g<-ZdJ[N,M+1]/R
   MODE<-1
   MODE.PWE<-3
}
kz<-sqrt(k^2-g^2)
#-------------------------------------------------------------------------------
# Geometry of the calculations
NP=200
x<-seq(-R,R,by=2*R/(NP-1))
y<-x
z<-x
lmax<-10
#lmax<-as.integer(max(c(10,max(k*R))))
#-------------------------------------------------------------------------------
# ORIGIN OF THE EXPANSION (THE ORIGIN MUST BE INSIDE THE CYLINDER: RADIUS = rho)
#-------------------------------------------------------------------------------
rho<-function(x,y){
   return(sqrt(x^2+y^2))
}
rhoo<-outer(x,y,FUN='rho')
rhoo[rhoo<=R]<--1
rhoo[rhoo>R]<-0
# Choosing a point inside the cylinder
choose.inside.radius<-function(x,y,R){
   xo<-sample(x,1)
   yo<-sample(y,1)
   ro<-rho(xo,yo)
   while(ro>R){
      xo<-sample(x,1)
      yo<-sample(y,1)
      ro<-rho(xo,yo)
   }
   return(c(xo,yo))
}
po<-choose.inside.radius(x,y,R) #Origin of the new system
pe<-choose.inside.radius(x,y,R) #point inside the cylinder
# New system of coordinates (expansion reference frame)
xo<-po[1]
yo<-po[2]
zo<-0
#xo<-yo<-zo<-0    # Origin at the center of the waveguide
# Points at which the expansion will be performed
xe<-pe[1]
ye<-pe[2]
ze<-0
\dontrun{
#x11();image(x,y,rhoo,asp=1,main='Cylinder in the WG Reference frame')
#abline(v=xo,h=yo)
#abline(v=xe,h=ye,col='green')
}
# Changing the reference frame
x<-x-xo
y<-y-yo
z<-z-zo
xe<-xe-xo
ye<-ye-yo
\dontrun{
x11();image(x,y,rhoo,asp=1,main='Cylinder in the Expansion Reference frame')
abline(v=xo-xo,h=yo-yo)
abline(v=xe,h=ye,col='green')
}
#-------------------------------------------------------------------------------
# CALCULATIONS OF THE WAVEFIELDS (EXACT VALUES)
#-------------------------------------------------------------------------------
# MODO TRANSVERSAL MAGNETICO s=1
# psi(gama,kz,x,y,z,m,s=1)
E.TM.H.TE.m<-vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M-S,s=S)*( 1i*S*kz/g)/sqrt(2)
E.TM.H.TE.z<-vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M,s=S)
E.TM.H.TE.p<-vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M+S,s=S)*(-1i*S*kz/g)/sqrt(2)
#
H.TM.E.TE.m<--vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M-S,s=S)*(k/g)*S/sqrt(2)
H.TM.E.TE.z<-0                                              
H.TM.E.TE.p<--vswf.psi(g,kz,xe+xo,ye+yo,ze+zo,M+S,s=S)*(k/g)*S/sqrt(2)
#
E.TM.H.TE<-c(E.TM.H.TE.m,E.TM.H.TE.z,E.TM.H.TE.p)
H.TM.E.TE<-MODE*c(H.TM.E.TE.m,H.TM.E.TE.z,H.TM.E.TE.p)
#-------------------------------------------------------------------------------
if(TM){
   cat("WAVEFIELD DEFINITIONS TM MODE\n")
   print(data.frame(E.TM=E.TM.H.TE,H.TM=H.TM.E.TE))
}else{
   cat("WAVEFIELD DEFINITIONS TE MODE\n")
   print(data.frame(E.TE=H.TM.E.TE,H.TE=E.TM.H.TE))
}
#-------------------------------------------------------------------------------
# BSC CALCULATIONS
#-------------------------------------------------------------------------------
u<-vswf.cwg(TM,g,kz,xo,yo,zo,lmax,m=M,s=S)
v<-vswf.hmp(k,xe,ye,ze,lmax)
#
Em.PWE<-sum(u$GTE*v$M.m-u$GTM*v$N.m)
Ez.PWE<-sum(u$GTE*v$M.z-u$GTM*v$N.z)
Ep.PWE<-sum(u$GTE*v$M.p-u$GTM*v$N.p)

Hm.PWE<-sum(u$GTM*v$M.m+u$GTE*v$N.m)
Hz.PWE<-sum(u$GTM*v$M.z+u$GTE*v$N.z)
Hp.PWE<-sum(u$GTM*v$M.p+u$GTE*v$N.p)
#-------------------------------------------------------------------------------
# COMPARING
#-------------------------------------------------------------------------------
E.PWE<-c(Em.PWE,Ez.PWE,Ep.PWE)
H.PWE<-c(Hm.PWE,Hz.PWE,Hp.PWE)
print(cbind(E.PWE,H.PWE))
\dontrun{
#-------------------------------------------------------------------------------
# LARGE SCALE CALCULATIONS (OVER ALL x,y,z POINTS) 
#-------------------------------------------------------------------------------
CWG<-vwfd.cwg(TE=!TM,M,S,g,kz,x+xo,y+yo,ze+zo)
#
tem.CWG<-array(CWG$Em,c(NP,NP,1))[,,1]
#tez.CWG<-array(CWG$Ez,c(NP,NP,1))[,,1]
#tep.CWG<-array(CWG$Ep,c(NP,NP,1))[,,1]
#thm.CWG<-array(CWG$Hm,c(NP,NP,1))[,,1]
#thz.CWG<-array(CWG$Hz,c(NP,NP,1))[,,1]
#thp.CWG<-array(CWG$Hp,c(NP,NP,1))[,,1]
# Real part
x11();image(x+xo,y+yo,Re(tem.CWG),main="CWG Em",col=cm.colors(1024));grid()
abline(v=xo,h=yo)
abline(v=xe+xo,h=ye+yo,col='green')
#x11();image(x+xo,y+yo,Re(tez.CWG),main="CWG Ez",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(tep.CWG),main="CWG Ep",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(thm.CWG),main="CWG Hm",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(thz.CWG),main="CWG Hz",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(thp.CWG),main="CWG Hp",col=cm.colors(1024));grid()
# Imaginary part
#x11();image(x+xo,y+yo,Im(tem.CWG),main="CWG Em",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(tez.CWG),main="CWG Ez",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(tep.CWG),main="CWG Ep",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(thm.CWG),main="CWG Hm",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(thz.CWG),main="CWG Hz",col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(thp.CWG),main="CWG Hp",col=cm.colors(1024));grid()
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
# Real part
x11();image(x+xo,y+yo,Re(tem.PWE),main=paste("PWE Em",lmax),col=cm.colors(1024));grid()
abline(v=xo,h=yo)
abline(v=xe+xo,h=ye+yo,col='green')
#x11();image(x+xo,y+yo,Re(tez.PWE),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(tep.PWE),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(thm.PWE),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(thz.PWE),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Re(thp.PWE),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid()
# Imaginary part
#x11();image(x+xo,y+yo,Im(tem.PWE),main=paste("PWE Em",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(tez.PWE),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(tep.PWE),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(thm.PWE),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(thz.PWE),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid()
#x11();image(x+xo,y+yo,Im(thp.PWE),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid()
#-------------------------------------------------------------------------------
}
}

