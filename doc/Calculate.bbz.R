#-------------------------------------------------------------------------------
# WAVE GUIDE PARAMETERS
#-------------------------------------------------------------------------------
rm(list=ls())
#-------------------------------------------------------------------------------
# BEAM PARAMETERS
#-------------------------------------------------------------------------------
# We take as base the Rectangular Wave Guide parameters.
# Basic Parameters
lambda=.5e-6            # Propagating wavelength
a<-8*lambda             # Size x of the waveguide (Rectangular Wave Guide)
b<-8*lambda             # Size y of the waveguide (Rectangular Wave Guide)
l<-6*lambda             # Size z of the waveguide (Rectangular Wave Guide)
M<-3                    # x wavefield mode
N<-5                    # y wavefield mode
S<--1                   # Chirality
TM<-FALSE               # Mode (TM-TRUE/TE-FALSE)
# Wave Field Parameters
k=2*pi/lambda           # Propagating wavenumber
kx<-M*pi/a              # x component of the wavevector
ky<-N*pi/b              # y component of the wavevector
gama<-sqrt(kx^2+ky^2)   # gama component of the wavevector
kz<-sqrt(k^2-gama^2)    # z component of the wavevector
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
z<-seq(-l,l,by=dz)
#-------------------------------------------------------------------------------
lmax<-2 # Number of partial waves to sum (for some l we have 2l+1 values of m)
#-------------------------------------------------------------------------------
# POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#-------------------------------------------------------------------------------
# ARBITRARY
set.seed(512)
xo<-sample(x,1)
yo<-sample(y,1)
zo<-sample(z,1)
# FIXED
xo<-x[NPX%/%4+1]
yo<-y[NPY%/%4+1]
zo<-0
xo<-yo<-zo<-0
#-------------------------------------------------------------------------------
# CHANGE THE REFERENCE SYSTEM TO THE NEW POSITIONS
#-------------------------------------------------------------------------------
x<-x-(xo+dx/2)
y<-y-(yo+dy/2)
z<-z-(zo+dz/2)
z<-0;NPZ<-1
#-------------------------------------------------------------------------------
# BSC CALCULATIONS
#-------------------------------------------------------------------------------
BBZ<-vwfd.bbz(TM,M,S,gama,kz,x+xo,y+yo,z+zo)
BSC<-vswf.bbz(TM,gama,kz,xo,yo,zo,lmax,M,S)
PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)
if(TM){ # TM implies Hz=0
   tez.BBZ<-array(BBZ$Ez,c(NPZ,NPY,NPX))[1,,]
   tez.PWE<-array(PWE$Ez,c(NPZ,NPY,NPX))[1,,]
}else{  # TE implies Ez=0
   thz.BBZ<-array(BBZ$Hz,c(NPZ,NPY,NPX))[1,,]
   thz.PWE<-array(PWE$Hz,c(NPZ,NPY,NPX))[1,,]
}
# #-------------------------------------------------------------------------------
# # IMAGE
# #-------------------------------------------------------------------------------
source("plots.bbz.R")
#source("PoyntingVectorBBZ.R")
# source("Image.R")
# u<-seq(0,2*pi,length.out=200)
# if(TM){
#    zl<-range(Re(tez.BBZ))
#    x11()
#    Image((y+yo)/lambda,(x+xo)/lambda,z=Re(tez.BBZ),nlevels=256,axes=TRUE,
#       color.palette=cm.colors,asp=1,#zlim=zl,
#       plot.axes={
#          axis(1);
#          axis(2);
#          abline(h=xo/lambda,v=yo/lambda,col='green')
#       },
#       xlab=expression(y/lambda),ylab=expression(x/lambda))
#    x11()
#    Image((y+yo)/lambda,(x+xo)/lambda,z=Re(tez.PWE),nlevels=256,axes=TRUE,
#       color.palette=cm.colors,asp=1,#zlim=zl,
#       plot.axes={
#          axis(1);
#          axis(2);
#          abline(h=xo/lambda,v=yo/lambda,col='green')
#       },
#       xlab=expression(y/lambda),ylab=expression(x/lambda))
#}else{
#    zl<-range(Re(thz.BBZ))
#    x11()
#    Image((y+yo)/lambda,(x+xo)/lambda,z=Re(thz.BBZ),nlevels=256,axes=TRUE,
#       color.palette=cm.colors,asp=1,#zlim=zl,
#       plot.axes={
#          axis(1);
#          axis(2);
#          abline(h=xo/lambda,v=yo/lambda,col='green')
#          },
#       xlab=expression(y/lambda),ylab=expression(x/lambda))
#    x11()
#    Image((y+yo)/lambda,(x+xo)/lambda,z=Re(thz.PWE),nlevels=256,axes=TRUE,
#    color.palette=cm.colors,asp=1,#zlim=zl,
#       plot.axes={
#          axis(1);
#          axis(2);
#          abline(h=xo/lambda,v=yo/lambda,col='green')
#       },
#       xlab=expression(y/lambda),ylab=expression(x/lambda))
# }
# 
