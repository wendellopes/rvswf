#-------------------------------------------------------------------------------
# WAVE GUIDE PARAMETERS
#-------------------------------------------------------------------------------
rm(list=ls())
#-------------------------------------------------------------------------------
# Basic Parameters
#-------------------------------------------------------------------------------
lambda<-.5e-6           # Propagating wavelength
a<-7*lambda             # Size x of the waveguide
b<-5*lambda             # Size y of the waveguide
l<-3*lambda
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
NPZ=2                   # Number of points in each direction (all equal)
#-------------------------------------------------------------------------------
# Vectors
#-------------------------------------------------------------------------------
dx<-a/(NPX-1)
dy<-b/(NPY-1)
dz<-l/(NPZ-1)
x<-seq(0,a,dx)          # x vector of positions
y<-seq(0,b,dy)          # y vector of positions
z<-seq(0,l,dz)          # z vector of positions
#-------------------------------------------------------------------------------
TM<-FALSE
lmax<- 4
#-------------------------------------------------------------------------------
# POSITION AT WHICH THE EXPANSION WILL BE PERFORMED  (REFERENCE SYSTEM)
#-------------------------------------------------------------------------------
# ARBITRARY
set.seed(512)
xo<-sample(x,1)
yo<-sample(y,1)
zo<-sample(z,1)
# FIXED
xo<-x[NPX%/%2+1]
yo<-y[NPY%/%2+1]
zo<-0
#-------------------------------------------------------------------------------
# CHANGE THE REFERENCE SYSTEM TO THE NEW POSITIONS
#-------------------------------------------------------------------------------
x<-x-(xo+dx/2)
y<-y-(yo+dy/2)
z<-z-(zo+dz/2)
z<-0;NPZ<-1
# #-------------------------------------------------------------------------------
# # BSC CALCULATIONS
# #-------------------------------------------------------------------------------
RWG<-vwfd.rwg(TE=!TM,kx,ky,kz,x+xo,y+yo,z+zo)
BSC<-vswf.rwg(TM,kx,ky,kz,xo,yo,zo,lmax)
PWE<-vswf.pwe(k,x,y,z,lmax,BSC$GTE,BSC$GTM)
if(TM){ # TM implies Hz=0
   tez.RWG<-array(RWG$Ez,c(NPZ,NPY,NPX))[1,,]
   tez.PWE<-array(PWE$Ez,c(NPZ,NPY,NPX))[1,,]
}else{  # TE implies Ez=0
   thz.RWG<-array(RWG$Hz,c(NPZ,NPY,NPX))[1,,]
   thz.PWE<-array(PWE$Hz,c(NPZ,NPY,NPX))[1,,]
}
#-------------------------------------------------------------------------------
# NAMES
#-------------------------------------------------------------------------------
nm.vwf<-"rwg.vwfd.tm.00"
md<-ifelse(TM,"tm","te")
nm.pwe<-ifelse(lmax<10, paste("rwg.vswf.",md,".0",lmax,sep=""),
                        paste("rwg.vswf.",md,"." ,lmax,sep=""))
nm.vwf.i<-paste(nm.vwf,".png",sep="")
nm.vwf.d<-paste(nm.vwf,".Rdata",sep="")
nm.pwe.i<-paste(nm.pwe,".png",sep="")
nm.pwe.d<-paste(nm.pwe,".Rdata",sep="")
#-------------------------------------------------------------------------------
# IMAGE
#-------------------------------------------------------------------------------
source("plots.rwg.R")
#-------------------------------------------------------------------------------
#source("Image.R")
#if(TM){
#   zl<-range(Re(tez.RWG))
#   #1
#   if(!file.exists(nm.vwf.i)){
#      png(nm.vwf.i)
#      Image((y+yo)/lambda,(x+xo)/lambda,z=Re(tez.RWG),nlevels=256,axes=TRUE,color.palette=cm.colors,#zlim=zl,
#         plot.axes={axis(1);axis(2);abline(h=xo/lambda,v=yo/lambda,col='green')},
#         xlab=expression(y/lambda),ylab=expression(x/lambda))
#      dev.off()
#   }
#   #2
#   png(nm.pwe.i)
#   Image((y+yo)/lambda,(x+xo)/lambda,z=Re(tez.PWE),nlevels=256,axes=TRUE,color.palette=cm.colors,#zlim=zl,
#      plot.axes={axis(1);axis(2);abline(h=xo/lambda,v=yo/lambda,col='green')},
#      xlab=expression(y/lambda),ylab=expression(x/lambda))
#   dev.off()
#}else{
#   zl<-range(Re(thz.RWG))
#   #1
#   if(!file.exists(nm.vwf.i)){
#      png(nm.vwf.i)
#      Image((y+yo)/lambda,(x+xo)/lambda,z=Re(thz.RWG),nlevels=256,axes=TRUE,color.palette=cm.colors,#zlim=zl,
#         plot.axes={axis(1);axis(2);abline(h=xo/lambda,v=yo/lambda,col='green')},
#         xlab=expression(y/lambda),ylab=expression(x/lambda))
#      dev.off()
#   }
#   #2
#   png(nm.pwe.i)
#   Image((y+yo)/lambda,(x+xo)/lambda,z=Re(thz.PWE),nlevels=256,axes=TRUE,color.palette=cm.colors,#zlim=zl,
#      plot.axes={axis(1);axis(2);abline(h=xo/lambda,v=yo/lambda,col='green')},
#      xlab=expression(y/lambda),ylab=expression(x/lambda))
#   dev.off()
#}
##-------------------------------------------------------------------------------
## DATASETS
##-------------------------------------------------------------------------------
#if(TM){
#   save(RWG, file=nm.vwf.d)
#   save(PWE, file=nm.pwe.d)
#}else{
#   save(RWG, file=nm.vwf.d)
#   save(PWE, file=nm.pwe.d)
#}
