\name{vswf.psi}
\alias{vswf.psi}
\title{Psi_m(x) is a recurrent function in cylindrical symetries}
\description{Used in Cylindrical Wave Guides and in Bessel Beams, is 
a solution of the wave equation in cylindrical coordinates}
\usage{
vswf.psi(gama, kz, x, y, z, m, s = 1)
}

\arguments{
   \item{gama,kz}{Components of the wave vector}
   \item{x,y,z}{Position}
   \item{m}{Order of the Cylindrical Bessel Function}
   \item{s}{Chirality}
}
\seealso{
   \code{
      \link{vswf.rwg},
      \link{vswf.cwg},
      \link{vwfd.rwg},
      \link{vwfd.cwg},
      \link{vwfd.bbz},
      \link{vwfd.bbp}
      }
}
\examples{
#-------------------------------------------------------------------------------
x<-seq(-10,10,by=.1)
y<-x
# random location
xo<-sample(x,1)
yo<-sample(y,1)
# fixed location
#xo=5
#yo=-10
# Changing the reference systema
x<-x-xo
y<-y-yo
#-------------------------------------------------------------------------------
# using outer product to do it faster
M<-3
psi<-function(x,y){
   return(vswf.psi(1,0,x,y,0,M,1))
}
U<-outer(x+xo,y+yo,FUN='psi')
#-------------------------------------------------------------------------------
# Creating the somation
jM<-20
U.MJ<-0*U
for(j in -jM:jM){
   Psi.o<-vswf.psi(1,0,xo,yo,0,M-j,1)
   psi<-function(x,y){
      return(vswf.psi(1,0,x,y,0,j,1))
   }
   U.j<-outer(x,y,FUN='psi')
   U.MJ<-U.MJ+Psi.o*U.j
}
# diference between two results
# one may perform calculations over a wide range of jM, depending on how far
#    the calculating point is from the origin.
DU<-U.MJ-U
#-------------------------------------------------------------------------------
\dontrun{
x11();image(x+xo,y+yo,Re(U),main="psi_M(x+xo)")
abline(v=xo,h=yo)
x11();image(x+xo,y+yo,Re(U.MJ),main="Sum_j psi_{M-j}(xo)*psi_j(x)")
abline(v=xo,h=yo)
x11();image(x+xo,y+yo,Re(DU),main="psi_M(x+xo)-Sum_j psi_{M-j}(xo)*psi_j(x)")
abline(v=xo,h=yo)
}
#-------------------------------------------------------------------------------
}
