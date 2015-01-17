#' Vector Wave Field Definition: Rectangular Wave Guide
#' 
#' @details Calculates the vector components of a rectangular wave guide of type 
#' Transversal Electric (TE), that means that \eqn{E_z=0}, or Transversal
#' Magnetic (TM), that means \eqn{H_z=0},  in terms of the components of the
#' wave vector and the position vector.
#' @param TE Type of the wave field.
#' @param kx Component \eqn{x} of the wave vector(vector).
#' @param ky Component \eqn{y} of the wave vector(vector).
#' @param kz Component \eqn{z} of the wave vector(vector).
#' @param kx Component \eqn{x} of the position vector(single value).
#' @param ky Component \eqn{y} of the position vector(single value).
#' @param kz Component \eqn{z} of the position vector(single value).
#' @return A list with the input values, the number of points and the complex
#' six values of components of the electromagnetic fields.
vwfd.rwg<-function(TE=TRUE,kx,ky,kz,x,y,z){
if(TE)
{
   te<-1
}else{
   te<-0
}
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("vwfd_rwg",
   TE=as.integer(te),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   kx=as.double(kx) ,ky=as.double(ky) ,kz=as.double(kz),
   x=as.double(x) ,y=as.double(y) ,z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
